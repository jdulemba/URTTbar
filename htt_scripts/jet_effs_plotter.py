'''
Jet_Effs Analyzer Plotter macro
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os, glob
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
from rootpy.plotting import views, Graph, Hist
from argparse import ArgumentParser
import URAnalysis.Utilities.roottools as urtools
import ROOT
from URAnalysis.Utilities.tables import latex_table
import argparse

parser = argparse.ArgumentParser(description='Create plots using files from jet_effs.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full_Analysis).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, or Combined).')
parser.add_argument('plot', help='Choose type of plots to generate (Gen_Plots, Gen_Ratios, Matched_Objects, Delta_Plots, BTagging_Eff, Merged_Perms).')
args = parser.parse_args()

if args.analysis == "Test":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../%s.jet_effs.test.root' % args.sample, 'read')
		normfile = views.NormalizeView(root_open('../%s.jet_effs.test.root' % args.sample, 'read'))#normalized file
		plotter = BasePlotter(
			'plots/jet_effs/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Sample_ttJetsM0 = root_open('../ttJetsM0.jet_effs.test.root', 'read')
	#	Sample_ttJetsM700 = root_open('../ttJetsM700.jet_effs.test.root', 'read')
	#	Sample_ttJetsM1000 = root_open('../ttJetsM1000.jet_effs.test.root', 'read')
		plotter = BasePlotter(
			'plots/jet_effs/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)

elif args.analysis == "Full_Analysis":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../results/%s/jet_effs/%s.root' % (jobid, args.sample), 'read')
		normfile = views.NormalizeView(root_open('../results/%s/jet_effs/%s.root' % (jobid, args.sample), 'read'))
		plotter = BasePlotter(
			'plots/jet_effs/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Full_ttJetsM0 = root_open('../results/%s/jet_effs/ttJets.root' % jobid, 'read')
	#	Full_ttJetsM700 = root_open('../results/%s/jet_effs/ttJetsM700.root' % jobid, 'read')
	#	Full_ttJetsM1000 = root_open('../results/%s/jet_effs/ttJetsM1000.root' % jobid, 'read')
		plotter = BasePlotter(
			'plots/jet_effs/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)


##### Comparing right/wrong combinations from TTBarSolver settings ####
##open files
#NS_MD_file = root_open('../ttJets.jet_effs_NS_MassDisc.test.root', 'read') #only NS and MassDis right used in solver
#NS_MD_ratio_file = root_open('../ttJets.jet_effs_NS_MassDisc_Ratio.test.root', 'read') #NS and MassDis right/wrong used in solver
#NS_MD_Ang_file = root_open('../ttJets.jet_effs_NS_MassDisc_Ang.test.root', 'read') #NS, MassDisc, Ang vars right used in solver
#NS_MD_Ang_ratio_file = root_open('../ttJets.jet_effs_NS_MassDisc_Ang_Ratio.test.root', 'read') #NS, MassDisc, Ang vars right/wrong used in solver
#
#files = [
#	(NS_MD_file, 'NS_MD', 'blue'),
#	(NS_MD_ratio_file, 'NS_MD_ratio', 'red'),
#	(NS_MD_Ang_file, 'NS_MD_Ang', 'cyan'),
#	(NS_MD_Ang_ratio_file, 'NS_MD_Ang_ratio', 'black')
#]
##files = [HR_LW_file]
#
#to_draw = []
#for infile, name, colors in files:
#	right_hist = asrootpy(infile.Get('njets_RIGHT')).Clone()
#	THAD = asrootpy(infile.Get('njets_RIGHT_THAD')).Clone()
#	TLEP = asrootpy(infile.Get('njets_RIGHT_TLEP')).Clone()
#	WRONG = asrootpy(infile.Get('njets_WRONG')).Clone()
#	OTHER = asrootpy(infile.Get('njets_OTHER')).Clone()
#	Wrong_hists = THAD+TLEP+WRONG+OTHER
#
#	Eff_hist = right_hist/(right_hist+Wrong_hists)
#	Eff_hist.GetXaxis().SetRangeUser(0, 10)
#	plotter.set_histo_style(Eff_hist, title=name, color=colors)
#	to_draw.append(Eff_hist)
#	plotter.plot(Eff_hist, xtitle='nJets', ytitle='Efficiency')
#	plotter.save('%s_Eff_Ang' % name)
#
#plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle='nJets', ytitle='Efficiency')
#plotter.save('Eff_Ang')
#####

##### Global Var. Definitions ####
defcol = 'black' # default color
defyax = 'a.u.' # yaxis title
if args.sample == "ttJetsM0":
	mass_min = 0
	mass_max = 2000
if args.sample == "ttJetsM700":
	mass_min = 700
	mass_max = 1000
if args.sample == "ttJetsM1000":
	mass_min = 1000
	mass_max = 2000

DRvals = {
		('DRP4', '#Delta R < 0.4', 'red'),
		('DRP5', '#Delta R < 0.5', 'blue'),
		('DRP6', '#Delta R < 0.6', 'green'),
		('DRP8', '#Delta R < 0.8', 'black')
	}

#lpos = 'NE' #legend position
#types = ['dtype', 'utype'] #variable is down- or up-type jet
#colors = ['red', 'blue'] #plot legend colors (down is red, up is blue)
#labels = ['down', 'up'] #legend variable labels
#####

########################################################################################
if args.plot == "Gen_Plots":
	# Gen Object Plots
	Gen_DR = [
		('DR_LepBHad', '#Delta R (l, b_{h})'),
		('DR_LepBLep', '#Delta R (l, b_{l})'),
		('DR_LepWJa', '#Delta R (l, WJa)'),
		('DR_LepWJb', '#Delta R (l, WJb)'),
		('DR_BHadBLep', '#Delta R (b_{h}, b_{l})'),
		('DR_BHadWJa', '#Delta R (b_{h}, WJa)'),
		('DR_BHadWJb', '#Delta R (b_{h} WJb)'),
		('DR_BLepWJa', '#Delta R (b_{l}, WJa)'),
		('DR_BLepWJb', '#Delta R (b_{l}, WJb)'),
		('DR_WJaWJb', '#Delta R (WJa, WJb)'),
		('DRmin_thad', 'min t_{h} #Delta R'),
		('DRmin_tlep', 'min t_{l} #Delta R'),
		('DRmax_thad', 'max t_{h} #Delta R')
	]
	
	for var, xaxis in Gen_DR:
	      hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	      plotter.save(var)
	
	Gen_Pt = [
		('Pt_Lep', 'l p_{t}'),
		('Pt_BHad', 'b_{h} p_{t}'),
		('Pt_BLep', 'b_{l} p_{t}'),
		('Pt_WJa', 'WJa p_{t}'),
		('Pt_WJb', 'WJb p_{t}'),
		('Pt_ttbar', 't#bar t p_{t}'),
		('Pt_thad', 't_{h} p_{t}'),
		('Pt_tlep', 't_{l} p_{t}')
	]
	
	for var, xaxis in Gen_Pt:
	      hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	      plotter.save(var)
	
	Gen_Eta = [
		('Eta_Lep', 'l #eta'),
		('Eta_BHad', 'b_{h} #eta'),
		('Eta_BLep', 'b_{l} #eta'),
		('Eta_WJa', 'WJa #eta'),
		('Eta_WJb', 'WJb #eta'),
	]
	
	for var, xaxis in Gen_Eta:
	      hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	      plotter.save(var)
	
	Gen_DR_2D = [
		('DR_LepBHad_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (l, b_{h})'),
		('DR_LepBLep_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (l, b_{l})'),
		('DR_LepWJa_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (l, WJa)'),
		('DR_LepWJb_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (l, WJb)'),
		('DR_BHadBLep_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (b_{h}, b_{l})'),
		('DR_BHadWJa_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (b_{h}, WJa)'),
		('DR_BHadWJb_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (b_{h}, WJb)'),
		('DR_BLepWJa_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (b_{l}, WJa)'),
		('DR_BLepWJb_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (b_{l}, WJb)'),
		('DR_WJaWJb_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (WJa, WJb)'),
		('DRmin_thad_vs_mttbar', 'm_{t#bar t} [GeV]', 'min t_{h} #Delta R'),
		('DRmin_tlep_vs_mttbar', 'm_{t#bar t} [GeV]', 'min t_{l} #Delta R'),
		('DRmax_thad_vs_mttbar', 'm_{t#bar t} [GeV]', 'max t_{h} #Delta R'),
		('DRmin_thad_vs_ptthad', 't_{h} p_{t}', 'min t_{h} #Delta R'),
		('DRmin_tlep_vs_pttlep', 't_{l} p_{t}', 'min t_{l} #Delta R'),
		('DRmax_thad_vs_ptthad', 't_{h} p_{t}', 'max t_{h} #Delta R')
	]
	
	for var, xaxis,yaxis in Gen_DR_2D:
		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	#	plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=yaxis, drawstyle='hist')
		hist.Draw('colz')
		hist.xaxis.set_title(xaxis)
		hist.yaxis.set_title(yaxis)
		if var == 'DRmin_thad_vs_ptthad_hist' or var == 'DRmin_tlep_vs_pttlep_hist' or var == 'DRmax_thad_vs_ptthad_hist':
			hist.xaxis.range_user = 0, 2000.
		else:
			hist.xaxis.range_user = mass_min, mass_max
		plotter.save(var)
	
	Gen_Pt_2D = [
		('Pt_Lep_vs_Mtt', 'm_{t#bar t} [GeV]', 'l p_{t}'),
		('Pt_BHad_vs_Mtt', 'm_{t#bar t} [GeV]', 'b_{h} p_{t}'),
		('Pt_BLep_vs_Mtt', 'm_{t#bar t} [GeV]', 'b_{l} p_{t}'),
		('Pt_WJa_vs_Mtt', 'm_{t#bar t} [GeV]', 'WJa p_{t}'),
		('Pt_WJb_vs_Mtt', 'm_{t#bar t} [GeV]', 'WJb p_{t}')
	]
	
	for var, xaxis,yaxis in Gen_Pt_2D:
		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	#	plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
		hist.Draw('colz')
		hist.xaxis.set_title(xaxis)
		hist.yaxis.set_title(yaxis)
		plotter.save(var)
	
	Gen_Eta_2D = [
		('Eta_Lep_vs_Mtt', 'm_{t#bar t} [GeV]', 'l #eta'),
		('Eta_BHad_vs_Mtt', 'm_{t#bar t} [GeV]', 'b_{h} #eta'),
		('Eta_BLep_vs_Mtt', 'm_{t#bar t} [GeV]', 'b_{l} #eta'),
		('Eta_WJa_vs_Mtt', 'm_{t#bar t} [GeV]', 'WJa #eta'),
		('Eta_WJb_vs_Mtt', 'm_{t#bar t} [GeV]', 'WJb #eta')
	]
	
	for var, xaxis, yaxis in Gen_Eta_2D:
		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	#	plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
		hist.Draw('colz')
		hist.xaxis.set_title(xaxis)
		hist.yaxis.set_title(yaxis)
		plotter.save(var)

	MassTTbar = [
		('Mass_ttbar', 'm_{t#bar t} [GeV]'),
		('Mass_thad', 'm_{t_{h}} [GeV]'),
		('Mass_tlep', 'm_{t_{l}} [GeV]')
	]
	
	for var, xaxis in MassTTbar:
		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	#	hist.xaxis.range_user = 150, 700
		mean = hist.GetMean()
		rms = hist.GetRMS()
		plotter.set_histo_style(hist, color=defcol)
	#	print("RMS = %f" % RMS)
		plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
		plotter.save(var)

	nJets = [
		('nJets')
	]

	for var in nJets:
		efficiencies = []
		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
		mean = hist.GetMean()
		rms = hist.GetRMS()
		total= hist.Integral()
		plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=var, ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
		plotter.save(var)
		for i in range(0, hist.GetXaxis().GetNbins()):
			efficiencies.append(format(hist.GetBinContent(i+1)/total, '.4f'))
		print(var, efficiencies)

##################################################################################################
if args.plot == "Gen_Ratios":
	DRmin_ratios = [
		('DRmin_thad_lDRP4_vs_mttbar', 'DRmin_thad_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', 't_{h} min #Delta R', 'DRmin_thad_lp4_ratio_vs_mttbar'),
		('DRmin_thad_lDRP4_vs_ptthad', 'DRmin_thad_gDRP4_vs_ptthad', 't_{h} p_{t}', '#Delta R < 0.4 Fraction', 't_{h} min #Delta R', 'DRmin_thad_lp4_ratio_vs_ptthad'),
		('DRmin_tlep_lDRP4_vs_mttbar', 'DRmin_tlep_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', 't_{l} min #Delta R', 'DRmin_tlep_lp4_ratio_vs_mttbar'),
		('DRmin_tlep_lDRP4_vs_pttlep', 'DRmin_tlep_gDRP4_vs_pttlep', 't_{l} p_{t}', '#Delta R < 0.4 Fraction', 't_{l} min #Delta R', 'DRmin_tlep_lp4_ratio_vs_pttlep'),
		
		('DR_LepBHad_lDRP4_vs_mttbar', 'DR_LepBHad_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', '#Delta R (l, b_{h})', 'DR_LepBHad_lp4_ratio_vs_mttbar'),
                ('DR_LepBLep_lDRP4_vs_mttbar', 'DR_LepBLep_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', '#Delta R (l, b_{l})', 'DR_LepBLep_lp4_ratio_vs_mttbar'),
                ('DR_LepWJa_lDRP4_vs_mttbar', 'DR_LepWJa_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', '#Delta R (l, WJa)', 'DR_LepWJa_lp4_ratio_vs_mttbar'),
                ('DR_LepWJb_lDRP4_vs_mttbar', 'DR_LepWJb_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', '#Delta R (l, WJb)', 'DR_LepWJb_lp4_ratio_vs_mttbar'),
                ('DR_BHadBLep_lDRP4_vs_mttbar', 'DR_BHadBLep_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', '#Delta R (b_{h}, b_{l})', 'DR_BHadBLep_lp4_ratio_vs_mttbar'),
                ('DR_BHadWJa_lDRP4_vs_mttbar', 'DR_BHadWJa_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', '#Delta R (b_{h}, WJa)', 'DR_BHadWJa_lp4_ratio_vs_mttbar'),
                ('DR_BHadWJb_lDRP4_vs_mttbar', 'DR_BHadWJb_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', '#Delta R (b_{h}, WJb)', 'DR_BHadWJb_lp4_ratio_vs_mttbar'),
                ('DR_BLepWJa_lDRP4_vs_mttbar', 'DR_BLepWJa_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', '#Delta R (b_{l}, WJa)', 'DR_BLepWJa_lp4_ratio_vs_mttbar'),
                ('DR_BLepWJb_lDRP4_vs_mttbar', 'DR_BLepWJb_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', '#Delta R (b_{l}, WJb)', 'DR_BLepWJb_lp4_ratio_vs_mttbar'),
                ('DR_WJaWJb_lDRP4_vs_mttbar', 'DR_WJaWJb_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', '#Delta R (WJa, WJb)', 'DR_WJaWJb_lp4_ratio_vs_mttbar'),
		('DRmax_thad_lDRP4_vs_mttbar', 'DRmax_thad_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', 't_{h} max #Delta R', 'DRmax_thad_lp4_ratio_vs_mttbar'),
		('DRmax_thad_lDRP4_vs_ptthad', 'DRmax_thad_gDRP4_vs_ptthad', 't_{h} p_{t}', '#Delta R < 0.4 Fraction', 't_{h} max #Delta R', 'DRmax_thad_lp4_ratio_vs_ptthad')
	]

	for lp4, gp4, xaxis, yaxis, title, name in DRmin_ratios:# < 0.4, > 0.4, xaxis, yaxis, png name
		efficiencies = []
		
		LP4 = asrootpy(myfile.Get('Gen_Plots/'+lp4)).Clone()
		GP4 = asrootpy(myfile.Get('Gen_Plots/'+gp4)).Clone()
		DR_Ratio = LP4/(LP4+GP4)
		plotter.set_histo_style(DR_Ratio, color=defcol, markerstyle=20)
		plotter.plot(DR_Ratio)
		DR_Ratio.Draw("P")
		DR_Ratio.xaxis.set_title(xaxis)
		DR_Ratio.yaxis.set_title(yaxis)
#		DR_Ratio.yaxis.range_user = 0, 1
		box = plotter.make_text_box(title, position='NE')
		plotter.save(name)
		
		for i in range(0, DR_Ratio.GetXaxis().GetNbins()):
			efficiencies.append(format(DR_Ratio.GetBinContent(i+1), '.4f'))
		#	mass_range.append(Mu_Eff.GetXaxis().GetBinLowEdge(i+1))
		print(yaxis, efficiencies)

	MassTTbar = [
		('Mass_ttbar', 'm_{t#bar t} [GeV]'),
		('Mass_thad', 'm_{t_{h}} [GeV]'),
		('Mass_tlep', 'm_{t_{l}} [GeV]')
	]
	
	for var, xaxis, in MassTTbar:
		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
		mean = hist.GetMean()
		rms = hist.GetRMS()
		plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
		plotter.save(var)

##############################################################################################
if args.plot == "Matched_Objects":

#	if args.sample == "ttJetsM0":
#		Matched_WJets = [
#			('Matched_perm_WJa_ttbarM700_pt', 'Matched WJa p_{t}', '700 <= m_{t#bar t} < 1000 [GeV]', 'Matched_WJa_ttbarM700_pt'),
#			('Matched_perm_WJa_ttbarM1000_pt', 'Matched WJa p_{t}', 'm_{t#bar t} >= 1000 [GeV]', 'Matched_WJa_ttbarM1000_pt'),
#			('Matched_perm_WJb_ttbarM700_pt', 'Matched WJb p_{t}', '700 <= m_{t#bar t} < 1000 [GeV]', 'Matched_WJb_ttbarM700_pt'),
#			('Matched_perm_WJb_ttbarM1000_pt', 'Matched WJb p_{t}', 'm_{t#bar t} >= 1000 [GeV]', 'Matched_WJb_ttbarM1000_pt')
#		]
#
#		for var, xaxis, titles, names in Matched_WJets:
#			hist = asrootpy(myfile.Get(var)).Clone()
#			mean = hist.GetMean()
#			rms = hist.GetRMS()
#			plotter.set_histo_style(hist, color=defcol)
#			plotter.plot(hist, xtitle=xaxis, ytitle=defyax)
#			box = plotter.make_text_box(titles+'\nMean = %f\nRMS = %f' % (mean,rms), position='NE')
#			plotter.save(names)
#
#		Matched_WJets_Frac = [
#			('Matched_perm_WJa_ttbarM700_frac_p', 'Matched WJa Fractional p', '700 <= m_{t#bar t} < 1000 [GeV]', 'Matched_WJa_ttbarM700_frac_p'),
#			('Matched_perm_WJa_ttbarM1000_frac_p', 'Matched WJa Fractional p', 'm_{t#bar t} >= 1000 [GeV]', 'Matched_WJa_ttbarM1000_frac_p'),
#			('Matched_perm_WJb_ttbarM700_frac_p', 'Matched WJb Fractional p', '700 <= m_{t#bar t} < 1000 [GeV]', 'Matched_WJb_ttbarM700_frac_p'),
#			('Matched_perm_WJb_ttbarM1000_frac_p', 'Matched WJb Fractional p', 'm_{t#bar t} >= 1000 [GeV]', 'Matched_WJb_ttbarM1000_frac_p')
#		]
#
#		for var, xaxis, titles, names in Matched_WJets_Frac:
#			hist = asrootpy(myfile.Get(var)).Clone()
#			mean = hist.GetMean()
#			rms = hist.GetRMS()
##			hist.xaxis.range_user = 0, 1.0
#			plotter.set_histo_style(hist, color=defcol)
#			plotter.plot(hist, xtitle=xaxis, ytitle=defyax)
#			box = plotter.make_text_box(titles+'\nMean = %f\nRMS = %f' % (mean,rms), position='NE')
#			plotter.save(names)

	Matched_Pt = [
		('Matched_perm_BHad_pt', 'Matched b_{h} p_{t}'),
		('Matched_perm_BLep_pt', 'Matched b_{l} p_{t}'),
		('Matched_perm_WJa_pt', 'Matched WJa p_{t}'),
		('Matched_perm_WJb_pt', 'Matched WJb p_{t}')
	]

	for var, xaxis in Matched_Pt:
		to_draw = []
		for folder, legend, colors in DRvals:
			hist = asrootpy(myfile.Get(folder+'/'+var)).Clone()
			mean = hist.GetMean()
			rms = hist.GetRMS()
	#		hist.xaxis.range_user = 25, 2000
			plotter.set_histo_style(hist, title=legend, color=colors)
			to_draw.append(hist)

		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=None, y_range=None, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
#		box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
		plotter.save(var)

	Matched_Eta = [
		('Matched_perm_BHad_eta', 'Matched b_{h} #eta'),
		('Matched_perm_BLep_eta', 'Matched b_{l} #eta'),
		('Matched_perm_WJa_eta', 'Matched WJa #eta'),
		('Matched_perm_WJb_eta', 'Matched WJb #eta')
	]
	
	for var, xaxis in Matched_Eta:
		to_draw = []
		for folder, legend, colors in DRvals:
			hist = asrootpy(myfile.Get(folder+'/'+var)).Clone()
			mean = hist.GetMean()
			rms = hist.GetRMS()
			hist.xaxis.range_user = -2.4, 2.4
			plotter.set_histo_style(hist, title=legend, color=colors)
			to_draw.append(hist)

		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle=xaxis, ytitle=defyax, drawstyle='hist')
#		box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
		plotter.save(var)

	Same_Matched_Objs = [
		('Matched_BHadWJa_ptthad_3J', 'BHad & WJa Objs Matched to Same Jet'),
		('Matched_BHadWJb_ptthad_3J', 'BHad & WJb Objs Matched to Same Jet'),
		('Matched_WJaWJb_ptthad_3J', 'WJa & WJb Objs Matched to Same Jet'),
		('Matched_BHadWJa_ptthad_4J', 'BHad & WJa Objs Matched to Same Jet'),
		('Matched_BHadWJb_ptthad_4J', 'BHad & WJb Objs Matched to Same Jet'),
		('Matched_WJaWJb_ptthad_4J', 'WJa & WJb Objs Matched to Same Jet'),
		('Matched_BHadWJa_ptthad_5PJ', 'BHad & WJa Objs Matched to Same Jet'),
		('Matched_BHadWJb_ptthad_5PJ', 'BHad & WJb Objs Matched to Same Jet'),
		('Matched_WJaWJb_ptthad_5PJ', 'WJa & WJb Objs Matched to Same Jet')
	]
	xaxis = 't_{h} p_{t}'

	for var, titles in Same_Matched_Objs:
		to_draw = []
		for folder, legend, colors in DRvals:
			hist = asrootpy(myfile.Get(folder+'/'+var)).Clone()
			mean = hist.GetMean()
			rms = hist.GetRMS()
			total= hist.Integral()
			#hist.xaxis.range_user = -2.4, 2.4
			plotter.set_histo_style(hist, title=legend, color=colors)
			to_draw.append(hist)

		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle=xaxis, ytitle=defyax, drawstyle='hist')
		#box = plotter.make_text_box(titles+'\nMean = %f\nRMS = %f' % (mean,rms), position='NE')
		plotter.save('Same_'+var)

	Matched_Objs_Ratios = [
		('Matched_BHadWJa_ptthad_3J', 'All_Matched_BHadWJa_ptthad_3J', 'b_{h} & WJa Objects'),
		('Matched_BHadWJb_ptthad_3J', 'All_Matched_BHadWJb_ptthad_3J', 'b_{h} & WJb Objects'),
		('Matched_WJaWJb_ptthad_3J', 'All_Matched_WJaWJb_ptthad_3J', 'WJa & WJb Objects'),
		('Matched_BHadWJa_ptthad_4J', 'All_Matched_BHadWJa_ptthad_4J', 'b_{h} & WJa Objects'),
		('Matched_BHadWJb_ptthad_4J', 'All_Matched_BHadWJb_ptthad_4J', 'b_{h} & WJb Objects'),
		('Matched_WJaWJb_ptthad_4J', 'All_Matched_WJaWJb_ptthad_4J', 'WJa & WJb Objects'),
		('Matched_BHadWJa_ptthad_5PJ', 'All_Matched_BHadWJa_ptthad_5PJ', 'b_{h} & WJa Objects'),
		('Matched_BHadWJb_ptthad_5PJ', 'All_Matched_BHadWJb_ptthad_5PJ', 'b_{h} & WJb Objects'),
		('Matched_WJaWJb_ptthad_5PJ', 'All_Matched_WJaWJb_ptthad_5PJ', 'WJa & WJb Objects')
	]

	yaxis = 'Same Jet Matching Ratio'

	for Same_Match, All_Match, titles in Matched_Objs_Ratios:
		to_draw = []
		for folder, legend, colors in DRvals:
			efficiencies = []
			SameMatch = asrootpy(myfile.Get(folder+'/'+Same_Match)).Clone()
			AllMatch = asrootpy(myfile.Get(folder+'/'+All_Match)).Clone()
			MatchRatio = SameMatch/AllMatch
			plotter.set_histo_style(MatchRatio, title=legend, color=colors, markerstyle=20)
			MatchRatio.Draw("P")
#			MatchRatio.xaxis.set_title(xaxis)
#			MatchRatio.yaxis.set_title(yaxis)
			MatchRatio.yaxis.range_user = 0, 1
			to_draw.append(MatchRatio)

#			for i in range(0, MatchRatio.GetXaxis().GetNbins()):
#				efficiencies.append(format(MatchRatio.GetBinContent(i+1), '.4f'))
#			#	mass_range.append(Mu_Eff.GetXaxis().GetBinLowEdge(i+1))
#			#print(xaxis, efficiencies)

		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), y_range=(0, 1), xtitle=xaxis, ytitle=yaxis)
		box = plotter.make_text_box(titles, position='SE')
		plotter.save('Ratio_'+Same_Match)

#
##	nJets = [
##		('nJets_hist', 'nJets', 'nJets'),
##		('nMatched_objects_3J_hist', 'Matched Objects for 3 jets', 'Matched_objects_3J'),
##		('nMatched_objects_4J_hist', 'Matched Objects for 4 jets', 'Matched_objects_4J'),
##		('nMatched_objects_5PJ_hist', 'Matched Objects for 5+ jets', 'Matched_objects_5PJ')
##	]
##
##	for var, xaxis, names in nJets:
##		efficiencies = []
##		hist = asrootpy(myfile.Get(var)).Clone()
##		mean = hist.GetMean()
##		rms = hist.GetRMS()
##		total= hist.Integral()
##		plotter.set_histo_style(hist, color=defcol)
##		plotter.plot(hist, xtitle=xaxis, ytitle=defyax)
##		box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
##		plotter.save(names)
##		for i in range(0, hist.GetXaxis().GetNbins()):
##			efficiencies.append(format(hist.GetBinContent(i+1)/total, '.4f'))
##		print(xaxis, efficiencies)

#############################################################################################
if args.plot == "BTagging_Eff":
	xaxis = 'm_{t#bar t} [GeV]'
	yaxis = 'BTagging Efficiency'

	BTag = [
		('BTag_less_BHad_loose_pass', 'BTag_less_BHad_loose_fail', 'BTag_great_BHad_loose_pass', 'BTag_great_BHad_loose_fail', 'b_{h} Loose', 'BTag_BHad_loose_eff'),
		('BTag_less_BHad_medium_pass', 'BTag_less_BHad_medium_fail', 'BTag_great_BHad_medium_pass', 'BTag_great_BHad_medium_fail', 'b_{h} Medium', 'BTag_BHad_medium_eff'),
		('BTag_less_BHad_tight_pass', 'BTag_less_BHad_tight_fail', 'BTag_great_BHad_tight_pass', 'BTag_great_BHad_tight_fail', 'b_{h} Tight', 'BTag_BHad_tight_eff'),
		('BTag_less_BLep_loose_pass', 'BTag_less_BLep_loose_fail', 'BTag_great_BLep_loose_pass', 'BTag_great_BLep_loose_fail', 'b_{l} Loose', 'BTag_BLep_loose_eff'),
		('BTag_less_BLep_medium_pass', 'BTag_less_BLep_medium_fail', 'BTag_great_BLep_medium_pass', 'BTag_great_BLep_medium_fail', 'b_{l} Medium', 'BTag_BLep_medium_eff'),
		('BTag_less_BLep_tight_pass', 'BTag_less_BLep_tight_fail', 'BTag_great_BLep_tight_pass', 'BTag_great_BLep_tight_fail', 'b_{l} Tight', 'BTag_BLep_tight_eff')
	]

	labels = ['Merged', 'Unmerged']
#	colors = ['red', 'blue'] 
	for Merged_Pass, Merged_Fail, Unmerged_Pass, Unmerged_Fail, title, name in BTag:
		Merged_DR = []
		Unmerged_DR = []

		for folder, legend, colors in DRvals:
			to_draw = []
			
		#	Merged_Effs = []
			Merged_Passing = asrootpy(myfile.Get(folder+'/'+Merged_Pass)).Clone()
			Merged_Failing = asrootpy(myfile.Get(folder+'/'+Merged_Fail)).Clone()
			Merged_BTag_Eff = Merged_Passing/(Merged_Passing+Merged_Failing)
			plotter.set_histo_style(Merged_BTag_Eff, title=labels[0], color='red', markerstyle=20)
			plotter.plot(Merged_BTag_Eff)
			Merged_BTag_Eff.Draw("P")
		#	Merged_BTag_Eff.xaxis.set_title(xaxis)
		#	Merged_BTag_Eff.yaxis.set_title(yaxis)
			Merged_BTag_Eff.yaxis.range_user = 0, 1
		#	for i in range(0, Merged_BTag_Eff.GetXaxis().GetNbins()):
		#		Merged_Effs.append(format(Merged_BTag_Eff.GetBinContent(i+1), '.4f'))
		#	print(title+' Merged', Merged_Effs)
			to_draw.append(Merged_BTag_Eff)	

		#	Unmerged_Effs = []
			Unmerged_Passing = asrootpy(myfile.Get(folder+'/'+Unmerged_Pass)).Clone()
			Unmerged_Failing = asrootpy(myfile.Get(folder+'/'+Unmerged_Fail)).Clone()
			Unmerged_BTag_Eff = Unmerged_Passing/(Unmerged_Passing+Unmerged_Failing)
			plotter.set_histo_style(Unmerged_BTag_Eff, title=labels[1], color='blue', markerstyle=20)
			plotter.plot(Unmerged_BTag_Eff)
			Unmerged_BTag_Eff.Draw("P")
	#		Unmerged_BTag_Eff.xaxis.set_title(xaxis)
	#		Unmerged_BTag_Eff.yaxis.set_title(yaxis)
			Unmerged_BTag_Eff.yaxis.range_user = 0, 1
	#		for i in range(0, Unmerged_BTag_Eff.GetXaxis().GetNbins()):
	#			Unmerged_Effs.append(format(Unmerged_BTag_Eff.GetBinContent(i+1), '.4f'))
	#		print(title+' Unmerged', Unmerged_Effs)
			to_draw.append(Unmerged_BTag_Eff)	

			plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), y_range=(0,1), xtitle=xaxis, ytitle=yaxis)
			box = plotter.make_text_box(title+'\n%s' % legend, position='NW')
			plotter.save(name+'_'+folder)

			plotter.set_histo_style(Merged_BTag_Eff, title=legend, color=colors, markerstyle=20)
			plotter.plot(Merged_BTag_Eff)
			Merged_BTag_Eff.Draw("P")
			Merged_DR.append(Merged_BTag_Eff)
		
			plotter.set_histo_style(Unmerged_BTag_Eff, title=legend, color=colors, markerstyle=20)
			plotter.plot(Unmerged_BTag_Eff)
			Unmerged_BTag_Eff.Draw("P")
			Unmerged_DR.append(Unmerged_BTag_Eff)

		plotter.overlay(Merged_DR, legend_def=LegendDefinition(position='NE'), y_range=(0, 1), xtitle=xaxis, ytitle=yaxis)
		box = plotter.make_text_box('Merged '+title, position='NW')
		plotter.save('Merged_'+name)

		plotter.overlay(Unmerged_DR, legend_def=LegendDefinition(position='NE'), y_range=(0, 1), xtitle=xaxis, ytitle=yaxis)
		box = plotter.make_text_box('Unmerged '+title, position='NW')
		plotter.save('Unmerged_'+name)


#####################################################################################################################
#for var, xaxis, title in helframe: #variable name in histogram, xaxis title, saved title of png
#	j = 0
#	to_draw = []
#	for i in types:
#		hist = asrootpy(myfile.Get(var % i)).Clone()
#		plotter.set_histo_style(hist, title=labels[j], color=colors[j])
#		to_draw.append(hist)
#		j = j+1
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position=lpos), xtitle=xaxis, ytitle=yaxis)
#	plotter.save(title)

#	BTag_Unmerged = [
#		('BTag_gp4_BHad_loose_pass', 'BTag_gp4_BHad_loose_fail', 'm_{t#bar t} [GeV]', 'b_{h} Unmerged Loose WP Efficiency', 'BTag_Unmerged_BHad_loose_eff'),
#		('BTag_gp4_BHad_medium_pass', 'BTag_gp4_BHad_medium_fail', 'm_{t#bar t} [GeV]', 'b_{h} Unmerged Medium WP Efficiency', 'BTag_Unmerged_BHad_medium_eff'),
#		('BTag_gp4_BHad_tight_pass', 'BTag_gp4_BHad_tight_fail', 'm_{t#bar t} [GeV]', 'b_{h} Unmerged Tight WP Efficiency', 'BTag_Unmerged_BHad_tight_eff'),
#		('BTag_gp4_BLep_loose_pass', 'BTag_gp4_BLep_loose_fail', 'm_{t#bar t} [GeV]', 'b_{l} Unmerged Loose WP Efficiency', 'BTag_Unmerged_BLep_loose_eff'),
#		('BTag_gp4_BLep_medium_pass', 'BTag_gp4_BLep_medium_fail', 'm_{t#bar t} [GeV]', 'b_{l} Unmerged Medium WP Efficiency', 'BTag_Unmerged_BLep_medium_eff'),
#		('BTag_gp4_BLep_tight_pass', 'BTag_gp4_BLep_tight_fail', 'm_{t#bar t} [GeV]', 'b_{l} Unmerged Tight WP Efficiency', 'BTag_Unmerged_BLep_tight_eff')
#	]
#	for Pass, Fail, xaxis, yaxis, name in BTag_Unmerged:
#		efficiencies = []
#		
#		Passing = asrootpy(myfile.Get(Pass)).Clone()
#		Failing = asrootpy(myfile.Get(Fail)).Clone()
#		BTag_Eff = Passing/(Passing+Failing)
#		plotter.set_histo_style(BTag_Eff, color=defcol, markerstyle=20)
#		plotter.plot(BTag_Eff)
#		BTag_Eff.Draw("P")
#		BTag_Eff.xaxis.set_title(xaxis)
#		BTag_Eff.yaxis.set_title(yaxis)
#		BTag_Eff.yaxis.range_user = 0, 1
#		plotter.save(name)
#		
#		for i in range(0, BTag_Eff.GetXaxis().GetNbins()):
#			efficiencies.append(format(BTag_Eff.GetBinContent(i+1), '.4f'))
#		print(yaxis, efficiencies)

#############################################################################################
if args.plot == "Merged_Perms": #WJet merged with BJet
	
	Norm_Merged_Combined_Masses = [
		('Merged_BHadWJa_perm_and_WJb_mass', 'm_{merged+WJb} [GeV]', 'Merged b_{h} & WJa', '_BHadWJa_'),
		('Merged_BLepWJa_perm_and_WJb_mass', 'm_{merged+WJb} [GeV]', 'Merged b_{l} & WJa', '_BLepWJa_'),
		('Merged_BHadWJb_perm_and_WJa_mass', 'm_{merged+WJa} [GeV]', 'Merged b_{h} & WJb', '_BHadWJb_'),
		('Merged_BLepWJb_perm_and_WJa_mass', 'm_{merged+WJa} [GeV]', 'Merged b_{l} & WJb', '_BLepWJb_')
	]

	for var, xtitles, merged, name in Norm_Merged_Combined_Masses:
		to_draw = []
		for folder, legend, colors in DRvals:
			hist = asrootpy(myfile.Get(folder+'/'+var)).Clone()
			plotter.set_histo_style(hist, title=legend, color=colors)
			hist.xaxis.range_user = 0, 300
			to_draw.append(hist)

		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(0, 300), xtitle=xtitles, ytitle='Normalized', drawstyle='hist')
		box = plotter.make_text_box('Combined Masses\n%s' % merged, position='NW')
		plotter.save('Merged_Norm'+name+'Combined_Masses')


	Merged_B_Masses = [
		('Merged_BHadWJa_perm_mass', 'Merged_BHadWJb_perm_mass', 'm_{merged} [GeV]', 'b_{h}', 'Merged_BHad_Masses_'),
		('Merged_BLepWJa_perm_mass', 'Merged_BLepWJb_perm_mass', 'm_{merged} [GeV]', 'b_{l}', 'Merged_BLep_Masses_')
	]	

	for wja, wjb, xaxis, b_type, name in Merged_B_Masses:
		for folder, legend, colors in DRvals:
			to_draw = []
			WJa = asrootpy(myfile.Get(folder+'/'+wja)).Clone()
			plotter.set_histo_style(WJa, title=b_type+' & WJa' , color='red')
			WJa.xaxis.range_user = 0, 300
			to_draw.append(WJa)
	
			WJb = asrootpy(myfile.Get(folder+'/'+wjb)).Clone()
			plotter.set_histo_style(WJb, title=b_type+' & WJb' , color='blue')
			WJb.xaxis.range_user = 0, 300
			to_draw.append(WJb)
	
			plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(0, 300), xtitle=xaxis, ytitle=defyax, drawstyle='hist')
			box = plotter.make_text_box('Merged Jet Masses\n%s' % legend, position='NW')
			plotter.save(name+folder)


	Other_Masses = [
		('Merged_WJb_mass', 'Merged_WJa_mass', 'm_{other W} [GeV]', 'Merged_OtherW_Masses_')
	]

	for wjb, wja, xaxis, name in Other_Masses:
		for folder, legend, colors in DRvals:
			to_draw = []
			WJb = asrootpy(myfile.Get(folder+'/'+wjb)).Clone()
			plotter.set_histo_style(WJb, title='WJb', color='red')
			WJb.xaxis.range_user = 0, 150
			to_draw.append(WJb)
	
			WJa = asrootpy(myfile.Get(folder+'/'+wja)).Clone()
			plotter.set_histo_style(WJa, title='WJa', color='blue')
			WJa.xaxis.range_user = 0, 150
			to_draw.append(WJa)
	
			plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(0, 150), xtitle=xaxis, ytitle=defyax, drawstyle='hist')
			box = plotter.make_text_box('Other W Masses\n%s' % legend, position='NW')
			plotter.save(name+folder)

	Merged_B_DR = [
		('Merged_BHadWJa_perm_DRLepBLep', 'Merged_BHadWJb_perm_DRLepBLep', '#Delta R (l, b_{l})', 'b_{h}', 'Merged_BHadWJ_LepBLep_'),
		('Merged_BLepWJa_perm_DRLepBLep', 'Merged_BLepWJb_perm_DRLepBLep', '#Delta R (l, b_{l})', 'b_{l}', 'Merged_BLepWJ_LepBLep_')
	]

	for wja, wjb, xaxis, b_type, name in Merged_B_DR:
		for folder, legend, colors in DRvals:
			to_draw = []
			WJa = asrootpy(myfile.Get(folder+'/'+wja)).Clone()
			plotter.set_histo_style(WJa, title=b_type+' & WJa', color='red')
			to_draw.append(WJa)
			
			WJb = asrootpy(myfile.Get(folder+'/'+wjb)).Clone()
			plotter.set_histo_style(WJb, title=b_type+' & WJb', color='blue')
			to_draw.append(WJb)
			
			plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle=xaxis, ytitle=defyax, drawstyle='hist')
			box = plotter.make_text_box('Merged Jet\n%s' % legend, position='NW')
			plotter.save(name+folder)


	Unmerged_Jets = [
		('Unmerged_BHad_DRLepBLep', 'Unmerged_BLep_DRLepBLep', 'Unmerged_WJa_DRLepBLep', 'Unmerged_WJb_DRLepBLep', '#Delta R (l, b_{l})', '#Delta R (l, b_{l})', '#Delta R (l, b_{l})' , '#Delta R (l, b_{l})'),
		('Unmerged_BHad_mass', 'Unmerged_BLep_mass', 'Unmerged_WJa_mass', 'Unmerged_WJb_mass', 'm_{b_{h}} [GeV]', 'm_{b_{l}} [GeV]', 'm_{WJa} [GeV]', 'm_{WJb} [GeV]')
	]
	
	for bhad, blep, wja, wjb, bhad_xtitles, blep_xtitles, wja_xtitles, wjb_xtitles in Unmerged_Jets:
		BHad_to_draw = []
		BLep_to_draw = []
		WJa_to_draw = []
		WJb_to_draw = []

		for folder, legend, colors in DRvals:
			BHad = asrootpy(myfile.Get(folder+'/'+bhad)).Clone()
			plotter.set_histo_style(BHad, title=legend, color=colors)
			BHad_to_draw.append(BHad)
	
			BLep = asrootpy(myfile.Get(folder+'/'+blep)).Clone()
			plotter.set_histo_style(BLep, title=legend, color=colors)
			BLep_to_draw.append(BLep)
	
			WJa = asrootpy(myfile.Get(folder+'/'+wja)).Clone()
			plotter.set_histo_style(WJa, title=legend, color=colors)
			WJa_to_draw.append(WJa)
	
			WJb = asrootpy(myfile.Get(folder+'/'+wjb)).Clone()
			plotter.set_histo_style(WJb, title=legend, color=colors)
			WJb_to_draw.append(WJb)
	
		plotter.overlay(BHad_to_draw, legend_def=LegendDefinition(position='NE'), xtitle=bhad_xtitles, ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Unmerged b_{h}', position='NW')
		plotter.save(bhad)
		
		plotter.overlay(BLep_to_draw, legend_def=LegendDefinition(position='NE'), xtitle=blep_xtitles, ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Unmerged b_{l}', position='NW')
		plotter.save(blep)
		
		plotter.overlay(WJa_to_draw, legend_def=LegendDefinition(position='NE'), xtitle=wja_xtitles, ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Unmerged WJa', position='NW')
		plotter.save(wja)
		
		plotter.overlay(WJb_to_draw, legend_def=LegendDefinition(position='NE'), xtitle=wjb_xtitles, ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Unmerged WJb', position='NW')
		plotter.save(wjb)


	M_U_B_Mass = [#merged/unmerged B mass
		('Unmerged_BHad_mass', 'Merged_BHadWJa_perm_mass', 'Merged_BHadWJb_perm_mass', 'm_{jet} [GeV]', 'b_{h}', 'BHad_jet_masses_'),
		('Unmerged_BLep_mass', 'Merged_BLepWJa_perm_mass', 'Merged_BLepWJb_perm_mass', 'm_{jet} [GeV]', 'b_{l}', 'BLep_jet_masses_')
	]

	for unmerged_b, merged_bwja, merged_bwjb, xaxis, b_type, name in M_U_B_Mass:
		for folder, legend, colors in DRvals:
			to_draw = []
			Merged_BWJa = asrootpy(myfile.Get(folder+'/'+merged_bwja)).Clone()
			plotter.set_histo_style(Merged_BWJa, title='Merged '+b_type+' & WJa', color='blue')
			Merged_BWJa.xaxis.range_user = 0, 200
			to_draw.append(Merged_BWJa)

			Merged_BWJb = asrootpy(myfile.Get(folder+'/'+merged_bwjb)).Clone()
			plotter.set_histo_style(Merged_BWJb, title='Merged '+b_type+' & WJb', color='black')
			Merged_BWJb.xaxis.range_user = 0, 200
			to_draw.append(Merged_BWJb)

			Unmerged_B = asrootpy(myfile.Get(folder+'/'+unmerged_b)).Clone()
			plotter.set_histo_style(Unmerged_B, title='Unmerged '+b_type, color='red')
			Unmerged_B.xaxis.range_user = 0, 200
			to_draw.append(Unmerged_B)

			plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(0, 200),  xtitle=xaxis, drawstyle='hist')
			box = plotter.make_text_box('%s Mass\n%s' % (b_type, legend), position='NW')
			plotter.save(name+folder)


	Mass_Div_Pt = [# Merged/Unmerged jet mass/pt
		('Unmerged_BHad_massDivpt', 'Unmerged_BLep_massDivpt', 'Unmerged_WJa_massDivpt', 'Unmerged_WJb_massDivpt', 'Merged_BHadWJa_massDivpt', 'Merged_BHadWJb_massDivpt', 'Merged_BLepWJa_massDivpt', 'Merged_BLepWJb_massDivpt')
	]

	Unmerged_BHad_to_draw = []
	Unmerged_BLep_to_draw = []
	Unmerged_WJa_to_draw = []
	Unmerged_WJb_to_draw = []
	Merged_BHadWJa_to_draw = []
	Merged_BHadWJb_to_draw = []
	Merged_BLepWJa_to_draw = []
	Merged_BLepWJb_to_draw = []
	for U_bhad, U_blep, U_wja, U_wjb, M_bhadwja, M_bhadwjb, M_blepwja, M_blepwjb in Mass_Div_Pt:
		for folder, legend, colors in DRvals:
			U_BHad = asrootpy(myfile.Get(folder+'/'+U_bhad)).Clone()
			plotter.set_histo_style(U_BHad, title=legend, color=colors)
			Unmerged_BHad_to_draw.append(U_BHad)
	
			U_BLep = asrootpy(myfile.Get(folder+'/'+U_blep)).Clone()
			plotter.set_histo_style(U_BLep, title=legend, color=colors)
			Unmerged_BLep_to_draw.append(U_BLep)
	
			U_WJa = asrootpy(myfile.Get(folder+'/'+U_wja)).Clone()
			plotter.set_histo_style(U_WJa, title=legend, color=colors)
			Unmerged_WJa_to_draw.append(U_WJa)
	
			U_WJb = asrootpy(myfile.Get(folder+'/'+U_wjb)).Clone()
			plotter.set_histo_style(U_WJb, title=legend, color=colors)
			Unmerged_WJb_to_draw.append(U_WJb)
	
			M_BHadWJa = asrootpy(myfile.Get(folder+'/'+M_bhadwja)).Clone()
			plotter.set_histo_style(M_BHadWJa, title=legend, color=colors)
			Merged_BHadWJa_to_draw.append(M_BHadWJa)
	
			M_BHadWJb = asrootpy(myfile.Get(folder+'/'+M_bhadwjb)).Clone()
			plotter.set_histo_style(M_BHadWJb, title=legend, color=colors)
			Merged_BHadWJb_to_draw.append(M_BHadWJb)
	
			M_BLepWJa = asrootpy(myfile.Get(folder+'/'+M_blepwja)).Clone()
			plotter.set_histo_style(M_BLepWJa, title=legend, color=colors)
			Merged_BLepWJa_to_draw.append(M_BLepWJa)
	
			M_BLepWJb = asrootpy(myfile.Get(folder+'/'+M_blepwjb)).Clone()
			plotter.set_histo_style(M_BLepWJb, title=legend, color=colors)
			Merged_BLepWJb_to_draw.append(M_BLepWJb)
	
		plotter.overlay(Unmerged_BHad_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='m/p_{t}', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Unmerged b_{h}', position='NW')
		plotter.save('Unmerged_BHad_massDivpt')
	
		plotter.overlay(Unmerged_BLep_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='m/p_{t}', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Unmerged b_{l}', position='NW')
		plotter.save('Unmerged_BLep_massDivpt')
	
		plotter.overlay(Unmerged_WJa_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='m/p_{t}', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Unmerged WJa', position='NW')
		plotter.save('Unmerged_WJa_massDivpt')
	
		plotter.overlay(Unmerged_WJb_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='m/p_{t}', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Unmerged WJb', position='NW')
		plotter.save('Unmerged_WJb_massDivpt')
	
		plotter.overlay(Merged_BHadWJa_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='m/p_{t}', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Merged b_{h} & WJa', position='NW')
		plotter.save('Merged_BHadWJa_massDivpt')
	
		plotter.overlay(Merged_BHadWJb_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='m/p_{t}', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Merged b_{h} & WJb', position='NW')
		plotter.save('Merged_BHadWJb_massDivpt')
	
		plotter.overlay(Merged_BLepWJa_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='m/p_{t}', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Merged b_{l} & WJa', position='NW')
		plotter.save('Merged_BLepWJa_massDivpt')
	
		plotter.overlay(Merged_BLepWJb_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='m/p_{t}', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Merged b_{l} & WJb', position='NW')
		plotter.save('Merged_BLepWJb_massDivpt')



	mBWJ_vs_mB = [ #mass of b+wjet vs mass of b
		('Merged_BHadWJet_mass_vs_BHad_mass', 'Unmerged_BHadWJet_mass_vs_BHad_mass', 'Unmerged_BHadWJet_highest_mass_vs_BHad_mass', 'b_{h}'),
		('Merged_BLepWJet_mass_vs_BLep_mass', 'Unmerged_BLepWJet_mass_vs_BLep_mass', 'Unmerged_BLepWJet_highest_mass_vs_BLep_mass', 'b_{l}')
	]
	
	for merged, unmerged, unmerged_highest, b_type in mBWJ_vs_mB:
		for folder, legend, colors in DRvals:
			Merged = asrootpy(myfile.Get(folder+'/'+merged)).Clone()
			plotter.plot(Merged)
			Merged.Draw('colz')
			Merged.xaxis.set_title('m_{%s} [GeV]' % b_type)
			Merged.yaxis.set_title('m_{%s+jet} [GeV]' % b_type)
			box = plotter.make_text_box('Merged '+legend, position='SE')
			plotter.save(merged+'_'+folder)
			
			Unmerged = asrootpy(myfile.Get(folder+'/'+unmerged)).Clone()
			plotter.plot(Unmerged)
			Unmerged.Draw('colz')
			Unmerged.xaxis.set_title('m_{%s} [GeV]' % b_type)
			Unmerged.yaxis.set_title('m_{%s+jet} [GeV]' % b_type)
			box = plotter.make_text_box('Unmerged '+legend, position='SE')
			plotter.save(unmerged+'_'+folder)
			
			Unmerged_highest = asrootpy(myfile.Get(folder+'/'+unmerged_highest)).Clone()
			plotter.plot(Unmerged_highest)
			Unmerged_highest.Draw('colz')
			Unmerged_highest.xaxis.set_title('m_{%s} [GeV]' % b_type)
			Unmerged_highest.yaxis.set_title('m_{%s+jet} [GeV]' % b_type)
			box = plotter.make_text_box('Unmerged Highest Mass\nWJet %s' % legend, position='SE')
			plotter.save(unmerged_highest+'_'+folder)

	Merged_AllEvents_Ratio = [
		('All_BHad_events_vs_mttbar', 'Merged_BHadBLep_vs_mttbar', 'b_{h} & b_{l}', 'Merged_BHadBLep_Frac'),
		('All_BHad_events_vs_mttbar', 'Merged_BHadWJa_vs_mttbar', 'b_{h} & WJa', 'Merged_BHadWJa_Frac'),
		('All_BHad_events_vs_mttbar', 'Merged_BHadWJb_vs_mttbar', 'b_{h} & WJb', 'Merged_BHadWJb_Frac'),
		('All_BLep_events_vs_mttbar', 'Merged_BHadBLep_vs_mttbar', 'b_{l} & b_{h}', 'Merged_BLepBHad_Frac'),
		('All_BLep_events_vs_mttbar', 'Merged_BLepWJa_vs_mttbar', 'b_{l} & WJa', 'Merged_BLepWJa_Frac'),
		('All_BLep_events_vs_mttbar', 'Merged_BLepWJb_vs_mttbar', 'b_{l} & WJb', 'Merged_BLepWJb_Frac'),
		('All_WJa_events_vs_mttbar', 'Merged_BHadWJa_vs_mttbar', 'WJa & b_{h}', 'Merged_WJaBHad_Frac'),
		('All_WJa_events_vs_mttbar', 'Merged_BLepWJa_vs_mttbar', 'WJa & b_{l}', 'Merged_WJaBLep_Frac'),
		('All_WJa_events_vs_mttbar', 'Merged_WJaWJb_vs_mttbar', 'WJa & WJb', 'Merged_WJaWJb_Frac'),
		('All_WJb_events_vs_mttbar', 'Merged_WJaWJb_vs_mttbar', 'WJb & WJa', 'Merged_WJbWJa_Frac'),
		('All_WJb_events_vs_mttbar', 'Merged_BHadWJb_vs_mttbar', 'WJb & b_{h}', 'Merged_WJbBHad_Frac'),
		('All_WJb_events_vs_mttbar', 'Merged_BLepWJb_vs_mttbar', 'WJb & b_{l}', 'Merged_WJbBLep_Frac')
	]

	for all_events, merged, objects, name in Merged_AllEvents_Ratio:

		to_draw = []
		
		for folder, legend, colors in DRvals:
			All = asrootpy(myfile.Get(folder+'/'+all_events)).Clone()
			Merged = asrootpy(myfile.Get(folder+'/'+merged)).Clone()
			Merged_Frac = Merged/All
			plotter.set_histo_style(Merged_Frac, title=legend, color=colors, markerstyle=20)
			plotter.plot(Merged_Frac)
			Merged_Frac.Draw("P")
		#	Merged_Frac.xaxis.set_title(xaxis)
		#	Merged_Frac.yaxis.set_title(yaxis)
			Merged_Frac.yaxis.range_user = 0, 1
			to_draw.append(Merged_Frac)

		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'),y_range=(0,1), xtitle='m_{t#bart}', ytitle='Fraction of Merged Events/ Total Events')
		box = plotter.make_text_box(objects, position='NW')
		plotter.save(name)


	Invar_mass_BWJet_Kin = [
		('Merged_BHadWJet_l100_pt', 'Merged_BHadWJet_l100_eta', 'Merged_BLepWJet_l100_pt', 'Merged_BLepWJet_l100_eta')
	]

	for bhad_pt, bhad_eta, blep_pt, blep_eta in Invar_mass_BWJet_Kin:

		BHad_Pt_to_draw = []
		BHad_Eta_to_draw = []
		BLep_Pt_to_draw = []
		BLep_Eta_to_draw = []

		for folder, legend, colors in DRvals:
			BHad_Pt = asrootpy(myfile.Get(folder+'/'+bhad_pt)).Clone()
			plotter.set_histo_style(BHad_Pt, title=legend, color=colors)
			BHad_Pt.xaxis.range_user = 0, 600
			BHad_Pt_to_draw.append(BHad_Pt)

			BHad_Eta = asrootpy(myfile.Get(folder+'/'+bhad_eta)).Clone()
			plotter.set_histo_style(BHad_Eta, title=legend, color=colors)
			BHad_Eta.xaxis.range_user = -2.4, 2.4
			BHad_Eta_to_draw.append(BHad_Eta)

			BLep_Pt = asrootpy(myfile.Get(folder+'/'+blep_pt)).Clone()
			plotter.set_histo_style(BLep_Pt, title=legend, color=colors)
			BLep_Pt.xaxis.range_user = 0, 600
			BLep_Pt_to_draw.append(BLep_Pt)

			BLep_Eta = asrootpy(myfile.Get(folder+'/'+blep_eta)).Clone()
			plotter.set_histo_style(BLep_Eta, title=legend, color=colors)
			BLep_Eta.xaxis.range_user = -2.4, 2.4
			BLep_Eta_to_draw.append(BLep_Eta)

		plotter.overlay(BHad_Pt_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='Merged b_{h} p_{T}', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('m_{merged+WJet} < 100 GeV', position='NW')
		plotter.save('Merged_BHadWJet_l100_pt')

		plotter.overlay(BHad_Eta_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='Merged b_{h} #eta', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('m_{merged+WJet} < 100 GeV', position='NW')
		plotter.save('Merged_BHadWJet_l100_eta')

		plotter.overlay(BLep_Pt_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='Merged b_{l} p_{T}', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('m_{merged+WJet} < 100 GeV', position='NW')
		plotter.save('Merged_BLepWJet_l100_pt')

		plotter.overlay(BLep_Eta_to_draw, legend_def=LegendDefinition(position='NE'), xtitle='Merged b_{l} #eta', ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('m_{merged+WJet} < 100 GeV', position='NW')
		plotter.save('Merged_BLepWJet_l100_eta')


#############################################################################################
if args.plot == "Delta_Plots":
	# Delta (Reco - Gen) Plots
	Pt = [
	      ('delta_pt_lep_hist', '#mu #Delta p_{t}', '#Delta p_{t} #mu', 'Pt_Mu'),
	      ('delta_pt_BHad_hist', 'b_{had} #Delta p_{t}', '#Delta p_{t} b_{had}', 'Pt_BHad'),
	      ('delta_pt_BLep_hist', 'b_{lep} #Delta p_{t}', '#Delta p_{t} b_{lep}', 'Pt_BLep'),
	      ('delta_pt_WJa_hist', 'WJa #Delta p_{t}', '#Delta p_{t} WJa', 'Pt_WJa'),
	      ('delta_pt_WJb_hist', 'WJb #Delta p_{t}', '#Delta p_{t} WJb', 'Pt_WJb')
	]
	
	for var, xaxis, titles, names in Pt:
	      hist = asrootpy(myfile.Get(var)).Clone()
	#	hist.xaxis.range_user = -250,250
	 #     height = hist.GetMaximum()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      hist.Fit("gaus")
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, title=titles, xtitle=xaxis, ytitle=defyax)
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	#	box.Draw()
	#	plotter.keep.append(box)
	      plotter.save(names)
	
	Phi = [
	      ('delta_phi_lep_hist', '#mu #Delta #phi', '#Delta#phi #mu', 'Phi_Mu'),
	      ('delta_phi_BHad_hist', 'b_{had} #Delta #phi', '#Delta#phi b_{had}', 'Phi_BHad'),
	      ('delta_phi_BLep_hist', 'b_{lep} #Delta #phi', '#Delta#phi b_{lep}', 'Phi_BLep'),
	      ('delta_phi_WJa_hist', 'WJa #Delta #phi', '#Delta#phi WJa', 'Phi_WJa'),
	      ('delta_phi_WJb_hist', 'WJb #Delta #phi', '#Delta#phi WJb', 'Phi_WJb')
	]
	
	for var, xaxis, titles, names in Phi:
	      hist = asrootpy(myfile.Get(var)).Clone()
	#      height = hist.GetMaximum()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	#	hist.xaxis.range_user = -3.2, 3.2
	      
	      hist.Fit("gaus")
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, title=titles, xtitle=xaxis, ytitle=defyax)
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	#	box.Draw()
	#	plotter.keep.append(box)
	      plotter.save(names)
	
	Rapidity = [
	      ('delta_y_lep_hist', '#mu #Delta y', '#Delta y #mu', 'Rapidity_Mu'),
	      ('delta_y_BHad_hist', 'b_{had} #Delta y', '#Delta y b_{had}', 'Rapidity_BHad'),
	      ('delta_y_BLep_hist', 'b_{lep} #Delta y', '#Delta y b_{lep}', 'Rapidity_BLep'),
	      ('delta_y_WJa_hist', 'WJa #Delta y', '#Delta y WJa', 'Rapidity_WJa'),
	      ('delta_y_WJb_hist', 'WJb #Delta y', '#Delta y WJb', 'Rapidity_WJb')
	]
	
	for var, xaxis, titles, names in Rapidity:
	      hist = asrootpy(myfile.Get(var)).Clone()
	#      height = hist.GetMaximum()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	#	hist.xaxis.range_user = -2, 2
	
	      hist.Fit("gaus")
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, title=titles, xtitle=xaxis, ytitle=defyax)
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	#	box.Draw()
	#	plotter.keep.append(box)
	      plotter.save(names)

#################################################################################################	
if args.plot == "Muon":
	## Muon Plots
	Lep_ISO = [
	#	('ElISO_hist', 'e ISO', 'Electron ISO', 'ISO_El'),
		('MuISO_hist', '#mu ISO', 'Muon ISO', 'ISO_Mu')
	]
	
	for var, xaxis, titles, names in Lep_ISO:
		hist = asrootpy(myfile.Get(var)).Clone()
	#	RMS = hist.GetRMS()
		plotter.set_histo_style(hist, color=defcol)
	#	print("RMS = %f" % RMS)
	#	plotter.make_text_box("RMS = %f" % RMS, position='NE')
		plotter.plot(hist, title=titles, xtitle=xaxis, ytitle=defyax)
		plotter.save(names)
	
	##Lep_ISO_Zoom = [
	###	('ElISO_hist', 'e ISO', 'Electron ISO', 'ISO_El_Zoom'),
	##	('MuISO_hist', '#mu ISO', 'Muon ISO', 'ISO_Mu_Zoom')
	##]
	##
	##for var, xaxis, titles, names in Lep_ISO_Zoom:
	##	hist = asrootpy(myfile.Get(var)).Clone()
	##	hist.xaxis.range_user = 0, 4
	###	RMS = hist.GetRMS()
	##	plotter.set_histo_style(hist, color=defcol)
	###	print("RMS = %f" % RMS)
	###	plotter.make_text_box("RMS = %f" % RMS, position='NE')
	##	plotter.plot(hist, title=titles, xtitle=xaxis, ytitle=defyax)
	##	plotter.save(names)
	#
	MassTTbar = [
		('Mu_Mttbar_hist', 'm_{t#bar t} [GeV]', 'Mass tt_bar for #mu Events', 'Muon_System_Mass')
	]
	
	for var, xaxis, titles, names in MassTTbar:
		hist = asrootpy(myfile.Get(var)).Clone()
	#	hist.xaxis.range_user = 150, 700
		mean = hist.GetMean()
		rms = hist.GetRMS()
		plotter.set_histo_style(hist, color=defcol)
	#	print("RMS = %f" % RMS)
		plotter.plot(hist, title=titles, xtitle=xaxis, ytitle=defyax)
		box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
		plotter.save(names)
	
	nJets = [
		('nJets_hist', 'nJets', 'nJets_hist')
	]

	for var, xaxis, names in nJets:
		hist = asrootpy(myfile.Get(var)).Clone()
		plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=defyax)
		plotter.save(names)

	Hists_2D = [
		('LepBLep_DR_vs_mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (l,b_{l})','#Delta R (l,b_{l}) vs m_{t#bar t}', 'LepBLep_DR_vs_mtt'),
		('MuISO_vs_mtt_hist', 'm_{t#bar t} [GeV]', '#mu ISO', '#mu ISO vs m_{t#bar t}', 'MuISO_vs_mtt'),
		('MuDR_vs_MuRelDPT_hist', '#mu #Delta R (gen,reco)', '#mu Rel #Delta p_{t}', '#mu #Delta R vs #mu Rel #Delta p_{t}', 'MuDR_vs_MuRelDPT')
	]
	
	for var, xaxis,yaxis, titles, names in Hists_2D:
		hist = asrootpy(myfile.Get(var)).Clone()
	#	plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, title=titles, xtitle=xaxis, ytitle=yaxis)
		hist.Draw('colz')
		hist.xaxis.set_title(xaxis)
		hist.yaxis.set_title(yaxis)
		if var=='MuDR_vs_MuRelDPT_hist':
			hist.yaxis.range_user = 0, 10
		plotter.save(names)
	
	Muon_Efficiencies = [
		('Unmatched_Muons_hist', 'Matched_Muons_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', 'Muon_Matching_Efficiencies'),
		('Non_ISO_Mu_hist', 'ISO_Mu_hist', 'm_{t#bar t} [GeV]', 'Matched #mu ISO Efficiency', 'Matched_Muon_ISO_Efficiencies')
	]
	
	for unmatched, matched, xaxis, yaxis, name in Muon_Efficiencies:
		efficiencies = []
		
		Unmatched = asrootpy(myfile.Get(unmatched)).Clone()
		Matched = asrootpy(myfile.Get(matched)).Clone()
		Match_Eff = Matched/(Matched+Unmatched)
		plotter.set_histo_style(Match_Eff, color=defcol, markerstyle=20)
		plotter.plot(Match_Eff)
		Match_Eff.Draw("P")
		Match_Eff.xaxis.set_title(xaxis)
		Match_Eff.yaxis.set_title(yaxis)
		Match_Eff.yaxis.range_user = 0, 1
		plotter.save(name)
		
		for i in range(0, Match_Eff.GetXaxis().GetNbins()):
			efficiencies.append(format(Match_Eff.GetBinContent(i+1), '.4f' ))
		#	mass_range.append(Mu_Eff.GetXaxis().GetBinLowEdge(i+1))
		print(yaxis, efficiencies)

###############################################################################################################
if args.plot == "Jets":
	## Jet Plots	
	Jet_Efficiencies = [
		('Unmatched_WJa_hist', 'Matched_WJa_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', 'WJa_Matching_Efficiencies'),
		('Unmatched_WJb_hist', 'Matched_WJb_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', 'WJb_Matching_Efficiencies'),
		('Unmatched_BHad_hist', 'Matched_BHad_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', 'BHad_Matching_Efficiencies'),
		('Unmatched_BLep_hist', 'Matched_BLep_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', 'BLep_Matching_Efficiencies'),
		('Unmatched_2Wjets_hist', 'Matched_2Wjets_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', '2Wjets_Matching_Efficiencies')
	]
	
	for unmatched, matched, xaxis, yaxis, name in Jet_Efficiencies:
		efficiencies = []
		
		Unmatched = asrootpy(myfile.Get(unmatched)).Clone()
		Matched = asrootpy(myfile.Get(matched)).Clone()
		Match_Eff = Matched/(Matched+Unmatched)
		plotter.set_histo_style(Match_Eff, color=defcol, markerstyle=20)
		plotter.plot(Match_Eff)
		Match_Eff.Draw("P")
		Match_Eff.xaxis.set_title(xaxis)
		Match_Eff.yaxis.set_title(yaxis)
		Match_Eff.yaxis.range_user = 0, 1
		plotter.save(name)
		
		for i in range(0, Match_Eff.GetXaxis().GetNbins()):
			efficiencies.append(format(Match_Eff.GetBinContent(i+1), '.4f' ))
		#	mass_range.append(Mu_Eff.GetXaxis().GetBinLowEdge(i+1))
		print(yaxis, efficiencies)
	
	
	Jet_Zoom_Efficiencies = [
		('Unmatched_BHad_zoom_hist', 'Matched_BHad_zoom_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', 'BHad_Matching_Zoom_Efficiencies'),
		('Unmatched_BLep_zoom_hist', 'Matched_BLep_zoom_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', 'BLep_Matching_Zoom_Efficiencies'),
		('Unmatched_WJa_zoom_hist', 'Matched_WJa_zoom_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', 'WJa_Matching_Zoom_Efficiencies'),
		('Unmatched_WJb_zoom_hist', 'Matched_WJb_zoom_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', 'WJb_Matching_Zoom_Efficiencies'),
		('Unmatched_2Wjets_zoom_hist', 'Matched_2Wjets_zoom_hist', 'm_{t#bar t} [GeV]', 'Matching Efficiency', '2Wjets_Matching_Zoom_Efficiencies')
	]
	
	for unmatched, matched, xaxis, yaxis, name in Jet_Zoom_Efficiencies:
		efficiencies = []
		
		Unmatched = asrootpy(myfile.Get(unmatched)).Clone()
		Matched = asrootpy(myfile.Get(matched)).Clone()
		Match_Eff = Matched/(Matched+Unmatched)
		plotter.set_histo_style(Match_Eff, color=defcol, markerstyle=20)
		plotter.plot(Match_Eff)
		Match_Eff.Draw("P")
		Match_Eff.xaxis.set_title(xaxis)
		Match_Eff.yaxis.set_title(yaxis)
		Match_Eff.yaxis.range_user = 0, 1
		plotter.save(name)
		
		for i in range(0, Match_Eff.GetXaxis().GetNbins()):
			efficiencies.append(format(Match_Eff.GetBinContent(i+1), '.4f' ))
		#	mass_range.append(Mu_Eff.GetXaxis().GetBinLowEdge(i+1))
		print(yaxis, efficiencies)
	
	Unmatched_Jet_Vars = [
		('Unmatched_BLep_pt_hist', 'b_{lep} p_{t}', 'Unmatched b_{lep}', 'Unmatched_BLep_Pt'),
		('Unmatched_BHad_pt_hist', 'b_{had} p_{t}', 'Unmatched b_{had}', 'Unmatched_BHad_Pt'),
		('Unmatched_WJa_pt_hist', 'WJa p_{t}', 'Unmatched WJa', 'Unmatched_WJa_Pt'),
		('Unmatched_WJb_pt_hist', 'WJb p_{t}', 'Unmatched WJb', 'Unmatched_WJb_Pt'),
		('Unmatched_BLep_y_hist', 'b_{lep} y', 'Unmatched b_{lep}', 'Unmatched_BLep_Rapidity'),
		('Unmatched_BHad_y_hist', 'b_{had} y', 'Unmatched b_{had}', 'Unmatched_BHad_Rapidity'),
		('Unmatched_WJa_y_hist', 'WJa y', 'Unmatched WJa', 'Unmatched_WJa_Rapidity'),
		('Unmatched_WJb_y_hist', 'WJb y', 'Unmatched WJb', 'Unmatched_WJb_Rapidity')
	] 
	
	for var, xaxis, titles, names in Unmatched_Jet_Vars:
		hist = asrootpy(myfile.Get(var)).Clone()
		mean = hist.GetMean()
		rms = hist.GetRMS()
		if (var=='Unmatched_BLep_pt_hist') or (var =='Unmatched_BHad_pt_hist') or (var=='Unmatched_WJa_pt_hist') or (var=='Unmatched_WJb_pt_hist'):
			hist.xaxis.range_user = 0, 150
		plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, title=titles, xtitle=xaxis, ytitle=defyax)
		box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
		plotter.save(names)

	nJets = [
		('nJets_hist', 'nJets', 'nJets')
	]

	for var, xaxis, names in nJets:
		hist = asrootpy(myfile.Get(var)).Clone()
		plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=defyax)
		plotter.save(names)

	MassTTbar = [
		('Jets_Mttbar_hist', 'm_{t#bar t} [GeV]', 'Mass tt_bar for Jets Events', 'Jets_System_Mass')
	]
	
	for var, xaxis, titles, names in MassTTbar:
		hist = asrootpy(myfile.Get(var)).Clone()
	#	hist.xaxis.range_user = 150, 700
		mean = hist.GetMean()
		rms = hist.GetRMS()
		plotter.set_histo_style(hist, color=defcol)
	#	print("RMS = %f" % RMS)
		plotter.plot(hist, title=titles, xtitle=xaxis, ytitle=defyax)
		box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
		plotter.save(names)


##Top CM frame
#helframe = [
#	('helframe_costheta_%s', 'cos #theta*', 'helframe_costheta'),
#	('helframe_cosdelta_lep_%s', 'cos #delta(l,q)', 'helframe_cosdelta_lep'),
#	('helframe_prodcosth_lep_%s', 'cos #theta*_{l} cos #theta*_{q}', 'helframe_prodcosth_lep'),
#]
#
#for var, xaxis, title in helframe: #variable name in histogram, xaxis title, saved title of png
#	j = 0
#	to_draw = []
#	for i in types:
#		hist = asrootpy(myfile.Get(var % i)).Clone()
#		plotter.set_histo_style(hist, title=labels[j], color=colors[j])
#		to_draw.append(hist)
#		j = j+1
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position=lpos), xtitle=xaxis, ytitle=yaxis)
#	plotter.save(title)
#####
#
##W CM frame
#Wframe = [
#	('Wframe_costheta_%s', 'cos #theta*', 'Wframe_costheta'),
#	('Wframe_cosdelta_lep_%s', 'cos #delta(l,q)', 'Wframe_cosdelta_lep'),
#	('Wframe_prodcosth_lep_%s', 'cos #theta*_{l} cos #theta*{q}', 'Wframe_prodcosth_lep')
#]
#
#for var, xaxis, title in Wframe: #variable name in histogram, xaxis title, saved title of png
#	j = 0
#	to_draw = []
#	for i in types:
#		hist = asrootpy(myfile.Get(var % i)).Clone()
#		plotter.set_histo_style(hist, title=labels[j], color=colors[j])
#		to_draw.append(hist)
#		j = j+1
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position=lpos), xtitle=xaxis, ytitle=yaxis)
#	plotter.save(title)
#####
#
##Lab CM frame
#labframe = [
#	('labframe_cosdeltaphi_lep_%s', 'cos #Delta#varphi(l, q)', 'labframe_cosdeltaphi_lep'),
#	('labframe_deltaphi_lep_%s', '#Delta#varphi(l, q)', 'labframe_deltaphi_lep'),
#]
#
#for var, xaxis, title in labframe: #variable name in histogram, xaxis title, saved title of png
#	j = 0
#	to_draw = []
#	for i in types:
#		hist = asrootpy(myfile.Get(var % i)).Clone()
#		plotter.set_histo_style(hist, title=labels[j], color=colors[j])
#		to_draw.append(hist)
#		j = j+1
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position=lpos), xtitle=xaxis, ytitle=yaxis)
#	plotter.save(title)
#####
#
##Top CM delta plots
#helframe_deltvars = [
#	('delta_helframe_costheta_hist', 'HF cos #theta*', 'blue', 100, 'HF_costheta'),
#        ('delta_helframe_cosdelta_hist', 'HF cos #delta(l,q)', 'violet', 100, 'HF_cosdelta'),
#        ('delta_helframe_prodcosth_hist', 'HF cos #theta*_{l} cos #theta*_{q}', 'black', 100, 'HF_prodcosth'),
#        ('delta_labframe_cosdeltaphi_hist', 'LF cos #Delta#varphi(l, q)', 'cyan', 100, 'LF_cosdeltaphi'),
#        ('delta_labframe_deltaphi_hist', 'LF #Delta#varphi(l, q)', 'red', 100, 'LF_deltaphi'),
#	('frac_delta_E_hist', 'Rel. #Delta E', '#00960a', 100, 'Rel_delta_E')
#]
#
#headers = ["\tVariable", "Total Number of Events", "Number of Events > 0", "Efficiency"]
#print "\t".join(headers)
#to_draw = []
#posev = [] # number of events > 0 in hists
#events = [] # number of events in hists
#efficiencies = [] #positive/total events
#for var, labels, colors, bins, tablename in helframe_deltvars: #variable name in histogram, legend label, variable color, number of bins in hist
#	hist = asrootpy(myfile.Get(var)).Clone()
#	plotter.set_histo_style(hist, title=labels, color=colors)
#	to_draw.append(hist)
#	total= hist.Integral()
#	events.append(total)
#	positive = hist.Integral((bins/2+1), bins)
#	posev.append(positive)
#	efficiency = positive/total
#	efficiencies.append(efficiency)
#	print "\t%s\t\t%i\t\t\t%i\t\t%f" % (tablename, total, positive, efficiency)
#
##print "Number of Events > 0: ", posev
##print "Number of Events: ", events
##print "Efficiencies: ", efficiencies
#plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle='Helframe #Delta', ytitle=yaxis)
#plotter.save('Top_RF_Delta')
#
####
#
##Top CM fractional delta plots
#helframe_frac_deltvars = [
#        ('frac_delta_helframe_costheta_hist', 'HF cos #theta*', 'blue'),
#        ('frac_delta_labframe_cosdeltaphi_hist', 'LF cos #Delta#varphi(l, q)', 'cyan'),
#        ('frac_delta_labframe_deltaphi_hist', 'LF #Delta#varphi(l, q)', 'red'),
#        ('frac_delta_helframe_cosdelta_hist', 'HF cos #delta(l,q)', 'violet'),
#        ('frac_delta_helframe_prodcosth_hist', 'HF cos #theta*_{l} cos #theta*_{q}', 'black'),
#]
#
#to_draw = []
#
#for var, labels, colors in helframe_frac_deltvars: #variable name in histogram, legend label, variable color
#        hist = asrootpy(myfile.Get(var)).Clone()
#        plotter.set_histo_style(hist, title=labels, color=colors)
#        to_draw.append(hist)
#
#plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle='Helframe Rel. #Delta', ytitle=yaxis)
#plotter.save('Top_RF_Rel_Delta')
#####
#
##W CM delta plots 
#Wframe_deltvars = [
#        ('delta_Wframe_costheta_hist', 'WF cos #theta*', 'blue', 100, 'WF_costheta'),
#        ('delta_Wframe_cosdelta_hist', 'WF cos #delta(l,q)', 'violet', 100, 'WF_cosdelta'),
#        ('delta_Wframe_prodcosth_hist', 'WF cos #theta*_{l} cos #theta*_{q}', 'black', 100, 'WF_prodcosth'),
#        ('delta_labframe_cosdeltaphi_hist', 'LF cos #Delta#varphi(l, q)', 'cyan', 100, 'LF_cosdeltaphi'),
#        ('delta_labframe_deltaphi_hist', 'LF #Delta#varphi(l, q)', 'red', 100, 'LF_deltaphi'),
#        ('frac_delta_E_hist', 'Rel. #Delta E', '#00960a', 100, 'Rel_delta_E')
#]
#
#headers = ["\tVariable", "Total Number of Events", "Number of Events > 0", "Efficiency"]
#print "\t".join(headers)
#to_draw = []
#posev = [] # number of events > 0 in hists
#events = [] # number of events in hists
#efficiencies = [] #positive/total events
#for var, labels, colors, bins, tablename in Wframe_deltvars: #variable name in histogram, legend label, variable color, number of bins in hist
#	hist = asrootpy(myfile.Get(var)).Clone()
#	plotter.set_histo_style(hist, title=labels, color=colors)
#	to_draw.append(hist)
#	total= hist.Integral()
#	events.append(total)
#	positive = hist.Integral((bins/2+1), bins)
#	posev.append(positive)
#	efficiency = positive/total
#	efficiencies.append(efficiency)
#	print "\t%s\t\t%i\t\t\t%i\t\t%f" % (tablename, total, positive, efficiency)
#
#plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle='Wframe #Delta', ytitle=yaxis)
#plotter.save('W_RF_Delta')
#####
#
##W CM fractional delta plots
#Wframe_frac_deltvars = [
#        ('frac_delta_Wframe_costheta_hist', 'WF cos #theta*', 'blue'),
#        ('frac_delta_labframe_cosdeltaphi_hist', 'LF cos #Delta#varphi(l, q)', 'cyan'),
#        ('frac_delta_labframe_deltaphi_hist', 'LF #Delta#varphi(l, q)', 'red'),
#        ('frac_delta_Wframe_cosdelta_hist', 'WF cos #delta(l,q)', 'violet'),
#        ('frac_delta_Wframe_prodcosth_hist', 'WF cos #theta*_{l} cos #theta*_{q}', 'black'),
#]
#
#to_draw = []
#
#for var, labels, colors in Wframe_frac_deltvars: #variable name in histogram, legend label, variable color
#        hist = asrootpy(myfile.Get(var)).Clone()
#        plotter.set_histo_style(hist, title=labels, color=colors)
#        to_draw.append(hist)
#
#plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle='Wframe Rel. #Delta', ytitle=yaxis)
#plotter.save('W_RF_Rel_Delta')
#####
#
##### best permutation variables
#best_perm_vars = ['TopHadronicframe_%s', 'TopLeptonicframe_%s', 'delta_Wframe_jcosth_%s']
#best_perm = [
#	('RIGHT', 'semilep_visible_right', 'violet', 100),
#	('RIGHT_THAD', 'semilep_right_thad', 'cyan', 100),
#	('RIGHT_TLEP', 'semilep_right_tlep', 'red', 100),
#	('WRONG', 'semilep_wrong', '#00960a', 100),
##	('OTHER', 'other_tt_decay', 'blue', 100)
#]
#best_perm_names = ['THfr', 'TLfr', 'Wfr_jcosth']
##THfr_vars = [
##	('THfr_RIGHT', 'semilep_visible_right', 'violet'),
##	('THfr_RIGHT_THAD','semilep_right_thad', 'cyan'),
##	('THfr_RIGHT_TLEP', 'semilep_right_tlep', 'red'),
##	('THfr_WRONG', 'semilep_wrong', '#00960a'),
##	('THfr_OTHER', 'other_tt_decay', 'blue')
##]
#
#j = 0
#for i in best_perm_vars:
##	print i
##	headers = ["\tVariable", "Total Number of Events", "Number of Events > 0", "Efficiency"]
##	print "\t".join(headers)
##	posev = [] # number of events > 0 in hists
##	events = [] # number of events in hists
##	efficiencies = [] #positive/total events
#	to_draw = []
#	for var, labels, colors, bins in best_perm:
#		hist = asrootpy(normttJets.Get(i % var)).Clone()
#		plotter.set_histo_style(hist, title=labels, color=colors)
#	       	to_draw.append(hist)
##		print var
##		total= hist.Integral()
##		events.append(total)
##		positive = hist.Integral((bins/2+1), bins)
##		posev.append(positive)
##		if total !=0:
##			efficiency = positive/total
##		else:
##			efficiency = 0
##		efficiencies.append(efficiency)
##		print "\t%s\t\t%i\t\t\t%i\t\t%f" % (labels, total, positive, efficiency)
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), ytitle=yaxis)
#	plotter.save(best_perm_names[j])
#	j = j+1
#
##Other variables
#lep_vars = [
#	('lep_%s', ['charge', 'pt'], ['lepton charge', 'lepton pt'])
#]
#for var, types, xaxis in lep_vars: #variable name in histogram, type of lepton variable, xaxis title
#	j = 0
#	for i in types:
#		hist = asrootpy(myfile.Get(var % i)).Clone()
#		plotter.set_histo_style(hist, color=defcol)
#		plotter.plot(hist, xtitle=xaxis[j], ytitle=yaxis)
#		plotter.save(var % i)
#		j = j+1
#
#lep_nu_vars = [ #other angular variables
#	('helframe_costheta_%s', ['lep', 'nu'], 'cos #theta*'),
#	('Wframe_costheta_%s', ['lep', 'nu'], 'cos #theta*')
#]
#
#for var, types, xaxis in lep_nu_vars: #variable name, type of ang. variable, xaxis title
#	for i in types:
#		hist = asrootpy(myfile.Get(var % i)).Clone()
#		plotter.set_histo_style(hist, color=defcol)
#                plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
#                plotter.save(var % i)
#
##variables that have no pairs
#SA_vars = [
#	('delta_E_hist', '#Delta E', 'Delta_E'),
#	('frac_delta_E_hist', 'Frac. #Delta E', 'Frac_Delta_E'),
#	('nJets', 'nJets', 'nJets'),
#	('cut_flow', 'Cut Flow', 'Cut_Flow')
#]
#
#for var, xaxis, name in SA_vars: #variable name, xaxis title, png saved name
#	hist = asrootpy(myfile.Get(var)).Clone()
#	plotter.set_histo_style(hist, color=defcol)
#	plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
#	plotter.save(name)


print('~/nobackup/CMSSW_8_0_7_patch2/src/Analyses/URTTbar/htt_scripts/plots/jet_effs/%s/%s/%s/%s/' % (jobid, args.analysis, args.plot, args.sample))
