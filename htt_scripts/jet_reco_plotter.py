'''
Jet_Reco Analyzer Plotter macro
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
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Create plots using files from jet_reco.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full_Analysis).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, or Combined).')
parser.add_argument('plot', help='Choose type of plots to generate (Gen_Plots, Gen_Ratios, Matched_Objects, Delta_Plots, BTagging_Eff, Merged_Perms).')
args = parser.parse_args()

if args.analysis == "Test":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../%s.jet_reco.test.root' % args.sample, 'read')
		normfile = views.NormalizeView(root_open('../%s.jet_reco.test.root' % args.sample, 'read'))#normalized file
		plotter = BasePlotter(
			'plots/jet_reco/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Sample_ttJetsM0 = root_open('../ttJetsM0.jet_reco.test.root', 'read')
	#	Sample_ttJetsM700 = root_open('../ttJetsM700.jet_reco.test.root', 'read')
	#	Sample_ttJetsM1000 = root_open('../ttJetsM1000.jet_reco.test.root', 'read')
		plotter = BasePlotter(
			'plots/jet_reco/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)

elif args.analysis == "Full_Analysis":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../results/%s/jet_reco/%s.root' % (jobid, args.sample), 'read')
		normfile = views.NormalizeView(root_open('../results/%s/jet_reco/%s.root' % (jobid, args.sample), 'read'))
		plotter = BasePlotter(
			'plots/jet_reco/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Full_ttJetsM0 = root_open('../results/%s/jet_reco/ttJets.root' % jobid, 'read')
	#	Full_ttJetsM700 = root_open('../results/%s/jet_reco/ttJetsM700.root' % jobid, 'read')
	#	Full_ttJetsM1000 = root_open('../results/%s/jet_reco/ttJetsM1000.root' % jobid, 'read')
		plotter = BasePlotter(
			'plots/jet_reco/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)


##### Comparing right/wrong combinations from TTBarSolver settings ####
##open files
#NS_MD_file = root_open('../ttJets.jet_reco_NS_MassDisc.test.root', 'read') #only NS and MassDis right used in solver
#NS_MD_ratio_file = root_open('../ttJets.jet_reco_NS_MassDisc_Ratio.test.root', 'read') #NS and MassDis right/wrong used in solver
#NS_MD_Ang_file = root_open('../ttJets.jet_reco_NS_MassDisc_Ang.test.root', 'read') #NS, MassDisc, Ang vars right used in solver
#NS_MD_Ang_ratio_file = root_open('../ttJets.jet_reco_NS_MassDisc_Ang_Ratio.test.root', 'read') #NS, MassDisc, Ang vars right/wrong used in solver
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

#DRvals = {
#		('DRP4', '#Delta R < 0.4', 'red'),
#		('DRP5', '#Delta R < 0.5', 'blue'),
#		('DRP6', '#Delta R < 0.6', 'green'),
#		('DRP8', '#Delta R < 0.8', 'black')
#	}

DRvals = {
		('DRP4', '#Delta R < 0.4', 'red'),
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
#	Gen_DR = [
#		('DR_LepBHad', '#Delta R (l, b_{h})'),
#		('DR_LepBLep', '#Delta R (l, b_{l})'),
#		('DR_LepWJa', '#Delta R (l, WJa)'),
#		('DR_LepWJb', '#Delta R (l, WJb)'),
#		('DR_BHadBLep', '#Delta R (b_{h}, b_{l})'),
#		('DR_BHadWJa', '#Delta R (b_{h}, WJa)'),
#		('DR_BHadWJb', '#Delta R (b_{h} WJb)'),
#		('DR_BLepWJa', '#Delta R (b_{l}, WJa)'),
#		('DR_BLepWJb', '#Delta R (b_{l}, WJb)'),
#		('DR_WJaWJb', '#Delta R (WJa, WJb)'),
#		('DRmin_thad', 'min t_{h} #Delta R'),
#		('DRmin_tlep', 'min t_{l} #Delta R'),
#		('DRmax_thad', 'max t_{h} #Delta R')
#	]
#	
#	for var, xaxis in Gen_DR:
#	      hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
#	      mean = hist.GetMean()
#	      rms = hist.GetRMS()
#	
#	      plotter.set_histo_style(hist, color=defcol)
#	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
#	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
#	      plotter.save(var)
	
	Gen_Pt = [
#		('Pt_Lep', 'l p_{t}'),
#		('Pt_BHad', 'b_{h} p_{t}'),
#		('Pt_BLep', 'b_{l} p_{t}'),
#		('Pt_WJa', 'WJa p_{t}'),
#		('Pt_WJb', 'WJb p_{t}'),
		('Pt_ttbar', 't#bar t p_{t}')
#		('Pt_thad', 't_{h} p_{t}'),
#		('Pt_tlep', 't_{l} p_{t}')
	]
	
	for var, xaxis in Gen_Pt:
	      hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	      plotter.save(var)
	
#	Gen_Eta = [
#		('Eta_Lep', 'l #eta'),
#		('Eta_BHad', 'b_{h} #eta'),
#		('Eta_BLep', 'b_{l} #eta'),
#		('Eta_WJa', 'WJa #eta'),
#		('Eta_WJb', 'WJb #eta'),
#	]
#	
#	for var, xaxis in Gen_Eta:
#	      hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
#	      mean = hist.GetMean()
#	      rms = hist.GetRMS()
#	
#	      plotter.set_histo_style(hist, color=defcol)
#	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
#	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
#	      plotter.save(var)
#	
#	Gen_DR_2D = [
#		('DR_LepBHad_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (l, b_{h})'),
#		('DR_LepBLep_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (l, b_{l})'),
#		('DR_LepWJa_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (l, WJa)'),
#		('DR_LepWJb_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (l, WJb)'),
#		('DR_BHadBLep_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (b_{h}, b_{l})'),
#		('DR_BHadWJa_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (b_{h}, WJa)'),
#		('DR_BHadWJb_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (b_{h}, WJb)'),
#		('DR_BLepWJa_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (b_{l}, WJa)'),
#		('DR_BLepWJb_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (b_{l}, WJb)'),
#		('DR_WJaWJb_vs_Mtt', 'm_{t#bar t} [GeV]', '#Delta R (WJa, WJb)'),
#		('DRmin_thad_vs_mttbar', 'm_{t#bar t} [GeV]', 'min t_{h} #Delta R'),
#		('DRmin_tlep_vs_mttbar', 'm_{t#bar t} [GeV]', 'min t_{l} #Delta R'),
#		('DRmax_thad_vs_mttbar', 'm_{t#bar t} [GeV]', 'max t_{h} #Delta R'),
#		('DRmin_thad_vs_ptthad', 't_{h} p_{t}', 'min t_{h} #Delta R'),
#		('DRmin_tlep_vs_pttlep', 't_{l} p_{t}', 'min t_{l} #Delta R'),
#		('DRmax_thad_vs_ptthad', 't_{h} p_{t}', 'max t_{h} #Delta R')
#	]
#	
#	for var, xaxis,yaxis in Gen_DR_2D:
#		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
#	#	plotter.set_histo_style(hist, color=defcol)
#		plotter.plot(hist, xtitle=xaxis, ytitle=yaxis, drawstyle='hist')
#		hist.Draw('colz')
#		hist.xaxis.set_title(xaxis)
#		hist.yaxis.set_title(yaxis)
#		if var == 'DRmin_thad_vs_ptthad_hist' or var == 'DRmin_tlep_vs_pttlep_hist' or var == 'DRmax_thad_vs_ptthad_hist':
#			hist.xaxis.range_user = 0, 2000.
#		else:
#			hist.xaxis.range_user = mass_min, mass_max
#		plotter.save(var)
	
	Gen_Pt_2D = [
#		('Pt_Lep_vs_Mtt', 'm_{t#bar t} [GeV]', 'l p_{t}'),
#		('Pt_BHad_vs_Mtt', 'm_{t#bar t} [GeV]', 'b_{h} p_{t}'),
#		('Pt_BLep_vs_Mtt', 'm_{t#bar t} [GeV]', 'b_{l} p_{t}'),
#		('Pt_WJa_vs_Mtt', 'm_{t#bar t} [GeV]', 'WJa p_{t}'),
#		('Pt_WJb_vs_Mtt', 'm_{t#bar t} [GeV]', 'WJb p_{t}')
		('ptttbar_vs_mttbar', 'm_{t#bar t} [GeV]', 'p_{T} t#bar t')
	]
	
	for var, xaxis,yaxis in Gen_Pt_2D:
		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	#	plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
		hist.Draw('colz')
		hist.xaxis.set_title(xaxis)
		hist.yaxis.set_title(yaxis)
		plotter.save(var)
	
#	Gen_Eta_2D = [
#		('Eta_Lep_vs_Mtt', 'm_{t#bar t} [GeV]', 'l #eta'),
#		('Eta_BHad_vs_Mtt', 'm_{t#bar t} [GeV]', 'b_{h} #eta'),
#		('Eta_BLep_vs_Mtt', 'm_{t#bar t} [GeV]', 'b_{l} #eta'),
#		('Eta_WJa_vs_Mtt', 'm_{t#bar t} [GeV]', 'WJa #eta'),
#		('Eta_WJb_vs_Mtt', 'm_{t#bar t} [GeV]', 'WJb #eta')
#	]
#	
#	for var, xaxis, yaxis in Gen_Eta_2D:
#		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
#	#	plotter.set_histo_style(hist, color=defcol)
#		plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
#		hist.Draw('colz')
#		hist.xaxis.set_title(xaxis)
#		hist.yaxis.set_title(yaxis)
#		plotter.save(var)

	MassTTbar = [
		('Mass_ttbar', 'm_{t#bar t} [GeV]')
#		('Mass_thad', 'm_{t_{h}} [GeV]'),
#		('Mass_tlep', 'm_{t_{l}} [GeV]')
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

#	nJets = [
#		('nJets')
#	]
#
#	for var in nJets:
#		efficiencies = []
#		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
#		mean = hist.GetMean()
#		rms = hist.GetRMS()
#		total= hist.Integral()
#		plotter.set_histo_style(hist, color=defcol)
#		plotter.plot(hist, xtitle=var, ytitle=defyax, drawstyle='hist')
#		box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
#		plotter.save(var)
#		for i in range(0, hist.GetXaxis().GetNbins()):
#			efficiencies.append(format(hist.GetBinContent(i+1)/total, '.4f'))
#		print(var, efficiencies)


    ### Hadronic side events with all jets
	Gen_Hadronic_Events = [
		('Gen_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
		('Gen_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
		('Gen_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
		('Gen_Merged_THad_vs_mttbar', 'All 3 Merged', 'green'),
		('Gen_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
	]

	for name, title, col in DRvals:

		to_draw = []
		Add_Ratio_Hists = []
		Ratio_Hists = []
	
		for var, legends, colors in Gen_Hadronic_Events:
			hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
			plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
			to_draw.append(hist)
	
		total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]
	
		### plot log of event occurrences
		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
		plotter.save('Gen_Hadronic_Events_vs_mttbar_log_'+name)
	
	
		### plot event ratios
		NonReco_Ratio = to_draw[4]/total
		Ratio_Hists.append(NonReco_Ratio)
	
		Resolved_Ratio = to_draw[0]/total
		Ratio_Hists.append(Resolved_Ratio)
	
		M_BHadW_Ratio = to_draw[1]/total
		Ratio_Hists.append(M_BHadW_Ratio)
	
		M_WJet_Ratio = to_draw[2]/total
		Ratio_Hists.append(M_WJet_Ratio)
	
		M_THad_Ratio = to_draw[3]/total
		Ratio_Hists.append(M_THad_Ratio)
	
		plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.15), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
		plotter.save('Gen_Hadronic_Events_vs_mttbar_ratio_'+name)
	
	
		### create stacked hists of ratios
		stack_hists = []
		for i in range(len(to_draw)):
			stack_hists.append(to_draw[i].Integral())
		stack_hists.sort()
	
		for i in range(len(stack_hists)):
			for j in range(len(to_draw)):
				if to_draw[j].Integral() == stack_hists[i]:
					to_draw[j].SetFillStyle(1001)
					Add_Ratio_Hists.append(to_draw[j])
	
		stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total, Add_Ratio_Hists[4]/total)
		plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title)
		plotter.save('Gen_Hadronic_Events_vs_mttbar_stack_'+name)



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

if args.plot == "Had_Comp":
### compares relative frequencies of:
### 1.   three jets on hadronic side being resolved
### 2.   only bhad and one w jet merging
### 3.   only jets from w merging
### 4.   all three hadronic jets merging
### 5.   non-reconstructable events


    ### Events with 3 jets
	Hadronic_Events_3J = [
		('3J_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
		('3J_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
		('3J_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
		('3J_Merged_THad_vs_mttbar', 'All 3 Merged', 'green'),
		('3J_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
	]

	for name, title, col in DRvals:

		to_draw = []
		Add_Ratio_Hists = []
		Ratio_Hists = []
	
		for var, legends, colors in Hadronic_Events_3J:
			hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
			plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
			to_draw.append(hist)
	
		total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]
	
		### plot log of event occurrences
		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
		plotter.save('3J_Hadronic_Events_vs_mttbar_log'+name)
	
	
		### plot event ratios
		NonReco_Ratio = to_draw[4]/total
		Ratio_Hists.append(NonReco_Ratio)
	
		Resolved_Ratio = to_draw[0]/total
		Ratio_Hists.append(Resolved_Ratio)
	
		M_BHadW_Ratio = to_draw[1]/total
		Ratio_Hists.append(M_BHadW_Ratio)
	
		M_WJet_Ratio = to_draw[2]/total
		Ratio_Hists.append(M_WJet_Ratio)
	
		M_THad_Ratio = to_draw[3]/total
		Ratio_Hists.append(M_THad_Ratio)
	
		plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.15), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
		plotter.save('3J_Hadronic_Events_vs_mttbar_ratio'+name)
	
	
		### create stacked hists of ratios
		stack_hists = []
		for i in range(len(to_draw)):
			stack_hists.append(to_draw[i].Integral())
		stack_hists.sort()
	
		for i in range(len(stack_hists)):
			for j in range(len(to_draw)):
				if to_draw[j].Integral() == stack_hists[i]:
					to_draw[j].SetFillStyle(1001)
					Add_Ratio_Hists.append(to_draw[j])
	
		stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total, Add_Ratio_Hists[4]/total)
		plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title)
		plotter.save('3J_Hadronic_Events_vs_mttbar_stack'+name)



   ### Events with 4 jets
	Hadronic_Events_4J = [
		('4J_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
		('4J_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
		('4J_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
		('4J_Merged_THad_vs_mttbar', 'All 3 Merged', 'green'),
		('4J_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
	]

	for name, title, col in DRvals:

		to_draw = []
		Add_Ratio_Hists = []
		Ratio_Hists = []
	
		for var, legends, colors in Hadronic_Events_4J:
			hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
			plotter.set_histo_style(hist, title=legends, color=colors, drawstyle='hist', linestyle='solid')
			to_draw.append(hist)
	
		total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]
	
		### plot log of event occurrences
		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
		plotter.save('4J_Hadronic_Events_vs_mttbar_log'+name)
	
	
		### plot event ratios
		NonReco_Ratio = to_draw[4]/total
		Ratio_Hists.append(NonReco_Ratio)
	
		Resolved_Ratio = to_draw[0]/total
		Ratio_Hists.append(Resolved_Ratio)
	
		M_BHadW_Ratio = to_draw[1]/total
		Ratio_Hists.append(M_BHadW_Ratio)
	
		M_WJet_Ratio = to_draw[2]/total
		Ratio_Hists.append(M_WJet_Ratio)
	
		M_THad_Ratio = to_draw[3]/total
		Ratio_Hists.append(M_THad_Ratio)
	
		plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,0.85), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
		plotter.save('4J_Hadronic_Events_vs_mttbar_ratio'+name)
	
	
		### create stacked hists of ratios
		stack_hists = []
		for i in range(len(to_draw)):
			stack_hists.append(to_draw[i].Integral())
		stack_hists.sort()
	
		for i in range(len(stack_hists)):
			for j in range(len(to_draw)):
				if to_draw[j].Integral() == stack_hists[i]:
					to_draw[j].SetFillStyle(1001)
					Add_Ratio_Hists.append(to_draw[j])
	
		stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total, Add_Ratio_Hists[4]/total)
		plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0, 1.3), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
		plotter.save('4J_Hadronic_Events_vs_mttbar_stack'+name)



    ##### Events with 5 or more jets
	Hadronic_Events_5PJ = [
		('5PJ_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
		('5PJ_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
		('5PJ_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
		('5PJ_Merged_THad_vs_mttbar', 'All 3 Merged', 'green'),
		('5PJ_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
	]

	for name, title, col in DRvals:

		to_draw = []
		Add_Ratio_Hists = []
		Ratio_Hists = []

		for var, legends, colors in Hadronic_Events_5PJ:
			hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
			plotter.set_histo_style(hist, title=legends, color=colors, drawstyle='hist', linestyle='solid')
			to_draw.append(hist)

		total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]

		### plot log of event occurrences
		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
		plotter.save('5PJ_Hadronic_Events_vs_mttbar_log'+name)


		### plot event ratios
		NonReco_Ratio = to_draw[4]/total
		Ratio_Hists.append(NonReco_Ratio)

		Resolved_Ratio = to_draw[0]/total
		Ratio_Hists.append(Resolved_Ratio)

		M_BHadW_Ratio = to_draw[1]/total
		Ratio_Hists.append(M_BHadW_Ratio)

		M_WJet_Ratio = to_draw[2]/total
		Ratio_Hists.append(M_WJet_Ratio)

		M_THad_Ratio = to_draw[3]/total
		Ratio_Hists.append(M_THad_Ratio)

		plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,0.75), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
		plotter.save('5PJ_Hadronic_Events_vs_mttbar_ratio'+name)


		### create stacked hists of ratios
		stack_hists = []
		for i in range(len(to_draw)):
			stack_hists.append(to_draw[i].Integral())
		stack_hists.sort()

		for i in range(len(stack_hists)):
			for j in range(len(to_draw)):
				if to_draw[j].Integral() == stack_hists[i]:
					to_draw[j].SetFillStyle(1001)
					Add_Ratio_Hists.append(to_draw[j])

		stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total, Add_Ratio_Hists[4]/total)
		plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0, 1.3), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
		plotter.save('5PJ_Hadronic_Events_vs_mttbar_stack'+name)



##############################################################################################

if args.plot == "Delta_Plots":
	# Delta (Gen - Reco) Plots

	Sig_Background_mass = [
		('3J_Preselection_Wrong_Perms_Delta_mass', 'Preselection Wrong', 'red', 'solid'),
		('3J_Square_Cut_Wrong_Perms_Delta_mass', ' Signal Wrong', 'red', 'dashed'),
		('3J_Preselection_Correct_Perms_Delta_mass', 'Preselection Correct', 'blue', 'solid'),
		('3J_Square_Cut_Correct_Perms_Delta_mass', 'Signal Correct', 'blue', 'dashed')
	]

	to_draw = []

	for var, legends, colors, linestyles in Sig_Background_mass:
		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
		plotter.set_histo_style(hist, title=legends, color=colors, linestyle=linestyles)
		to_draw.append(hist)

	Sig_correct_mass = asrootpy(myfile.Get('Delta_Plots/3J_Square_Cut_Correct_Perms_Delta_mass')).Clone()
	Presel_correct_mass = asrootpy(myfile.Get('Delta_Plots/3J_Preselection_Correct_Perms_Delta_mass')).Clone()
	Sig_wrong_mass = asrootpy(myfile.Get('Delta_Plots/3J_Square_Cut_Wrong_Perms_Delta_mass')).Clone()
	Presel_wrong_mass = asrootpy(myfile.Get('Delta_Plots/3J_Preselection_Wrong_Perms_Delta_mass')).Clone()
	print 'Gen-Reco Mass:'
	print 'Correct Preselection Purity = ', format(Presel_correct_mass.Integral()/(Presel_correct_mass.Integral()+Presel_wrong_mass.Integral()), '.4f')
	print 'Wrong Preselection Purity = ', format(Presel_wrong_mass.Integral()/(Presel_correct_mass.Integral()+Presel_wrong_mass.Integral()), '.4f')
	print 'Correct Signal Purity = ',  format(Sig_correct_mass.Integral()/(Sig_correct_mass.Integral()+Sig_wrong_mass.Integral()), '.4f')
	print 'Wrong Signal Purity = ',  format(Sig_wrong_mass.Integral()/(Sig_correct_mass.Integral()+Sig_wrong_mass.Integral()), '.4f')
	print 'Correct Efficiency = ', format(Sig_correct_mass.Integral()/Presel_correct_mass.Integral(), '.4f')
	print 'Wrong Efficiency = ', format(Sig_wrong_mass.Integral()/Presel_wrong_mass.Integral(), '.4f')


	#### plot log
	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l',logy=True, x_range=(-1000, 200), y_range=(10**(-1),10**6), xtitle='#Delta m_{t_{h}} [GeV] (Gen - Reco)', drawstyle='hist')
	plotter.save('3J_Sig_Background_Delta_mass_log')
	
	#### plot normal
	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', x_range=(-1000, 200), xtitle='#Delta m_{t_{h}} [GeV] (Gen - Reco)', drawstyle='hist')
	plotter.save('3J_Sig_Background_Delta_mass')
	

	Sig_Background_costh = [
		('3J_Preselection_Wrong_Perms_Delta_costh', 'Preselection Wrong', 'red', 'solid'),
		('3J_Square_Cut_Wrong_Perms_Delta_costh', 'Signal Wrong', 'red', 'dashed'),
		('3J_Preselection_Correct_Perms_Delta_costh', 'Preselection Correct', 'blue', 'solid'),
		('3J_Square_Cut_Correct_Perms_Delta_costh', 'Signal Correct', 'blue', 'dashed')
	]

	to_draw = []

	for var, legends, colors, linestyles in Sig_Background_costh:
		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
		plotter.set_histo_style(hist, title=legends, color=colors, linestyle=linestyles)
		to_draw.append(hist)

	Sig_correct_costh = asrootpy(myfile.Get('Delta_Plots/3J_Square_Cut_Correct_Perms_Delta_costh')).Clone()
	Presel_correct_costh = asrootpy(myfile.Get('Delta_Plots/3J_Preselection_Correct_Perms_Delta_costh')).Clone()
	Sig_wrong_costh = asrootpy(myfile.Get('Delta_Plots/3J_Square_Cut_Wrong_Perms_Delta_costh')).Clone()
	Presel_wrong_costh = asrootpy(myfile.Get('Delta_Plots/3J_Preselection_Wrong_Perms_Delta_costh')).Clone()
	print 'Gen-Reco Cos(#theta):'
	print 'Correct Presel/ Total Preselection = ', format(Presel_correct_costh.Integral()/(Presel_correct_costh.Integral()+Presel_wrong_costh.Integral()), '.4f')
	print 'Wrong Presel/ Total Preselection = ', format(Presel_wrong_costh.Integral()/(Presel_correct_costh.Integral()+Presel_wrong_costh.Integral()), '.4f')
	print 'Correct Signal/Total Signal = ',  format(Sig_correct_costh.Integral()/(Sig_correct_costh.Integral()+Sig_wrong_costh.Integral()), '.4f')
	print 'Wrong Signal/Total Signal = ',  format(Sig_wrong_costh.Integral()/(Sig_correct_costh.Integral()+Sig_wrong_costh.Integral()), '.4f')

	### plot log
	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', logy=True, x_range=(-2, 2), y_range=(10**(-1),10**7), xtitle='#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', drawstyle='hist')
	plotter.save('3J_Sig_Background_Delta_costh_log')
	
	### plot normal 
	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', x_range=(-2, 2), xtitle='#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', drawstyle='hist')
	plotter.save('3J_Sig_Background_Delta_costh')
	
		

	
#	CosTheta = [
#	      ('3J_LDA_less_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)'),
#	      ('3J_LDA_great_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)'),
#	      ('3J_Mass_Cut_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)'),
#	      ('3J_Square_Cut_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)'),
#	      ('3J_No_Cut_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)'),
#	      ('4J_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)')
#	]
#	
#	for var, xaxis in CosTheta:
#	      hist = asrootpy(normfile.Get('Delta_Plots/'+var)).Clone()
#	      mean = hist.GetMean()
#	      rms = hist.GetRMS()
#	
#	      plotter.set_histo_style(hist, color=defcol)
#	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
#	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
#	      plotter.save(var)
#	
#	Matched_CosTheta = [
#	      ('3J_LDA_less_Matched_BHad_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{h}'),
#	      ('3J_LDA_less_Matched_BLep_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{l}'),
#	      ('3J_LDA_great_Matched_BHad_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{h}'),
#	      ('3J_LDA_great_Matched_BLep_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{l}'),
#	      ('3J_Mass_Cut_Matched_BHad_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{h}'),
#	      ('3J_Mass_Cut_Matched_BLep_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{l}'),
#	      ('3J_Square_Cut_Matched_BHad_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{h}'),
#	      ('3J_Square_Cut_Matched_BLep_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{l}'),
#	      ('3J_No_Cut_Matched_BHad_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{h}'),
#	      ('3J_No_Cut_Matched_BLep_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{l}'),
#	      ('4J_Matched_BHad_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{h}'),
#	      ('4J_Matched_BLep_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', 'b_{l}')
#	]
#	
#	for var, xaxis, matched in Matched_CosTheta:
#		hist = asrootpy(normfile.Get('Delta_Plots/'+var)).Clone()
#		mean = hist.GetMean()
#		rms = hist.GetRMS()
#		
#		plotter.set_histo_style(hist, color=defcol)
#		if hist.Integral() != 0:
#			plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
#		box = plotter.make_text_box('Mean = %f\nRMS = %f\nMatched to %s' % (mean,rms,matched), position='NE')
#		plotter.save(var)
#
#	THad_Mass = [
#	      ('3J_LDA_less_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)'),
#	      ('3J_LDA_great_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)'),
#	      ('3J_Mass_Cut_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)'),
#	      ('3J_Square_Cut_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)'),
#	      ('3J_No_Cut_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)'),
#	      ('4J_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)')
#	]
#	
#	for var, xaxis in THad_Mass:
#	      hist = asrootpy(normfile.Get('Delta_Plots/'+var)).Clone()
#	      mean = hist.GetMean()
#	      rms = hist.GetRMS()
#	
#	      hist.xaxis.range_user = -1500, 500
#	      plotter.set_histo_style(hist, color=defcol)
#	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
#	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
#	      plotter.save(var)
#
#	Matched_THad_Mass = [
#	      ('3J_LDA_less_Matched_BHad_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{h}'),
#	      ('3J_LDA_less_Matched_BLep_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{l}'),
#	      ('3J_LDA_great_Matched_BHad_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{h}'),
#	      ('3J_LDA_great_Matched_BLep_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{l}'),
#	      ('3J_Mass_Cut_Matched_BHad_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{h}'),
#	      ('3J_Mass_Cut_Matched_BLep_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{l}'),
#	      ('3J_Square_Cut_Matched_BHad_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{h}'),
#	      ('3J_Square_Cut_Matched_BLep_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{l}'),
#	      ('3J_No_Cut_Matched_BHad_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{h}'),
#	      ('3J_No_Cut_Matched_BLep_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{l}'),
#	      ('4J_Matched_BHad_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{h}'),
#	      ('4J_Matched_BLep_THad_Delta_mass', '#Delta m_{t_{h}} [GeV] (Gen - Reco)', 'b_{l}')
#	]
#	
#	for var, xaxis, matched in Matched_THad_Mass:
#		hist = asrootpy(normfile.Get('Delta_Plots/'+var)).Clone()
#		mean = hist.GetMean()
#		rms = hist.GetRMS()
#		
#		hist.xaxis.range_user = -1500, 500
#		plotter.set_histo_style(hist, color=defcol)
#		if hist.Integral() != 0:		
#			plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
#		box = plotter.make_text_box('Mean = %f\nRMS = %f\nMatched to %s' % (mean,rms,matched), position='NE')
#		plotter.save(var)
#
#	Jet3_LDA_No_Match_Comparison_THad_Mass = [
#		('3J_LDA_less_THad_Delta_mass', 'No Matching <', 'red'),
#		('3J_LDA_great_THad_Delta_mass', 'No Matching >', 'blue')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	to_draw = []
#
#	for var, label, colors in Jet3_LDA_No_Match_Comparison_THad_Mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=label, color=colors)
#		to_draw.append(hist)
#
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#	plotter.save('3J_Comb_LDA_No_Matching_THad_Delta_mass')
#
#	Jet3_LDA_Matched_BHad_Comparison_THad_Mass = [
#		('3J_LDA_less_Matched_BHad_THad_Delta_mass', 'Matched to b_{h} <', 'red'),
#		('3J_LDA_great_Matched_BHad_THad_Delta_mass', 'Matched to b_{h} >', 'blue')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	to_draw = []
#
#	for var, label, colors in Jet3_LDA_Matched_BHad_Comparison_THad_Mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=label, color=colors)
#		to_draw.append(hist)
#
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#	plotter.save('3J_Comb_LDA_Matched_BHad_THad_Delta_mass')
#
#	Jet3_LDA_Matched_BLep_Comparison_THad_Mass = [
#		('3J_LDA_less_Matched_BLep_THad_Delta_mass', 'Matched to b_{l} <', 'red'),
#		('3J_LDA_great_Matched_BLep_THad_Delta_mass', 'Matched to b_{l} >', 'blue')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	to_draw = []
#
#	for var, label, colors in Jet3_LDA_Matched_BLep_Comparison_THad_Mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=label, color=colors)
#		to_draw.append(hist)
#
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#	plotter.save('3J_Comb_LDA_Matched_BLep_THad_Delta_mass')
#
#	Jet3_All_LDA_THad_Mass = [
#	      ('3J_LDA_less_THad_Delta_mass', 'No Matching <', 'black'),
#	      ('3J_LDA_less_Matched_BHad_THad_Delta_mass', 'Matched to b_{h} <', 'red'),
#	      ('3J_LDA_less_Matched_BLep_THad_Delta_mass', 'Matched to b_{l} <', 'blue'),
#	      ('3J_LDA_great_THad_Delta_mass', 'No Matching >', 'green'),
#	      ('3J_LDA_great_Matched_BHad_THad_Delta_mass', 'Matched to b_{h} >', 'cyan'),
#	      ('3J_LDA_great_Matched_BLep_THad_Delta_mass', 'Matched to b_{l} >', 'magenta')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	to_draw = []
#
#	for var, matched, colors in Jet3_All_LDA_THad_Mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=matched, color=colors)
#		to_draw.append(hist)
#
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#	plotter.save('3J_Comb_LDA_all_THad_Delta_mass')
#
#	Jet3_LDA_Less_THad_Mass = [
#	      ('3J_LDA_less_THad_Delta_mass', 'No Matching <', 'black'),
#	      ('3J_LDA_less_Matched_BHad_THad_Delta_mass', 'Matched to b_{h} <', 'red'),
#	      ('3J_LDA_less_Matched_BLep_THad_Delta_mass', 'Matched to b_{l} <', 'blue')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	to_draw = []
#
#	for var, matched, colors in Jet3_LDA_Less_THad_Mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=matched, color=colors)
#		to_draw.append(hist)
#
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#	plotter.save('3J_Comb_LDA_less_THad_Delta_mass')
#
#	Jet3_LDA_great_THad_Mass = [
#	      ('3J_LDA_great_THad_Delta_mass', 'No Matching >', 'black'),
#	      ('3J_LDA_great_Matched_BHad_THad_Delta_mass', 'Matched to b_{h} >', 'red'),
#	      ('3J_LDA_great_Matched_BLep_THad_Delta_mass', 'Matched to b_{l} >', 'blue')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	to_draw = []
#
#	for var, matched, colors in Jet3_LDA_great_THad_Mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=matched, color=colors)
#		to_draw.append(hist)
#
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#	plotter.save('3J_Comb_LDA_great_THad_Delta_mass')
#
#	Jet3_No_Cut_THad_Mass = [
#	      ('3J_No_Cut_THad_Delta_mass', 'No Matching Restrictions', 'black'),
#	      ('3J_No_Cut_Matched_BHad_THad_Delta_mass', 'Matched to b_{h}', 'red'),
#	      ('3J_No_Cut_Matched_BLep_THad_Delta_mass', 'Matched to b_{l}', 'blue')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	to_draw = []
#
#	for var, matched, colors in Jet3_No_Cut_THad_Mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=matched, color=colors)
#		to_draw.append(hist)
#
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#	plotter.save('3J_Comb_No_Cut_THad_Delta_mass')
#
#	Jet3_Mass_Cut_THad_Mass = [
#	      ('3J_Mass_Cut_THad_Delta_mass', 'No Matching Restrictions', 'black'),
#	      ('3J_Mass_Cut_Matched_BHad_THad_Delta_mass', 'Matched to b_{h}', 'red'),
#	      ('3J_Mass_Cut_Matched_BLep_THad_Delta_mass', 'Matched to b_{l}', 'blue')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	to_draw = []
#
#	for var, matched, colors in Jet3_Mass_Cut_THad_Mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=matched, color=colors)
#		to_draw.append(hist)
#
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#	plotter.save('3J_Comb_Mass_Cut_THad_Delta_mass')
#
#	Jet3_Square_Cut_THad_Mass = [
#	      ('3J_Square_Cut_THad_Delta_mass', 'No Matching Restrictions', 'black'),
#	      ('3J_Square_Cut_Matched_BHad_THad_Delta_mass', 'Matched to b_{h}', 'red'),
#	      ('3J__Cut_Matched_BLep_THad_Delta_mass', 'Matched to b_{l}', 'blue')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	to_draw = []
#
#	for var, matched, colors in Jet3_Mass_Cut_THad_Mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=matched, color=colors)
#		to_draw.append(hist)
#
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#	box = plotter.make_text_box('m_{btag} > 50 GeV\nm_{btag+jet} < 180 GeV', position='NW')
#	plotter.save('3J_Comb_Square_Cut_THad_Delta_mass')
#
#	Jet3_Cut_No_Cut_THad_Mass = [
#		('3J_No_Cut_THad_Delta_mass', '3J_LDA_less_THad_Delta_mass', '3J_LDA_great_THad_Delta_mass', '3J_Mass_Cut_THad_Delta_mass', '3J_Square_Cut_THad_Delta_mass', 'No Matching', '3J_Comb_No_Matching_THad_Delta_mass'),
#		('3J_No_Cut_Matched_BHad_THad_Delta_mass', '3J_LDA_less_Matched_BHad_THad_Delta_mass', '3J_LDA_great_Matched_BHad_THad_Delta_mass', '3J_Mass_Cut_Matched_BHad_THad_Delta_mass', '3J_Square_Cut_Matched_BHad_THad_Delta_mass', 'Matched to b_{h}', '3J_Comb_Matched_BHad_THad_Delta_mass'),
#		('3J_No_Cut_Matched_BLep_THad_Delta_mass', '3J_LDA_less_Matched_BLep_THad_Delta_mass', '3J_LDA_great_Matched_BLep_THad_Delta_mass', '3J_Mass_Cut_Matched_BLep_THad_Delta_mass', '3J_Square_Cut_Matched_BLep_THad_Delta_mass', 'Matched to b_{l}', '3J_Comb_Matched_BLep_THad_Delta_mass')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	for nc, lda_less, lda_great, mc, square, titles, name in Jet3_Cut_No_Cut_THad_Mass:
#		to_draw = []
#
#		No_Cut = asrootpy(myfile.Get('Delta_Plots/'+nc)).Clone()
#		plotter.set_histo_style(No_Cut, title='No Cut', color='black')
#		to_draw.append(No_Cut)
#
#		LDA_less = asrootpy(myfile.Get('Delta_Plots/'+lda_less)).Clone()
#		plotter.set_histo_style(LDA_less, title='LDA <', color='red')
#		to_draw.append(LDA_less)
#
#		LDA_great = asrootpy(myfile.Get('Delta_Plots/'+lda_great)).Clone()
#		plotter.set_histo_style(LDA_great, title='LDA >', color='cyan')
#		to_draw.append(LDA_great)
#
#		Mass_Cut = asrootpy(myfile.Get('Delta_Plots/'+mc)).Clone()
#		plotter.set_histo_style(Mass_Cut, title='Mass Cut', color='blue')
#		to_draw.append(Mass_Cut)
#
#		Square_Cut = asrootpy(myfile.Get('Delta_Plots/'+square)).Clone()
#		plotter.set_histo_style(Square_Cut, title='Square Cut', color='magenta')
#		to_draw.append(Square_Cut)
#
#		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#		box = plotter.make_text_box(titles, position='NW')
#		plotter.save(name)
#
#
#	Jet4_THad_Mass = [
#	      ('4J_THad_Delta_mass', 'No Matching Restrictions', 'black'),
#	      ('4J_Matched_BHad_THad_Delta_mass', 'Matched to b_{h}', 'red'),
#	      ('4J_Matched_BLep_THad_Delta_mass', 'Matched to b_{l}', 'blue')
#	]
#
#	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'
#
#	to_draw = []
#
#	for var, matched, colors in Jet4_THad_Mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=matched, color=colors)
#		to_draw.append(hist)
#
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1500, 500), xtitle=xaxis, drawstyle='hist')
#	plotter.save('4J_Comb_THad_Delta_mass')


###### merged best_jets

	Merged_BHad_THad_Mass = [
		('3J_Merged_BHad_Delta_THad_mass', 'b_{h}'),
		('4J_Merged_BHad_Delta_THad_mass', 'b_{h}')
	]

	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'

	for var, merged_b in Merged_BHad_THad_Mass:
	      hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      hist.xaxis.range_user = -200, 200
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
	      box = plotter.make_text_box('Mean = %f\nRMS = %f\nMerged %s' % (mean,rms,merged_b), position='NE')
	      plotter.save(var)

	Merged_BLep_THad_Mass = [
		('3J_Merged_BLep_Delta_THad_mass', 'b_{l}'),
		('4J_Merged_BLep_Delta_THad_mass', 'b_{l}')
	]

	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'

	for var, merged_b in Merged_BLep_THad_Mass:
	      hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      hist.xaxis.range_user = -1000, 200
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
	      box = plotter.make_text_box('Mean = %f\nRMS = %f\nMerged %s' % (mean,rms,merged_b), position='NE')
	      plotter.save(var)


	Comb_Merged_BHad_THad_Mass = [
		('3J_Merged_BHad_Delta_THad_mass', 'njets = 3', 'black'),
		('4J_Merged_BHad_Delta_THad_mass', 'njets = 4', 'red')
	]

	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'

	to_draw = []

	for var, njets, colors in Comb_Merged_BHad_THad_Mass:
		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
		plotter.set_histo_style(hist, title=njets, color=colors)
		to_draw.append(hist)

	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-200, 200), xtitle=xaxis, drawstyle='hist')
	box = plotter.make_text_box('Merged b_{h}', position='NW')
	plotter.save('Comb_Merged_BHad_THad_Delta_mass')

	Comb_Merged_BLep_THad_Mass = [
		('3J_Merged_BLep_Delta_THad_mass', 'njets = 3', 'black'),
		('4J_Merged_BLep_Delta_THad_mass', 'njets = 4', 'red')
	]

	xaxis = '#Delta m_{t_{h}} [GeV] (Gen - Reco)'

	to_draw = []

	for var, njets, colors in Comb_Merged_BLep_THad_Mass:
		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
		plotter.set_histo_style(hist, title=njets, color=colors)
		to_draw.append(hist)

	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(-1000, 200), xtitle=xaxis, drawstyle='hist')
	box = plotter.make_text_box('Merged b_{l}', position='NW')
	plotter.save('Comb_Merged_BLep_THad_Delta_mass')


##############################################
print('~/nobackup/CMSSW_8_0_7_patch2/src/Analyses/URTTbar/htt_scripts/plots/jet_reco/%s/%s/%s/%s/' % (jobid, args.analysis, args.plot, args.sample))
