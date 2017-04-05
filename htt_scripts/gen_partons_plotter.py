'''
Gen_Partons Analyzer Plotter macro
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

parser = argparse.ArgumentParser(description='Create plots using files from gen_partons.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full_Analysis).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, or Combined).')
parser.add_argument('plot', help='Choose type of plots to generate (DR, Pt, Eta, Costh, Mass, nJets, Ratios, Had_Comp).')
args = parser.parse_args()

if args.analysis == "Test":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../%s.gen_partons.test.root' % args.sample, 'read')
		normfile = views.NormalizeView(root_open('../%s.gen_partons.test.root' % args.sample, 'read'))#normalized file
		plotter = BasePlotter(
			'plots/gen_partons/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Sample_ttJetsM0 = root_open('../ttJetsM0.gen_partons.test.root', 'read')
	#	Sample_ttJetsM700 = root_open('../ttJetsM700.gen_partons.test.root', 'read')
	#	Sample_ttJetsM1000 = root_open('../ttJetsM1000.gen_partons.test.root', 'read')
		plotter = BasePlotter(
			'plots/gen_partons/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)

elif args.analysis == "Full_Analysis":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../results/%s/gen_partons/%s.root' % (jobid, args.sample), 'read')
		normfile = views.NormalizeView(root_open('../results/%s/gen_partons/%s.root' % (jobid, args.sample), 'read'))
		plotter = BasePlotter(
			'plots/gen_partons/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Full_ttJetsM0 = root_open('../results/%s/gen_partons/ttJets.root' % jobid, 'read')
	#	Full_ttJetsM700 = root_open('../results/%s/gen_partons/ttJetsM700.root' % jobid, 'read')
	#	Full_ttJetsM1000 = root_open('../results/%s/gen_partons/ttJetsM1000.root' % jobid, 'read')
		plotter = BasePlotter(
			'plots/gen_partons/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)


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
		('DRP8', '#Delta R < 0.8', 'black')
	}
#DRvals = {
#		('DRP4', '#Delta R < 0.4', 'red'),
#		('DRP5', '#Delta R < 0.5', 'blue'),
#		('DRP6', '#Delta R < 0.6', 'green'),
#		('DRP8', '#Delta R < 0.8', 'black')
#	}

#lpos = 'NE' #legend position
#types = ['dtype', 'utype'] #variable is down- or up-type jet
#colors = ['red', 'blue'] #plot legend colors (down is red, up is blue)
#labels = ['down', 'up'] #legend variable labels
#####

########################################################################################
if args.plot == "DR":
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
	
########################################################################################
if args.plot == "Pt":
	# Gen Object Plots
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
	
	Gen_Pt_2D = [
		('Pt_Lep_vs_Mtt', 'm_{t#bar t} [GeV]', 'l p_{T}'),
		('Pt_BHad_vs_Mtt', 'm_{t#bar t} [GeV]', 'b_{h} p_{T}'),
		('Pt_BLep_vs_Mtt', 'm_{t#bar t} [GeV]', 'b_{l} p_{T}'),
		('Pt_WJa_vs_Mtt', 'm_{t#bar t} [GeV]', 'WJa p_{T}'),
		('Pt_WJb_vs_Mtt', 'm_{t#bar t} [GeV]', 'WJb p_{T}'),
		('Pt_thad_vs_Mtt', 'm_{t#bar t} [GeV]', 't_{h} p_{T}')
	]
	
	for var, xaxis,yaxis in Gen_Pt_2D:
		hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	#	plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
		hist.Draw('colz')
		hist.xaxis.set_title(xaxis)
		hist.yaxis.set_title(yaxis)
		plotter.save(var)
	
########################################################################################
if args.plot == "Eta":
	# Gen Object Plots
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

	
	
########################################################################################
if args.plot == "Costh":
	# Gen Object Plots
	Gen_Costh = [
		('Costh_Lep', 'l Cos #theta'),
		('Costh_BLep', 'b_{l} Cos #theta'),
		('Costh_BHad', 'b_{h} Cos #theta'),
		('Costh_WJa', 'WJa Cos #theta'),
		('Costh_WJb', 'WJb Cos #theta'),
		('Costh_ttbar', 'WJb Cos #theta'),
	]
	
	for var, xaxis in Gen_Costh:
	      hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	      plotter.save(var)
	
########################################################################################
if args.plot == "Mass":
	# Gen Object Plots
	MassTTbar = [
		('Mass_Lep', 'm_{l} [GeV]'),
		('Mass_BLep', 'm_{b_{l}} [GeV]'),
		('Mass_BHad', 'm_{b_{h}} [GeV]'),
		('Mass_WJa', 'm_{WJa}} [GeV]'),
		('Mass_WJb', 'm_{WJb}} [GeV]'),
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

	
########################################################################################
if args.plot == "nJets":
	# Gen Object Plots
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
if args.plot == "Ratios":
	DRmin_ratios = [
#		('DRmin_thad_lDRP4_vs_mttbar', 'DRmin_thad_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', 't_{h} min #Delta R', 'DRmin_thad_lp4_ratio_vs_mttbar'),
#		('DRmin_thad_lDRP4_vs_ptthad', 'DRmin_thad_gDRP4_vs_ptthad', 't_{h} p_{t}', '#Delta R < 0.4 Fraction', 't_{h} min #Delta R', 'DRmin_thad_lp4_ratio_vs_ptthad'),
#		('DRmin_tlep_lDRP4_vs_mttbar', 'DRmin_tlep_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', 't_{l} min #Delta R', 'DRmin_tlep_lp4_ratio_vs_mttbar'),
#		('DRmin_tlep_lDRP4_vs_pttlep', 'DRmin_tlep_gDRP4_vs_pttlep', 't_{l} p_{t}', '#Delta R < 0.4 Fraction', 't_{l} min #Delta R', 'DRmin_tlep_lp4_ratio_vs_pttlep'),
		
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
#		('DRmax_thad_lDRP4_vs_mttbar', 'DRmax_thad_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', 't_{h} max #Delta R', 'DRmax_thad_lp4_ratio_vs_mttbar'),
#		('DRmax_thad_lDRP4_vs_ptthad', 'DRmax_thad_gDRP4_vs_ptthad', 't_{h} p_{t}', '#Delta R < 0.4 Fraction', 't_{h} max #Delta R', 'DRmax_thad_lp4_ratio_vs_ptthad')
	]

	for lp4, gp4, xaxis, yaxis, title, name in DRmin_ratios:# < 0.4, > 0.4, xaxis, yaxis, png name
		efficiencies = []
		
		LP4 = asrootpy(myfile.Get('Gen_Plots/'+lp4)).Clone()
		if LP4.Integral() != 0:
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


##################################################################################################
if args.plot == "Had_Comp":
### compares relative frequencies of:
### 1.   three jets on hadronic side being resolved
### 2.   only bhad and one w jet merging
### 3.   only jets from w merging
### 4.   all three hadronic jets merging
### 5.   non-reconstructable events


    ### Events with all jets
        Hadronic_Events = [
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

                for var, legends, colors in Hadronic_Events:
                        hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
                        plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
                        to_draw.append(hist)

                total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]

                ### plot log of event occurrences
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
                plotter.save('AllJets_Hadronic_Events_vs_mttbar_log_'+name)


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
                plotter.save('AllJets_Hadronic_Events_vs_mttbar_ratio_'+name)


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
                plotter.save('AllJets_Hadronic_Events_vs_mttbar_stack_'+name)



#    ### Events with 3 jets
#        Hadronic_Events_3J = [
#                ('3J_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
#                ('3J_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
#                ('3J_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
#                ('3J_Merged_THad_vs_mttbar', 'All 3 Merged', 'green'),
#                ('3J_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
#        ]
#
#        for name, title, col in DRvals:
#
#                to_draw = []
#                Add_Ratio_Hists = []
#                Ratio_Hists = []
#
#                for var, legends, colors in Hadronic_Events_3J:
#                        hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
#                        plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
#                        to_draw.append(hist)
#
#                total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]
#
#                ### plot log of event occurrences
#                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
#                plotter.save('3J_Hadronic_Events_vs_mttbar_log_'+name)
#
#
#                ### plot event ratios
#                NonReco_Ratio = to_draw[4]/total
#                Ratio_Hists.append(NonReco_Ratio)
#
#                Resolved_Ratio = to_draw[0]/total
#                Ratio_Hists.append(Resolved_Ratio)
#
#                M_BHadW_Ratio = to_draw[1]/total
#                Ratio_Hists.append(M_BHadW_Ratio)
#
#                M_WJet_Ratio = to_draw[2]/total
#                Ratio_Hists.append(M_WJet_Ratio)
#
#                M_THad_Ratio = to_draw[3]/total
#                Ratio_Hists.append(M_THad_Ratio)
#
#                plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.15), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
#                plotter.save('3J_Hadronic_Events_vs_mttbar_ratio_'+name)
#
#
#                ### create stacked hists of ratios
#                stack_hists = []
#                for i in range(len(to_draw)):
#                        stack_hists.append(to_draw[i].Integral())
#                stack_hists.sort()
#
#                for i in range(len(stack_hists)):
#                        for j in range(len(to_draw)):
#                                if to_draw[j].Integral() == stack_hists[i]:
#                                        to_draw[j].SetFillStyle(1001)
#                                        Add_Ratio_Hists.append(to_draw[j])
#
#                stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total, Add_Ratio_Hists[4]/total)
#                plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title)
#                plotter.save('3J_Hadronic_Events_vs_mttbar_stack_'+name)
#
#
#
#   ### Events with 4 jets
#        Hadronic_Events_4J = [
#                ('4J_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
#                ('4J_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
#                ('4J_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
#                ('4J_Merged_THad_vs_mttbar', 'All 3 Merged', 'green'),
#                ('4J_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
#        ]
#
#        for name, title, col in DRvals:
#
#                to_draw = []
#                Add_Ratio_Hists = []
#                Ratio_Hists = []
#
#                for var, legends, colors in Hadronic_Events_4J:
#                        hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
#                        plotter.set_histo_style(hist, title=legends, color=colors, drawstyle='hist', linestyle='solid')
#                        to_draw.append(hist)
#
#                total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]
#
#                ### plot log of event occurrences
#                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
#                plotter.save('4J_Hadronic_Events_vs_mttbar_log_'+name)
#
#
#                ### plot event ratios
#                NonReco_Ratio = to_draw[4]/total
#                Ratio_Hists.append(NonReco_Ratio)
#
#                Resolved_Ratio = to_draw[0]/total
#                Ratio_Hists.append(Resolved_Ratio)
#
#                M_BHadW_Ratio = to_draw[1]/total
#                Ratio_Hists.append(M_BHadW_Ratio)
#
#                M_WJet_Ratio = to_draw[2]/total
#                Ratio_Hists.append(M_WJet_Ratio)
#
#                M_THad_Ratio = to_draw[3]/total
#                Ratio_Hists.append(M_THad_Ratio)
#
#                plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,0.85), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
#                plotter.save('4J_Hadronic_Events_vs_mttbar_ratio_'+name)
#
#
#                ### create stacked hists of ratios
#                stack_hists = []
#                for i in range(len(to_draw)):
#                        stack_hists.append(to_draw[i].Integral())
#                stack_hists.sort()
#
#                for i in range(len(stack_hists)):
#                        for j in range(len(to_draw)):
#                                if to_draw[j].Integral() == stack_hists[i]:
#                                        to_draw[j].SetFillStyle(1001)
#                                        Add_Ratio_Hists.append(to_draw[j])
#
#                stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total, Add_Ratio_Hists[4]/total)
#                plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0, 1.3), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
#                plotter.save('4J_Hadronic_Events_vs_mttbar_stack_'+name)
#
#
#
#    ##### Events with 5 or more jets
#        Hadronic_Events_5PJ = [
#                ('5PJ_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
#                ('5PJ_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
#                ('5PJ_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
#                ('5PJ_Merged_THad_vs_mttbar', 'All 3 Merged', 'green'),
#                ('5PJ_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
#        ]
#
#        for name, title, col in DRvals:
#
#                to_draw = []
#                Add_Ratio_Hists = []
#                Ratio_Hists = []
#
#                for var, legends, colors in Hadronic_Events_5PJ:
#                        hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
#                        plotter.set_histo_style(hist, title=legends, color=colors, drawstyle='hist', linestyle='solid')
#                        to_draw.append(hist)
#
#                total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]
#
#                ### plot log of event occurrences
#                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
#                plotter.save('5PJ_Hadronic_Events_vs_mttbar_log_'+name)
#
#
#                ### plot event ratios
#                NonReco_Ratio = to_draw[4]/total
#                Ratio_Hists.append(NonReco_Ratio)
#
#                Resolved_Ratio = to_draw[0]/total
#                Ratio_Hists.append(Resolved_Ratio)
#
#                M_BHadW_Ratio = to_draw[1]/total
#                Ratio_Hists.append(M_BHadW_Ratio)
#
#                M_WJet_Ratio = to_draw[2]/total
#                Ratio_Hists.append(M_WJet_Ratio)
#
#                M_THad_Ratio = to_draw[3]/total
#                Ratio_Hists.append(M_THad_Ratio)
#
#                plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,0.75), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
#                plotter.save('5PJ_Hadronic_Events_vs_mttbar_ratio_'+name)
#
#
#                ### create stacked hists of ratios
#                stack_hists = []
#                for i in range(len(to_draw)):
#                        stack_hists.append(to_draw[i].Integral())
#                stack_hists.sort()
#
#                for i in range(len(stack_hists)):
#                        for j in range(len(to_draw)):
#                                if to_draw[j].Integral() == stack_hists[i]:
#                                        to_draw[j].SetFillStyle(1001)
#                                        Add_Ratio_Hists.append(to_draw[j])
#
#                stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total, Add_Ratio_Hists[4]/total)
#                plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0, 1.3), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
#                plotter.save('5PJ_Hadronic_Events_vs_mttbar_stack_'+name)
#

#############################################################

print('~/nobackup/CMSSW_8_0_7_patch2/src/Analyses/URTTbar/htt_scripts/plots/gen_partons/%s/%s/%s/%s/' % (jobid, args.analysis, args.plot, args.sample))
