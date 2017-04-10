'''
Test Analyzer Plotter macro
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

parser = argparse.ArgumentParser(description='Create plots using files from Test_Analyzer.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full_Analysis).')
parser.add_argument('sample', help='Choose a file (ttJets, ttJetsM700, ttJetsM1000, or Combined).')
parser.add_argument('plot', help='Choose type of plots to generate (Gen_Plots, Delta_Plots).')
args = parser.parse_args()

if args.analysis == "Test":
	if args.sample == "ttJets" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../%s.Test_Analyzer.test.root' % args.sample, 'read')
		plotter = BasePlotter(
			'plots/test_analyzer/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Sample_ttJets = root_open('../ttJets.Test_Analyzer.test.root', 'read')
	#	Sample_ttJetsM700 = root_open('../ttJetsM700.Test_Analyzer.test.root', 'read')
	#	Sample_ttJetsM1000 = root_open('../ttJetsM1000.Test_Analyzer.test.root', 'read')
		plotter = BasePlotter(
			'plots/test_analyzer/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)

elif args.analysis == "Full_Analysis":
	if args.sample == "ttJets" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../results/%s/Test_Analyzer/%s.root' % (jobid, args.sample), 'read')
		plotter = BasePlotter(
			'plots/test_analyzer/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Full_ttJets = root_open('../results/%s/Test_Analyzer/ttJets.root' % jobid, 'read')
	#	Full_ttJetsM700 = root_open('../results/%s/Test_Analyzer/ttJetsM700.root' % jobid, 'read')
	#	Full_ttJetsM1000 = root_open('../results/%s/Test_Analyzer/ttJetsM1000.root' % jobid, 'read')
		plotter = BasePlotter(
			'plots/test_analyzer/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)


##### Comparing right/wrong combinations from TTBarSolver settings ####
##open files
#NS_MD_file = root_open('../ttJets.Test_Analyzer_NS_MassDisc.test.root', 'read') #only NS and MassDis right used in solver
#NS_MD_ratio_file = root_open('../ttJets.Test_Analyzer_NS_MassDisc_Ratio.test.root', 'read') #NS and MassDis right/wrong used in solver
#NS_MD_Ang_file = root_open('../ttJets.Test_Analyzer_NS_MassDisc_Ang.test.root', 'read') #NS, MassDisc, Ang vars right used in solver
#NS_MD_Ang_ratio_file = root_open('../ttJets.Test_Analyzer_NS_MassDisc_Ang_Ratio.test.root', 'read') #NS, MassDisc, Ang vars right/wrong used in solver
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
#lpos = 'NE' #legend position
#types = ['dtype', 'utype'] #variable is down- or up-type jet
#colors = ['red', 'blue'] #plot legend colors (down is red, up is blue)
#labels = ['down', 'up'] #legend variable labels
#####

if args.plot == "Gen_Plots":
	# Gen Object Plots
	Gen_DR = [
		('DR_LepBHad_hist', '#Delta R (l, b_{h})', 'DR_LepBHad_hist'),
		('DR_LepBLep_hist', '#Delta R (l, b_{l})', 'DR_LepBLep_hist'),
		('DR_LepWJa_hist', '#Delta R (l, WJa)', 'DR_LepWJa_hist'),
		('DR_LepWJb_hist', '#Delta R (l, WJb)', 'DR_LepWJb_hist'),
		('DR_BHadBLep_hist', '#Delta R (b_{h}, b_{l})', 'DR_BHadBLep_hist'),
		('DR_BHadWJa_hist', '#Delta R (b_{h}, WJa)', 'DR_BHadWJa_hist'),
		('DR_BHadWJb_hist', '#Delta R (b_{h} WJb)', 'DR_BHadWJb_hist'),
		('DR_BLepWJa_hist', '#Delta R (b_{l}, WJa)', 'DR_BLepWJa_hist'),
		('DR_BLepWJb_hist', '#Delta R (b_{l}, WJb)', 'DR_BLepWJb_hist'),
		('DR_WJaWJb_hist', '#Delta R (WJa, WJb)', 'DR_WJaWJb_hist')
	]
	
	for var, xaxis, names in Gen_DR:
	      hist = asrootpy(myfile.Get(var)).Clone()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax)
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	      plotter.save(names)
	
	Gen_Pt = [
		('Pt_Lep_hist', 'l p_{t}', 'Pt_Lep_hist'),
		('Pt_BHad_hist', 'b_{h} p_{t}', 'Pt_BHad_hist'),
		('Pt_BLep_hist', 'b_{l} p_{t}', 'Pt_BLep_hist'),
		('Pt_WJa_hist', 'WJa p_{t}', 'Pt_WJa_hist'),
		('Pt_WJb_hist', 'WJb p_{t}', 'Pt_WJb_hist'),
	]
	
	for var, xaxis, names in Gen_Pt:
	      hist = asrootpy(myfile.Get(var)).Clone()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax)
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	      plotter.save(names)
	
	Gen_Eta = [
		('Eta_Lep_hist', 'l #eta', 'Eta_Lep_hist'),
		('Eta_BHad_hist', 'b_{h} #eta', 'Eta_BHad_hist'),
		('Eta_BLep_hist', 'b_{l} #eta', 'Eta_BLep_hist'),
		('Eta_WJa_hist', 'WJa #eta', 'Eta_WJa_hist'),
		('Eta_WJb_hist', 'WJb #eta', 'Eta_WJb_hist'),
	]
	
	for var, xaxis, names in Gen_Eta:
	      hist = asrootpy(myfile.Get(var)).Clone()
	      mean = hist.GetMean()
	      rms = hist.GetRMS()
	
	      plotter.set_histo_style(hist, color=defcol)
	      plotter.plot(hist, xtitle=xaxis, ytitle=defyax)
	      box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
	      plotter.save(names)
	
	Gen_DR_2D = [
		('DR_LepBHad_vs_Mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (l, b_{h})', 'DR_LepBHad_vs_Mtt_hist'),
		('DR_LepBLep_vs_Mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (l, b_{l})', 'DR_LepBLep_vs_Mtt_hist'),
		('DR_LepWJa_vs_Mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (l, WJa)', 'DR_LepWJa_vs_Mtt_hist'),
		('DR_LepWJb_vs_Mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (l, WJb)', 'DR_LepWJb_vs_Mtt_hist'),
		('DR_BHadBLep_vs_Mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (b_{h}, b_{l})', 'DR_BHadBLep_vs_Mtt_hist'),
		('DR_BHadWJa_vs_Mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (b_{h}, WJa)', 'DR_BHadWJa_vs_Mtt_hist'),
		('DR_BHadWJb_vs_Mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (b_{h}, WJb)', 'DR_BHadWJb_vs_Mtt_hist'),
		('DR_BLepWJa_vs_Mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (b_{l}, WJa)', 'DR_BLepWJa_vs_Mtt_hist'),
		('DR_BLepWJb_vs_Mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (b_{l}, WJb)', 'DR_BLepWJb_vs_Mtt_hist'),
		('DR_WJaWJb_vs_Mtt_hist', 'm_{t#bar t} [GeV]', '#Delta R (WJa, WJb)', 'DR_WJaWJb_vs_Mtt_hist')
	]
	
	for var, xaxis,yaxis, names in Gen_DR_2D:
		hist = asrootpy(myfile.Get(var)).Clone()
	#	plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
		hist.Draw('colz')
		hist.xaxis.set_title(xaxis)
		hist.yaxis.set_title(yaxis)
		plotter.save(names)
	
	Gen_Pt_2D = [
		('Pt_Lep_vs_Mtt_hist', 'm_{t#bar t} [GeV]', 'l p_{t}', 'Pt_Lep_vs_Mtt_hist'),
		('Pt_BHad_vs_Mtt_hist', 'm_{t#bar t} [GeV]', 'b_{h} p_{t}', 'Pt_BHad_vs_Mtt_hist'),
		('Pt_BLep_vs_Mtt_hist', 'm_{t#bar t} [GeV]', 'b_{l} p_{t}', 'Pt_BLep_vs_Mtt_hist'),
		('Pt_WJa_vs_Mtt_hist', 'm_{t#bar t} [GeV]', 'WJa p_{t}', 'Pt_WJa_vs_Mtt_hist'),
		('Pt_WJb_vs_Mtt_hist', 'm_{t#bar t} [GeV]', 'WJb p_{t}', 'Pt_WJb_vs_Mtt_hist')
	]
	
	for var, xaxis,yaxis, names in Gen_Pt_2D:
		hist = asrootpy(myfile.Get(var)).Clone()
	#	plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
		hist.Draw('colz')
		hist.xaxis.set_title(xaxis)
		hist.yaxis.set_title(yaxis)
		plotter.save(names)
	
	Gen_Eta_2D = [
		('Eta_Lep_vs_Mtt_hist', 'm_{t#bar t} [GeV]', 'l #eta', 'Eta_Lep_vs_Mtt_hist'),
		('Eta_BHad_vs_Mtt_hist', 'm_{t#bar t} [GeV]', 'b_{h} #eta', 'Eta_BHad_vs_Mtt_hist'),
		('Eta_BLep_vs_Mtt_hist', 'm_{t#bar t} [GeV]', 'b_{l} #eta', 'Eta_BLep_vs_Mtt_hist'),
		('Eta_WJa_vs_Mtt_hist', 'm_{t#bar t} [GeV]', 'WJa #eta', 'Eta_WJa_vs_Mtt_hist'),
		('Eta_WJb_vs_Mtt_hist', 'm_{t#bar t} [GeV]', 'WJb #eta', 'Eta_WJb_vs_Mtt_hist')
	]
	
	for var, xaxis,yaxis, names in Gen_Eta_2D:
		hist = asrootpy(myfile.Get(var)).Clone()
	#	plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
		hist.Draw('colz')
		hist.xaxis.set_title(xaxis)
		hist.yaxis.set_title(yaxis)
		plotter.save(names)

	MassTTbar = [
	#	('Mu_Mttbar_hist', 'm_{t#bar t} [GeV]', 'Mass tt_bar for #mu Events', 'Muon_System_Mass'),
		('Jets_Mttbar_hist', 'm_{t#bar t} [GeV]', 'Mass tt_bar for Jets Events', 'Jets_System_Mass'),
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
			efficiencies.append(Match_Eff.GetBinContent(i+1))
		#	mass_range.append(Mu_Eff.GetXaxis().GetBinLowEdge(i+1))
		print(yaxis, efficiencies)

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
			efficiencies.append(Match_Eff.GetBinContent(i+1))
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
			efficiencies.append(Match_Eff.GetBinContent(i+1))
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
		('nJets_hist', 'nJets', 'nJets_hist')
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
