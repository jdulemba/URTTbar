'''
Jet_Match Analyzer Plotter macro
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

parser = argparse.ArgumentParser(description='Create plots using files from jet_match.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full_Analysis).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, or Combined).')
parser.add_argument('plot', help='Choose type of plots to generate (Matched_Objects, Delta_Plots, BTagging_Eff, Merged_Perms).')
args = parser.parse_args()

if args.analysis == "Test":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../%s.jet_match.test.root' % args.sample, 'read')
		normfile = views.NormalizeView(root_open('../%s.jet_match.test.root' % args.sample, 'read'))#normalized file
		plotter = BasePlotter(
			'plots/jet_match/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Sample_ttJetsM0 = root_open('../ttJetsM0.jet_match.test.root', 'read')
	#	Sample_ttJetsM700 = root_open('../ttJetsM700.jet_match.test.root', 'read')
	#	Sample_ttJetsM1000 = root_open('../ttJetsM1000.jet_match.test.root', 'read')
		plotter = BasePlotter(
			'plots/jet_match/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)

elif args.analysis == "Full_Analysis":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../results/%s/jet_match/%s.root' % (jobid, args.sample), 'read')
		normfile = views.NormalizeView(root_open('../results/%s/jet_match/%s.root' % (jobid, args.sample), 'read'))
		plotter = BasePlotter(
			'plots/jet_match/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Full_ttJetsM0 = root_open('../results/%s/jet_match/ttJets.root' % jobid, 'read')
	#	Full_ttJetsM700 = root_open('../results/%s/jet_match/ttJetsM700.root' % jobid, 'read')
	#	Full_ttJetsM1000 = root_open('../results/%s/jet_match/ttJetsM1000.root' % jobid, 'read')
		plotter = BasePlotter(
			'plots/jet_match/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
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

##############################################################################################
if args.plot == "Matched_Objects":

	if args.sample == "ttJetsM0":
		Matched_WJets = [
			('Matched_perm_WJa_ttbarM700_pt', 'Matched WJa p_{t}', '700 <= m_{t#bar t} < 1000 [GeV]', 'Matched_WJa_ttbarM700_pt'),
			('Matched_perm_WJa_ttbarM1000_pt', 'Matched WJa p_{t}', 'm_{t#bar t} >= 1000 [GeV]', 'Matched_WJa_ttbarM1000_pt'),
			('Matched_perm_WJb_ttbarM700_pt', 'Matched WJb p_{t}', '700 <= m_{t#bar t} < 1000 [GeV]', 'Matched_WJb_ttbarM700_pt'),
			('Matched_perm_WJb_ttbarM1000_pt', 'Matched WJb p_{t}', 'm_{t#bar t} >= 1000 [GeV]', 'Matched_WJb_ttbarM1000_pt')
		]

		for var, xaxis, titles, names in Matched_WJets:
			hist = asrootpy(myfile.Get(var)).Clone()
			mean = hist.GetMean()
			rms = hist.GetRMS()
			plotter.set_histo_style(hist, color=defcol)
			plotter.plot(hist, xtitle=xaxis, ytitle=defyax)
			box = plotter.make_text_box(titles+'\nMean = %f\nRMS = %f' % (mean,rms), position='NE')
			plotter.save(names)

		Matched_WJets_Frac = [
			('Matched_perm_WJa_ttbarM700_frac_p', 'Matched WJa Fractional p', '700 <= m_{t#bar t} < 1000 [GeV]', 'Matched_WJa_ttbarM700_frac_p'),
			('Matched_perm_WJa_ttbarM1000_frac_p', 'Matched WJa Fractional p', 'm_{t#bar t} >= 1000 [GeV]', 'Matched_WJa_ttbarM1000_frac_p'),
			('Matched_perm_WJb_ttbarM700_frac_p', 'Matched WJb Fractional p', '700 <= m_{t#bar t} < 1000 [GeV]', 'Matched_WJb_ttbarM700_frac_p'),
			('Matched_perm_WJb_ttbarM1000_frac_p', 'Matched WJb Fractional p', 'm_{t#bar t} >= 1000 [GeV]', 'Matched_WJb_ttbarM1000_frac_p')
		]

		for var, xaxis, titles, names in Matched_WJets_Frac:
			hist = asrootpy(myfile.Get(var)).Clone()
			mean = hist.GetMean()
			rms = hist.GetRMS()
#			hist.xaxis.range_user = 0, 1.0
			plotter.set_histo_style(hist, color=defcol)
			plotter.plot(hist, xtitle=xaxis, ytitle=defyax)
			box = plotter.make_text_box(titles+'\nMean = %f\nRMS = %f' % (mean,rms), position='NE')
			plotter.save(names)

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
#		('All_BHad_events_vs_mttbar', 'Merged_BHadBLep_vs_mttbar', 'b_{h} & b_{l}', 'Merged_BHadBLep_Frac'),
		('All_BHad_events_vs_mttbar', 'Merged_BHadWJa_vs_mttbar', 'b_{h} & WJa', 'Merged_BHadWJa_Frac'),
		('All_BHad_events_vs_mttbar', 'Merged_BHadWJb_vs_mttbar', 'b_{h} & WJb', 'Merged_BHadWJb_Frac'),
#		('All_BLep_events_vs_mttbar', 'Merged_BHadBLep_vs_mttbar', 'b_{l} & b_{h}', 'Merged_BLepBHad_Frac'),
#		('All_BLep_events_vs_mttbar', 'Merged_BLepWJa_vs_mttbar', 'b_{l} & WJa', 'Merged_BLepWJa_Frac'),
#		('All_BLep_events_vs_mttbar', 'Merged_BLepWJb_vs_mttbar', 'b_{l} & WJb', 'Merged_BLepWJb_Frac'),
		('All_WJa_events_vs_mttbar', 'Merged_BHadWJa_vs_mttbar', 'WJa & b_{h}', 'Merged_WJaBHad_Frac'),
#		('All_WJa_events_vs_mttbar', 'Merged_BLepWJa_vs_mttbar', 'WJa & b_{l}', 'Merged_WJaBLep_Frac'),
		('All_WJa_events_vs_mttbar', 'Merged_WJaWJb_vs_mttbar', 'WJa & WJb', 'Merged_WJaWJb_Frac'),
		('All_WJb_events_vs_mttbar', 'Merged_WJaWJb_vs_mttbar', 'WJb & WJa', 'Merged_WJbWJa_Frac'),
		('All_WJb_events_vs_mttbar', 'Merged_BHadWJb_vs_mttbar', 'WJb & b_{h}', 'Merged_WJbBHad_Frac'),
#		('All_WJb_events_vs_mttbar', 'Merged_BLepWJb_vs_mttbar', 'WJb & b_{l}', 'Merged_WJbBLep_Frac'),
		('All_BHad_events_vs_mttbar', 'Merged_BHadWJaWJb_vs_mttbar', 'b_{h} & WJa & WJb', 'Merged_BHadWJaWJb_Frac'),
		('All_WJa_events_vs_mttbar', 'Merged_BHadWJaWJb_vs_mttbar', 'WJa & WJb & b_{h}', 'Merged_WJaWJbBHad_Frac'),
		('All_WJb_events_vs_mttbar', 'Merged_BHadWJaWJb_vs_mttbar', 'WJb & b_{h} & WJa', 'Merged_WJbBHadWJa_Frac')
	]

	jet_name = ['_3J', '_4J', '_5PJ']
	njets = ['nJets = 3', 'nJets = 4', 'nJets > 4']

	for i in range(0,3):
		for all_events, merged, objects, name in Merged_AllEvents_Ratio:
	
			to_draw = []
			
			for folder, legend, colors in DRvals:
				All = asrootpy(myfile.Get(folder+'/'+all_events+jet_name[i])).Clone()
				Merged = asrootpy(myfile.Get(folder+'/'+merged+jet_name[i])).Clone()
				Merged_Frac = Merged/All
				plotter.set_histo_style(Merged_Frac, title=legend, color=colors, markerstyle=20)
				if Merged_Frac.Integral() != 0:
					plotter.plot(Merged_Frac)
					Merged_Frac.Draw("P")
				#	Merged_Frac.xaxis.set_title(xaxis)
				#	Merged_Frac.yaxis.set_title(yaxis)
					Merged_Frac.yaxis.range_user = 0, 1
				to_draw.append(Merged_Frac)
	
			plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle='m_{t#bart}', ytitle='Fraction of Merged Events/ Total Events')
			box = plotter.make_text_box(njets[i]+'\n'+objects, position='NW')
			plotter.save(name+jet_name[i])

	Merged_AllEvents_Ratio = [
#		('All_BHad_events_vs_mttbar', 'Merged_BHadBLep_vs_mttbar', 'b_{h} & b_{l}', 'Merged_BHadBLep_Frac'),
		('All_BHad_events_vs_mttbar', 'Merged_BHadWJa_vs_mttbar', 'b_{h} & WJa', 'Merged_BHadWJa_Frac'),
		('All_BHad_events_vs_mttbar', 'Merged_BHadWJb_vs_mttbar', 'b_{h} & WJb', 'Merged_BHadWJb_Frac'),
#		('All_BLep_events_vs_mttbar', 'Merged_BHadBLep_vs_mttbar', 'b_{l} & b_{h}', 'Merged_BLepBHad_Frac'),
#		('All_BLep_events_vs_mttbar', 'Merged_BLepWJa_vs_mttbar', 'b_{l} & WJa', 'Merged_BLepWJa_Frac'),
#		('All_BLep_events_vs_mttbar', 'Merged_BLepWJb_vs_mttbar', 'b_{l} & WJb', 'Merged_BLepWJb_Frac'),
		('All_WJa_events_vs_mttbar', 'Merged_BHadWJa_vs_mttbar', 'WJa & b_{h}', 'Merged_WJaBHad_Frac'),
#		('All_WJa_events_vs_mttbar', 'Merged_BLepWJa_vs_mttbar', 'WJa & b_{l}', 'Merged_WJaBLep_Frac'),
		('All_WJa_events_vs_mttbar', 'Merged_WJaWJb_vs_mttbar', 'WJa & WJb', 'Merged_WJaWJb_Frac'),
		('All_WJb_events_vs_mttbar', 'Merged_WJaWJb_vs_mttbar', 'WJb & WJa', 'Merged_WJbWJa_Frac'),
		('All_WJb_events_vs_mttbar', 'Merged_BHadWJb_vs_mttbar', 'WJb & b_{h}', 'Merged_WJbBHad_Frac'),
#		('All_WJb_events_vs_mttbar', 'Merged_BLepWJb_vs_mttbar', 'WJb & b_{l}', 'Merged_WJbBLep_Frac'),
		('All_BHad_events_vs_mttbar', 'Merged_BHadWJaWJb_vs_mttbar', 'b_{h} & WJa & WJb', 'Merged_BHadWJaWJb_Frac'),
		('All_WJa_events_vs_mttbar', 'Merged_BHadWJaWJb_vs_mttbar', 'WJa & WJb & b_{h}', 'Merged_WJaWJbBHad_Frac'),
		('All_WJb_events_vs_mttbar', 'Merged_BHadWJaWJb_vs_mttbar', 'WJb & b_{h} & WJa', 'Merged_WJbBHadWJa_Frac')
	]

	jet_name = ['_3J', '_4J', '_5PJ']
	njets = ['nJets = 3', 'nJets = 4', 'nJets > 4']

	for all_events, merged, objects, name in Merged_AllEvents_Ratio:
		to_draw = []

		All_3J = asrootpy(myfile.Get(folder+'/'+all_events+jet_name[0])).Clone()
		Merged_3J = asrootpy(myfile.Get(folder+'/'+merged+jet_name[0])).Clone()
		Merged_Frac_3J = Merged_3J/All_3J
		plotter.set_histo_style(Merged_Frac_3J, title=njets[0], color='black', markerstyle=20)
		if Merged_Frac_3J.Integral() != 0:
			plotter.plot(Merged_Frac_3J)
			Merged_Frac_3J.Draw("P")
		to_draw.append(Merged_Frac_3J)

		All_4J = asrootpy(myfile.Get(folder+'/'+all_events+jet_name[1])).Clone()
		Merged_4J = asrootpy(myfile.Get(folder+'/'+merged+jet_name[1])).Clone()
		Merged_Frac_4J = Merged_4J/All_4J
		plotter.set_histo_style(Merged_Frac_4J, title=njets[1], color='red', markerstyle=20)
		if Merged_Frac_4J.Integral() != 0:
			plotter.plot(Merged_Frac_4J)
			Merged_Frac_4J.Draw("P")
		to_draw.append(Merged_Frac_4J)

		All_5PJ = asrootpy(myfile.Get(folder+'/'+all_events+jet_name[2])).Clone()
		Merged_5PJ = asrootpy(myfile.Get(folder+'/'+merged+jet_name[2])).Clone()
		Merged_Frac_5PJ = Merged_5PJ/All_5PJ
		plotter.set_histo_style(Merged_Frac_5PJ, title=njets[2], color='blue', markerstyle=20)
		if Merged_Frac_5PJ.Integral() != 0:
			plotter.plot(Merged_Frac_5PJ)
			Merged_Frac_5PJ.Draw("P")
		to_draw.append(Merged_Frac_5PJ)

		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle='m_{t#bart}', ytitle='Fraction of Merged Events/ Total Events')
		box = plotter.make_text_box(objects+'\n #Delta R < 0.4', position='NW')
		plotter.save(name+'_njets')


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


	Merged_BWJet_DRBotherW = [ #DR between merged BWJet and other WJet
		('Merged_BHadWJa_DRBHadWJb', 'Merged b_{h} & WJa', '#Delta R (merged, WJb)'),
		('Merged_BHadWJb_DRBHadWJa', 'Merged b_{h} & WJb', '#Delta R (merged, WJa)')
	]

	for merged, legend_title, xaxis in Merged_BWJet_DRBotherW:
		to_draw = []
		for folder, legend, colors in DRvals:
			Merged = asrootpy(myfile.Get(folder+'/'+merged)).Clone()
			plotter.set_histo_style(Merged, title=legend, color=colors)
			to_draw.append(Merged)
			
#			WJb = asrootpy(myfile.Get(folder+'/'+wjb)).Clone()
#			plotter.set_histo_style(WJb, title=b_type+' & WJb', color='blue')
#			to_draw.append(WJb)
			
		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle=xaxis, ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box(legend_title, position='NW')
		plotter.save(merged)

	Merged_BHadWJaWJb = [ # 3 jets from hadronic side are merged
		('Merged_BHadWJaWJb_mass', 'm_{merged} [GeV]'),
		('Merged_BHadWJaWJb_pt', 'Merged jet p_{T}'),
		('Merged_BHadWJaWJb_eta', 'Merged jet #eta')
	]

	for kin_var, xaxis in Merged_BHadWJaWJb:
		to_draw = []
		for folder, legend, colors in DRvals:
			hist = asrootpy(myfile.Get(folder+'/'+kin_var)).Clone()
			plotter.set_histo_style(hist, title=legend, color=colors)
			to_draw.append(hist)
			
		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle=xaxis, ytitle=defyax, drawstyle='hist')
		box = plotter.make_text_box('Merged b_{h} & WJets', position='NW')
		plotter.save(kin_var)


#############################################################################################
if args.plot == "Delta_Plots":
	# Delta (Gen - Reco) Plots
	
	THad_Delta = [
	      ('Merged_THad_Delta_costh', 'Unmerged_THad_Delta_costh', '#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', -0.5, 0.5, 'THad_Delta_costh'),
	      ('Merged_THad_Delta_mass', 'Unmerged_THad_Delta_mass', '#Delta m_{t_{h}} (Gen - Reco) [GeV]', -200, 200, 'THad_Delta_mass')
	]
	
	for merged, unmerged, xaxis, xmin, xmax, name in THad_Delta:
		to_draw = []
		
		Merged = asrootpy(normfile.Get('Delta_Plots/'+merged)).Clone()
		plotter.set_histo_style(Merged, title='Merged', color='black')
		to_draw.append(Merged)
		
		Unmerged = asrootpy(normfile.Get('Delta_Plots/'+unmerged)).Clone()
		plotter.set_histo_style(Unmerged, title='Unmerged', color='red')
		to_draw.append(Unmerged)
		
		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), x_range=(xmin,xmax), xtitle=xaxis, drawstyle='hist')
#		box = plotter.make_text_box(xaxis, position='NW')
		plotter.save(name)


#############################################################

print('~/nobackup/CMSSW_8_0_7_patch2/src/Analyses/URTTbar/htt_scripts/plots/jet_match/%s/%s/%s/%s/' % (jobid, args.analysis, args.plot, args.sample))
