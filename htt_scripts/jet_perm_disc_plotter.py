'''
Jet_Perm_Disc Analyzer Plotter macro
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

parser = argparse.ArgumentParser(description='Create plots using files from jet_perm_disc.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000).')
parser.add_argument('plot', help='Choose type of plots to generate (Had_Comp, Delta_Plots, Perm_Disc, NS_Chi, NS_Disc, Combined_Disc, Combine_Discs).')
args = parser.parse_args()

if args.analysis == "Test":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../%s.jet_perm_disc.test.root' % args.sample, 'read')
		normfile = views.NormalizeView(root_open('../%s.jet_perm_disc.test.root' % args.sample, 'read'))#normalized file
		plotter = BasePlotter(
			'plots/jet_perm_disc/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Sample_ttJetsM0 = root_open('../ttJetsM0.jet_perm_disc.test.root', 'read')
	#	Sample_ttJetsM700 = root_open('../ttJetsM700.jet_perm_disc.test.root', 'read')
	#	Sample_ttJetsM1000 = root_open('../ttJetsM1000.jet_perm_disc.test.root', 'read')
		plotter = BasePlotter(
			'plots/jet_perm_disc/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)

elif args.analysis == "Full":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../results/%s/jet_perm_disc/%s.root' % (jobid, args.sample), 'read')
		normfile = views.NormalizeView(root_open('../results/%s/jet_perm_disc/%s.root' % (jobid, args.sample), 'read'))
		plotter = BasePlotter(
			'plots/jet_perm_disc/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)
	elif args.sample == "Combined":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
	#	Full_ttJetsM0 = root_open('../results/%s/jet_perm_disc/ttJets.root' % jobid, 'read')
	#	Full_ttJetsM700 = root_open('../results/%s/jet_perm_disc/ttJetsM700.root' % jobid, 'read')
	#	Full_ttJetsM1000 = root_open('../results/%s/jet_perm_disc/ttJetsM1000.root' % jobid, 'read')
		plotter = BasePlotter(
			'plots/jet_perm_disc/%s/%s/%s/%s' % (jobid, args.analysis, args.plot, args.sample),
			defaults = {'save' : {'png' : True, 'pdf' : False}}
		)


def stack_plots(lists):
    lists.sort()
    total = 0
#    ratio_hists = []
    Stack_hists = []

    for i in lists:
        total += i

    for i in lists:
#        ratio_hists.append(i/total)
        i.SetFillStyle(1001)
        Stack_hists.append(i)

    if len(Stack_hists) == 4:
        stack = plotter.create_stack(Stack_hists[0], Stack_hists[1], Stack_hists[2], Stack_hists[3])
        norm_stack = plotter.create_stack(Stack_hists[0]/total, Stack_hists[1]/total, Stack_hists[2]/total, Stack_hists[3]/total)
    elif len(Stack_hists) == 3:
        stack = plotter.create_stack(Stack_hists[0], Stack_hists[1], Stack_hists[2])
        norm_stack = plotter.create_stack(Stack_hists[0]/total, Stack_hists[1]/total, Stack_hists[2]/total)
        
    return stack, norm_stack
    

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
	

##############################################################################################

if args.plot == "Mass_Disc":


        #### 3J Perm Disc  values
    Likelihood_All_3J = [
        ('3J_Permdisc_All', 'All Mass Disc Values from Event', 'black', '-Log(likelihood)', '3J_Massdisc_Event_All'),
        ('3J_Permdisc_Lowest', 'Lowest Mass Disc Value from Event', 'black', '-Log(likelihood)', '3J_Massdisc_Event_Lowest'),
        ('3J_Permdisc_Best_Perm', 'Perm Mass Value from Event Best Perm', 'black', '-Log(likelihood)', '3J_Massdisc_Event_Best_Perm')
    ]
   
    for var, legends, colors, xaxis, names in Likelihood_All_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
        plotter.plot(hist, drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(names)

#        #from pdb import set_trace; set_trace()
#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, drawstyle='hist')
#        plotter.save(names+'_Norm')

    ### Plots for lowest likelihood values from event

            ### 3 jets

    Matched_Perm_Likelihood_3J = [
        ('3J_Permdisc_Best_Perm_RIGHT', 'RIGHT', 'black', '-Log(likelihood)'),
        ('3J_Permdisc_Best_Perm_MERGE_SWAP', 'MERGE_SWAP', 'red', '-Log(likelihood)'),
        ('3J_Permdisc_Best_Perm_MERGE', 'MERGE', 'green', '-Log(likelihood)'),
        ('3J_Permdisc_Best_Perm_WRONG', 'WRONG', 'blue', '-Log(likelihood)')
    ]
   
    to_draw = [] 
    for var, legends, colors, xaxis in Matched_Perm_Likelihood_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(var+'_Categories')
        to_draw.append(hist)

#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
##        plotter.save(var+'_Matched_Perm_Norm')
#        norm_to_draw.append(norm_hist)

    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Massdisc_Best_Perm_Categories')

#    plotter.overlay(norm_to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#    plotter.save('3J_Permdisc_Best_Perm_Matched_Perm_Norm')
   
        
    stack, norm_stack = stack_plots(to_draw)

    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Massdisc_Best_Perm_Categories_Stack')

    plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Massdisc_Best_Perm_Categories_Stack_Norm')


##############################################################################################

if args.plot == "NS_Chi":


        #### 3J NS chi^2  values
    Likelihood_All_3J = [
        ('3J_NSchi_All', 'All Values from Event', 'black', '#chi^{2}', '3J_NSchi_Event_All'),
        ('3J_NSchi_Lowest', 'Lowest Value from Event', 'black', '#chi^{2}', '3J_NSchi_Event_Lowest'),
        ('3J_NSchi_Best_Perm', 'Value from Event Best Perm', 'black', '#chi^{2}', '3J_NSchi_Event_Best_Perm')
    ]
   
    for var, legends, colors, xaxis, names in Likelihood_All_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
        plotter.plot(hist, drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(names)

#        #from pdb import set_trace; set_trace()
#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, drawstyle='hist')
#        plotter.save(names+'_Norm')

    ### Plots for lowest likelihood values from event

            ### 3 jets

    Matched_Perm_Likelihood_3J = [
        ('3J_NSchi_Best_Perm_RIGHT', 'RIGHT', 'black', '#chi^{2}'),
        ('3J_NSchi_Best_Perm_MERGE_SWAP', 'MERGE_SWAP', 'red', '#chi^{2}'),
        ('3J_NSchi_Best_Perm_MERGE', 'MERGE', 'green', '#chi^{2}'),
        ('3J_NSchi_Best_Perm_WRONG', 'WRONG', 'blue', '#chi^{2}')
    ]
   
    to_draw = [] 
    for var, legends, colors, xaxis in Matched_Perm_Likelihood_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(var+'_Categories')
        to_draw.append(hist)

#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#        plotter.save(var+'_Matched_Perm_Norm')
#        norm_to_draw.append(norm_hist)

    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_NSchi_Best_Perm_Categories')

#    plotter.overlay(norm_to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#    plotter.save('3J_NSchi_Best_Perm_Matched_Perm_Norm')

        
    stack, norm_stack = stack_plots(to_draw)

    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_NSchi_Best_Perm_Categories_Stack')

    plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_NSchi_Best_Perm_Categories_Stack_Norm')


##############################################################################################

if args.plot == "NS_Disc":


        #### 3J NS Disc  values
    Likelihood_All_3J = [
        ('3J_NSdisc_All', 'All NS Disc Values from Event', 'black', '-Log(likelihood)', '3J_NSdisc_Event_All'),
        ('3J_NSdisc_Lowest', 'Lowest NS Disc Value from Event', 'black', '-Log(likelihood)', '3J_NSdisc_Event_Lowest'),
        ('3J_NSdisc_Best_Perm', 'NS Disc Value from Best Perm', 'black', '-Log(likelihood)', '3J_NSdisc_Event_Best_Perm')
    ]
   
    for var, legends, colors, xaxis, names in Likelihood_All_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
        plotter.plot(hist, drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(names)

#        #from pdb import set_trace; set_trace()
#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, drawstyle='hist')
#        plotter.save(names+'_Norm')

    ### Plots for lowest likelihood values from event

            ### 3 jets

    Matched_Perm_Likelihood_3J = [
        ('3J_NSdisc_Best_Perm_RIGHT', 'RIGHT', 'black', '-Log(likelihood)'),
        ('3J_NSdisc_Best_Perm_MERGE_SWAP', 'MERGE_SWAP', 'red', '-Log(likelihood)'),
        ('3J_NSdisc_Best_Perm_MERGE', 'MERGE', 'green', '-Log(likelihood)'),
        ('3J_NSdisc_Best_Perm_WRONG', 'WRONG', 'blue', '-Log(likelihood)')
    ]
   
    to_draw = [] 
    for var, legends, colors, xaxis in Matched_Perm_Likelihood_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(var+'_Categories')
        to_draw.append(hist)

#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#        plotter.save(var+'_Matched_Perm_Norm')
#        norm_to_draw.append(norm_hist)

    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_NSdisc_Best_Perm_Categories')

#    plotter.overlay(norm_to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#    plotter.save('3J_Permdisc_Best_Perm_Matched_Perm_Norm')
        
    stack, norm_stack = stack_plots(to_draw)

    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_NSdisc_Best_Perm_Categories_Stack')

    plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_NSdisc_Best_Perm_Categories_Stack_Norm')


##############################################################################################

if args.plot == "Combined_Disc":


        #### 3J Combined Disc  values
    Likelihood_All_3J = [
        ('3J_Totaldisc_All', 'All Combined Disc Values from Event', 'black', '-Log(likelihood)', '3J_Combined_Disc_Event_All'),
        ('3J_Totaldisc_Lowest', 'Lowest Combined Disc Value from Event', 'black', '-Log(likelihood)', '3J_Combined_Disc_Event_Lowest'),
        ('3J_Totaldisc_Best_Perm', 'Combined Disc Value from Event Best Perm', 'black', '-Log(likelihood)', '3J_Combined_Disc_Event_Best_Perm')
    ]
   
    for var, legends, colors, xaxis, names in Likelihood_All_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
        plotter.plot(hist, drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(names)

#        #from pdb import set_trace; set_trace()
#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, drawstyle='hist')
#        plotter.save(names+'_Norm')

    ### Plots for lowest likelihood values from event

            ### 3 jets

    Matched_Perm_Likelihood_3J = [
        ('3J_Totaldisc_Best_Perm_RIGHT', 'RIGHT', 'black', '-Log(likelihood)', '3J_Combined_Disc_Best_Perm_RIGHT'),
        ('3J_Totaldisc_Best_Perm_MERGE_SWAP', 'MERGE_SWAP', 'red', '-Log(likelihood)', '3J_Combined_Disc_Best_Perm_MERGE_SWAP'),
        ('3J_Totaldisc_Best_Perm_MERGE', 'MERGE', 'green', '-Log(likelihood)', '3J_Combined_Disc_Best_Perm_MERGE'),
        ('3J_Totaldisc_Best_Perm_WRONG', 'WRONG', 'blue', '-Log(likelihood)', '3J_Combined_Disc_Best_Perm_WRONG')
    ]
   
    to_draw = [] 
    for var, legends, colors, xaxis, names in Matched_Perm_Likelihood_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(names+'_Categories')
        to_draw.append(hist)

#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#        plotter.save(var+'_Matched_Perm_Norm')
#        norm_to_draw.append(norm_hist)

    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Combined_Disc_Best_Perm_Categories')

#    plotter.overlay(norm_to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#    plotter.save('3J_Permdisc_Best_Perm_Matched_Perm_Norm')

        
    stack, norm_stack = stack_plots(to_draw)

    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Combined_Disc_Best_Perm_Categories_Stack')

    plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Combined_Disc_Best_Perm_Categories_Stack_Norm')

    

##############################################################################################

if args.plot == "Combine_Discs":

    Likelihood_All_3J = [
        ('3J_Totaldisc_All', 'Combined Disc', 'black'),
        ('3J_Permdisc_All', 'Mass', 'red'),
        ('3J_NSdisc_All', 'Nu Solver', 'blue')
    ]

    to_draw = []

    for var, legends, colors in Likelihood_All_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle='-Log(likelihood)')
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        to_draw.append(hist)


    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Alldiscs_Event_All')


    Likelihood_Lowest_3J = [
        ('3J_Totaldisc_Lowest', 'Combined Disc', 'black'),
        ('3J_Permdisc_Lowest', 'Mass', 'red'),
        ('3J_NSdisc_Lowest', 'Nu Solver', 'blue')
    ]

    to_draw = []

    for var, legends, colors in Likelihood_All_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle='-Log(likelihood)')
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        to_draw.append(hist)


    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Alldiscs_Event_Lowest')


    Likelihood_Best_Perm_3J = [
        ('3J_Totaldisc_Best_Perm', 'Combined Disc', 'black'),
        ('3J_Permdisc_Best_Perm', 'Mass', 'red'),
        ('3J_NSdisc_Best_Perm', 'Nu Solver', 'blue')
    ]

    to_draw = []

    for var, legends, colors in Likelihood_Best_Perm_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        plotter.set_histo_style(hist, title=legends, color=colors, xtitle='-Log(likelihood)')
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        to_draw.append(hist)


    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Alldiscs_Event_Best_Perm')



##############################################
print('cp ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/%s/*.png .' % (jobid, args.analysis, args.plot, args.sample))
