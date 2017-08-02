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
parser.add_argument('plot', help='Choose type of plots to generate (Had_Comp, Delta_Plots, Discs, NS_Chi, Combine_Discs, TP_Tagger, MP_Tagger, Best_Perm, Best_Perm_Combined).')
args = parser.parse_args()

if args.analysis == "Test":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../%s.jet_perm_disc.test.root' % args.sample, 'read')
		normfile = views.NormalizeView(root_open('../%s.jet_perm_disc.test.root' % args.sample, 'read'))#normalized file
		plotter = BasePlotter(
			'plots/jet_perm_disc/%s/%s/%s' % (jobid, args.analysis, args.sample),
			defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}}
#			defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
		)

elif args.analysis == "Full":
	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
		print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
		myfile = root_open('../results/%s/jet_perm_disc/%s.root' % (jobid, args.sample), 'read')
		normfile = views.NormalizeView(root_open('../results/%s/jet_perm_disc/%s.root' % (jobid, args.sample), 'read'))
		plotter = BasePlotter(
			'plots/jet_perm_disc/%s/%s/%s' % (jobid, args.analysis, args.sample),
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
    elif len(Stack_hists) == 2:
        stack = plotter.create_stack(Stack_hists[0], Stack_hists[1])
        norm_stack = plotter.create_stack(Stack_hists[0]/total, Stack_hists[1]/total)
        
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

if args.plot == "Discs":

    disc_dict = {'Massdisc' : ['Mass Disc', 'black', 'Mass_Disc'], 'NSdisc' : ['NS Disc', 'black', 'NS_Disc'], 'Totaldisc' : ['Total Disc', 'black', 'Total_Disc']}
    for disc in disc_dict:
        disc_lab = disc_dict[disc][0]
        disc_col = disc_dict[disc][1]
        disc_dir = disc_dict[disc][2]

        a_dir = 'All_Perms/%s' % disc_dir
        plotter.set_subdir(a_dir)

        #### 3J Disc  values
        Likelihood_All_3J = [
            ('3J_Event_All_', 'All %s Values from Event' % disc_lab, 'black', '-Log(likelihood)', '3J_Event_All_'),
            ('3J_Event_Lowest_', 'Lowest %s Value from Event' % disc_lab, 'black', '-Log(likelihood)', '3J_Event_Lowest_'),
            ('3J_Event_Best_Perm_', 'Perm %s Value from Event Best Perm' % disc_lab, 'black', '-Log(likelihood)', '3J_Event_Best_Perm_')
        ]
       
        for var, legends, colors, xaxis, names in Likelihood_All_3J:
            hist = asrootpy(myfile.Get('Likelihood_Plots/'+var+disc)).Clone()
            if hist.Integral() == 0:
                print var
                continue
    
            mean = hist.GetMean()
            rms = hist.GetRMS()
    
            plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
            plotter.plot(hist, drawstyle='hist')
            box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
            box.Draw()
            plotter.save(names+disc)
    
#        #from pdb import set_trace; set_trace()
#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, drawstyle='hist')
#        plotter.save(names+'_Norm')

    ### Plots for lowest likelihood values from event

        m_dir = 'Matched_Perms/%s' % disc_dir
        plotter.set_subdir(m_dir)

            ### 3 jets

        Matched_Perm_Likelihood_3J = [
            ('3J_Best_Perm_RIGHT_', 'RIGHT', 'black', '-Log(likelihood)'),
            ('3J_Best_Perm_MERGE_SWAP_', 'MERGE_SWAP', 'red', '-Log(likelihood)'),
            ('3J_Best_Perm_MERGE_', 'MERGE', 'green', '-Log(likelihood)'),
            ('3J_Best_Perm_WRONG_', 'WRONG', 'blue', '-Log(likelihood)')
        ]
       
        to_draw = [] 
        for var, legends, colors, xaxis in Matched_Perm_Likelihood_3J:
            hist = asrootpy(myfile.Get('Likelihood_Plots/'+var+disc)).Clone()
            if hist.Integral() == 0:
                print var
                continue
    
            mean = hist.GetMean()
            rms = hist.GetRMS()
    
            plotter.set_histo_style(hist, title=legends, color=disc_col, xtitle=xaxis)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
            box.Draw()
            plotter.save(var+disc+'_Categories')

            plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(hist)
    
    #        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
    #        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
    #        plotter.plot(norm_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    ##        plotter.save(var+'_Matched_Perm_Norm')
    #        norm_to_draw.append(norm_hist)
    
        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        plotter.save('3J_Best_Perm_Categories_'+disc)
    
    #    plotter.overlay(norm_to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    #    plotter.save('3J_Permdisc_Best_Perm_Matched_Perm_Norm')
       
            
        stack, norm_stack = stack_plots(to_draw)
    
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        plotter.save('3J_Best_Perm_Categories_Stack_'+disc)
    
        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)



    print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/All_Perms .\n' % (jobid, args.analysis, args.sample)

    print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/Matched_Perms .\n' % (jobid, args.analysis, args.sample)



##############################################################################################


if args.plot == "NS_Chi":

    a_dir = 'All_Perms/%s' % args.plot
    plotter.set_subdir(a_dir)

        #### 3J NS chi^2  values
    Likelihood_All_3J = [
        ('3J_Event_All_NSchi', 'All Values from Event', 'black', '#chi^{2}'),
        ('3J_Event_Lowest_NSchi', 'Lowest Value from Event', 'black', '#chi^{2}'),
        ('3J_Event_Best_Perm_NSchi', 'Value from Event Best Perm', 'black', '#chi^{2}')
    ]
   
    for var, legends, colors, xaxis in Likelihood_All_3J:
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
        plotter.save(var)

#        #from pdb import set_trace; set_trace()
#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, drawstyle='hist')
#        plotter.save(names+'_Norm')

    ### Plots for lowest likelihood values from event

            ### 3 jets

    m_dir = 'Matched_Perms/%s' % args.plot
    plotter.set_subdir(m_dir)

    Matched_Perm_Likelihood_3J = [
        ('3J_Best_Perm_RIGHT_NSchi', 'RIGHT', 'black', '#chi^{2}'),
        ('3J_Best_Perm_MERGE_SWAP_NSchi', 'MERGE_SWAP', 'red', '#chi^{2}'),
        ('3J_Best_Perm_MERGE_NSchi', 'MERGE', 'green', '#chi^{2}'),
        ('3J_Best_Perm_WRONG_NSchi', 'WRONG', 'blue', '#chi^{2}')
    ]
   
    to_draw = [] 
    for var, legends, colors, xaxis in Matched_Perm_Likelihood_3J:
        hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, title=legends, color='black', xtitle=xaxis)
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(var+'_Categories')
        
        plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis)
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        to_draw.append(hist)
    
    
#        norm_hist = asrootpy(normfile.Get('Likelihood_Plots/'+var)).Clone()
#        plotter.set_histo_style(norm_hist, title=legends, color=colors, xtitle=xaxis)
#        plotter.plot(norm_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#        plotter.save(var+'_Matched_Perm_Norm')
#        norm_to_draw.append(norm_hist)

    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Best_Perm_Categories_NSchi')

#    plotter.overlay(norm_to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#    plotter.save('3J_NSchi_Best_Perm_Matched_Perm_Norm')

        
    stack, norm_stack = stack_plots(to_draw)

    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Best_Perm_Categories_Stack_NSchi')

    plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    plotter.save('3J_Best_Perm_Categories_Stack_Norm_NSchi')



    print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/All_Perms .\n' % (jobid, args.analysis, args.sample)

    print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/Matched_Perms .\n' % (jobid, args.analysis, args.sample)



##############################################################################################

if args.plot == "Combine_Discs":

    a_dir = 'All_Perms/'
    plotter.set_subdir(a_dir)

    Likelihood_All_3J = [
        ('3J_Event_All_Totaldisc', 'Total Disc', 'black'),
        ('3J_Event_All_Massdisc', 'Mass', 'red'),
        ('3J_Event_All_NSdisc', 'Nu Solver', 'blue')
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
    plotter.save('3J_Event_All_Disc_Comp')


    Likelihood_Lowest_3J = [
        ('3J_Event_Lowest_Totaldisc', 'Total Disc', 'black'),
        ('3J_Event_Lowest_Massdisc', 'Mass', 'red'),
        ('3J_Event_Lowest_NSdisc', 'Nu Solver', 'blue')
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
    plotter.save('3J_Event_Lowest_Disc_Comp')


    Likelihood_Best_Perm_3J = [
        ('3J_Event_Best_Perm_Totaldisc', 'Total Disc', 'black'),
        ('3J_Event_Best_Perm_Massdisc', 'Mass', 'red'),
        ('3J_Event_Best_Perm_NSdisc', 'Nu Solver', 'blue')
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
    plotter.save('3J_Event_Best_Perm_Disc_Comp')

    print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/All_Perms .\n' % (jobid, args.analysis, args.sample)



##############################################################################################

if args.plot == "TP_Tagger":

    s_dir = 'Tagger_Vars/Test_Perm/Mass'
    plotter.set_subdir(s_dir)

    Test_Perm_Mass_3J = [
        ('3J_bj1wj1_Mass', 'M_{b_{1}+w_{1}} [GeV]', 'black'),
        ('3J_bj2wj1_Mass', 'M_{b_{2}+w_{1}} [GeV]', 'black')
    ]

    for var, xaxis, colors in Test_Perm_Mass_3J:
        hist = asrootpy(myfile.Get('Test_Perm_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, color=colors, xtitle=xaxis)
        plotter.plot(hist, drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(var+'_Test_Perm')

    tp_dict = {'Mass' : [' Mass [GeV]', 'black'], 'CSV' : [' CSV', 'black'], 'CvsL' : [' CvsL', 'black'], 'Mult' : [' Multiplicity', 'black']}

    Test_Perm_3J = [
        ('3J_bj1_', 'b_{1}', 'black'),
        ('3J_bj2_', 'b_{2}', 'blue'),
        ('3J_wj1_', 'w_{1}', 'red')
    ]

    for key in tp_dict:
        dict_name = tp_dict[key][0]
        dict_color = tp_dict[key][1]

        s_dir = 'Tagger_Vars/Test_Perm/%s' % key
        plotter.set_subdir(s_dir)

        to_draw = []
    
        for var, objects, colors in Test_Perm_3J:
            hist = asrootpy(myfile.Get('Test_Perm_Plots/'+var+key)).Clone()
            if hist.Integral() == 0:
                print var
                continue
 
            mean = hist.GetMean()
            rms = hist.GetRMS()
    
            plotter.set_histo_style(hist, color=dict_color, xtitle=objects+dict_name)
            plotter.plot(hist, drawstyle='hist')
            box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
            box.Draw()
            plotter.save(var+key+'_Test_Perm')

            plotter.set_histo_style(hist, title=objects, color=colors, xtitle=dict_name)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(hist)


        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        plotter.save('3J_'+key+'_Test_Perm_Comp')


    print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/Tagger_Vars/Test_Perm .\n' % (jobid, args.analysis, args.sample)




##############################################################################################

if args.plot == "MP_Tagger":

    hist_dir = [
        ('Partial_Merge_Matched_Perm_Plots/', '_Partial_Merge_Matched_Perm', 'Merged_Matched_Perm'),
        ('Unmerged_Matched_Perm_Plots/', '_Unmerged_Matched_Perm', 'Unmerged_Matched_Perm')
    ]

    for dir_name, hist_name, sdir_name in hist_dir:

        s_dir = 'Tagger_Vars/%s/Mass' % sdir_name
        plotter.set_subdir(s_dir)

        Matched_Perm_Mass_3J = [
            ('3J_BHadWJet_Mass', 'M_{b_{h}+wjet} [GeV]', 'black'),
            ('3J_BLepWJet_Mass', 'M_{b_{l}+wjet} [GeV]', 'black')
        ]
    
        for var, xaxis, colors in Matched_Perm_Mass_3J:
            hist = asrootpy(myfile.Get(dir_name+var)).Clone()
            if hist.Integral() == 0:
                print var
                continue
    
            mean = hist.GetMean()
            rms = hist.GetRMS()
    
            plotter.set_histo_style(hist, color=colors, xtitle=xaxis)
            plotter.plot(hist, drawstyle='hist')
            box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
            box.Draw()
            plotter.save(var+hist_name)
    
        mp_dict = {'Mass' : [' Mass [GeV]', 'black'], 'CSV' : [' CSV', 'black'], 'CvsL' : [' CvsL', 'black'], 'Mult' : [' Multiplicity', 'black']}
    
        Matched_Perm_3J = [
            ('3J_BHad_', 'b_{h}', 'black'),
            ('3J_BLep_', 'b_{l}', 'blue'),
            ('3J_WJet_', 'wjet', 'red')
        ]
    
        for key in mp_dict:
            dict_name = mp_dict[key][0]
            dict_color = mp_dict[key][1]
    
            s_dir = 'Tagger_Vars/%s/%s' % (sdir_name, key)
            plotter.set_subdir(s_dir)

            to_draw = []
        
            for var, objects, colors in Matched_Perm_3J:
                hist = asrootpy(myfile.Get(dir_name+var+key)).Clone()
                if hist.Integral() == 0:
                    print var
                    continue
     
                mean = hist.GetMean()
                rms = hist.GetRMS()
        
                plotter.set_histo_style(hist, color=dict_color, xtitle=objects+dict_name)
                plotter.plot(hist, drawstyle='hist')
                box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
                box.Draw()
                plotter.save(var+key+hist_name)
    
                plotter.set_histo_style(hist, title=objects, color=colors, xtitle=dict_name)
                plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(hist)
    
    
            plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            plotter.save('3J_'+key+hist_name+'_Comp')


        print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/Tagger_Vars/%s .\n' % (jobid, args.analysis, args.sample, sdir_name)



##############################################################################################

if args.plot == "Best_Perm":

    # best perm categories pt, mass, eta, costh, mult

    permcat = ['RIGHT', 'MERGE_SWAP', 'MERGE', 'WRONG'] # perm category


    for cat in permcat:
        m_dir = 'Matched_Perms/Mass/%s' % cat
        plotter.set_subdir(m_dir)
    
        Invar_Mass = [
            ('3J_BHadWJet_Mass_%s' % cat, 'M_{b_{h}+wjet} [GeV]'),
            ('3J_BLepWJet_Mass_%s' % cat, 'M_{b_{l}+wjet} [GeV]')
        ]
        
        for var, xaxis in Invar_Mass:
            to_draw = [] 
            L_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+var+'_LCut')).Clone()
            if L_hist.Integral() == 0:
                print var 
                continue
        
            G_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+var+'_GCut')).Clone()
            if G_hist.Integral() == 0:
                print var 
                continue
        
            plotter.set_histo_style(L_hist, title='< 2', color='black')
            plotter.plot(L_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(L_hist)
    
            plotter.set_histo_style(G_hist, title='#geq 2', color='red')
            plotter.plot(G_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(G_hist)
    
            stack, norm_stack = stack_plots(to_draw)
        
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xaxis, drawstyle='hist')
            box = plotter.make_text_box('%s' % cat, position='NE')
            box.Draw()
            plotter.save(var+'_Stack')
    #    
    #        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    #        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)
    
    kinvar = {'Mass' : 'Mass [GeV]', 'Pt' : 'p_{T}', 'Eta' : '#eta', 'Costh' : 'Cos(#theta)', 'Mult' : 'Mult'}

    for kvar in kinvar:
        xlabel = kinvar[kvar]

        for cat in permcat:
            m_dir = 'Matched_Perms/%s/%s' % (kvar, cat)
            
            plotter.set_subdir(m_dir)
            
            VAR = [
                ('3J_BHad_%s_%s' % (kvar, cat), 'b_{h} %s' % xlabel ),
                ('3J_BLep_%s_%s' % (kvar, cat), 'b_{l} %s' % xlabel ),
                ('3J_WJet_%s_%s' % (kvar, cat), 'wjet %s' % xlabel )
            ]
            
            for var, xaxis in VAR:
                to_draw = [] 
                L_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+var+'_LCut')).Clone()
                if L_hist.Integral() == 0:
                    print var 
                    continue
            
                G_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+var+'_GCut')).Clone()
                if G_hist.Integral() == 0:
                    print var 
                    continue
            
                plotter.set_histo_style(L_hist, title='< 2', color='black')
                plotter.plot(L_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(L_hist)
        
                plotter.set_histo_style(G_hist, title='#geq 2', color='red')
                plotter.plot(G_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(G_hist)
        
                stack, norm_stack = stack_plots(to_draw)
            
                plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xaxis, drawstyle='hist')
                box = plotter.make_text_box('%s' % cat, position='NE')
                box.Draw()
                plotter.save(var+'_Stack')
        #    
        #        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        #        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)
    
    
    
    
        print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/Matched_Perms/%s .\n' % (jobid, args.analysis, args.sample, kvar)



##############################################################################################


if args.plot == "Best_Perm_Combined":

    # best perm categories pt, mass, eta, costh, mult

#        kinvar = ['Mass', 'Pt', 'Eta', 'Costh', 'Mult'] # kinematic variables
    permcat = {'RIGHT' : 'black', 'MERGE_SWAP': 'red', 'MERGE' : 'green', 'WRONG' : 'blue'} # perm category
    Cut = {'LCut' : '< 2', 'GCut' : '#geq 2'}

    for cut in Cut:
        m_dir = 'Matched_Perms/%s/Mass' % cut
        plotter.set_subdir(m_dir)
        
        bhad_to_draw = [] 
        blep_to_draw = [] 
            
        for cat in permcat:
            Invar_Mass = [
                ('3J_BHadWJet_Mass_%s_%s' % (cat, cut), '3J_BLepWJet_Mass_%s_%s' % (cat,cut))
            ]

            for bhad, blep in Invar_Mass:
                had_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+bhad)).Clone()
                if had_hist.Integral() == 0:
                    print bhad 
                    continue
            
                plotter.set_histo_style(had_hist, title=cat, color=permcat[cat])
                plotter.plot(had_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                bhad_to_draw.append(had_hist)
    
                lep_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+blep)).Clone()
                if lep_hist.Integral() == 0:
                    print blep 
                    continue
            
                plotter.set_histo_style(lep_hist, title=cat, color=permcat[cat])
                plotter.plot(lep_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                blep_to_draw.append(lep_hist)
    
        had_stack, norm_had_stack = stack_plots(bhad_to_draw)
        
        plotter.plot(had_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='M_{b_{h}+wjet} [GeV]', drawstyle='hist')
        box = plotter.make_text_box('Perm Value %s' % Cut[cut], position='NE')
        box.Draw()
        plotter.save('3J_BHadWJet_Mass_%s_Stack' % cut)

        lep_stack, norm_lep_stack = stack_plots(blep_to_draw)
        
        plotter.plot(lep_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='M_{b_{l}+wjet} [GeV]', drawstyle='hist')
        box = plotter.make_text_box('Perm Value %s' % Cut[cut], position='NE')
        box.Draw()
        plotter.save('3J_BLepWJet_Mass_%s_Stack' % cut)

        #    
        #        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        #        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)

    
    kinvar = {'Mass' : 'Mass [GeV]', 'Pt' : 'p_{T}', 'Eta' : '#eta', 'Costh' : 'Cos(#theta)', 'Mult' : 'Mult'} #kinematic variables
#    kinvar = {'Mass' : 'Mass [GeV]'} #kinematic variables
    permcat = {'RIGHT' : 'black', 'MERGE_SWAP': 'red', 'MERGE' : 'green', 'WRONG' : 'blue'} # perm category
    Cut = {'LCut' : '< 2', 'GCut' : '#geq 2'}

    for cut in Cut:
        for kvar in kinvar:
            xlabel = kinvar[kvar]
            m_dir = 'Matched_Perms/%s/%s' % (cut, kvar)
            plotter.set_subdir(m_dir)
            
            bhad_to_draw = [] 
            blep_to_draw = [] 
            wjet_to_draw = [] 
                
            for cat in permcat:
                VAR = [
                    ('3J_BHad_%s_%s_%s' % (kvar, cat, cut), '3J_BLep_%s_%s_%s' % (kvar,cat,cut), '3J_WJet_%s_%s_%s' % (kvar, cat, cut))
                ]
    
                for bhad, blep, wjet in VAR:
                    had_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+bhad)).Clone()
                    if had_hist.Integral() == 0:
                        print bhad 
                        continue
                
                    plotter.set_histo_style(had_hist, title=cat, color=permcat[cat])
                    plotter.plot(had_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                    bhad_to_draw.append(had_hist)
        
                    lep_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+blep)).Clone()
                    if lep_hist.Integral() == 0:
                        print blep 
                        continue
                
                    plotter.set_histo_style(lep_hist, title=cat, color=permcat[cat])
                    plotter.plot(lep_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                    blep_to_draw.append(lep_hist)
        
                    w_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+wjet)).Clone()
                    if w_hist.Integral() == 0:
                        print wjet 
                        continue
                
                    plotter.set_histo_style(w_hist, title=cat, color=permcat[cat])
                    plotter.plot(w_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                    wjet_to_draw.append(w_hist)
        
            had_stack, norm_had_stack = stack_plots(bhad_to_draw)
            
            plotter.plot(had_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='b_{h} %s' % xlabel , drawstyle='hist')
            box = plotter.make_text_box('Perm Value %s' % Cut[cut], position='NE')
            box.Draw()
            plotter.save('3J_BHad_%s_%s_Stack' % (kvar,cut))
    
            lep_stack, norm_lep_stack = stack_plots(blep_to_draw)
            
            plotter.plot(lep_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='b_{l} %s' % xlabel, drawstyle='hist')
            box = plotter.make_text_box('Perm Value %s' % Cut[cut], position='NE')
            box.Draw()
            plotter.save('3J_BLep_%s_%s_Stack' % (kvar, cut))
    
            w_stack, norm_w_stack = stack_plots(wjet_to_draw)
            
            plotter.plot(w_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='wjet %s' % xlabel, drawstyle='hist')
            box = plotter.make_text_box('Perm Value %s' % Cut[cut], position='NE')
            box.Draw()
            plotter.save('3J_WJet_%s_%s_Stack' % (kvar, cut))
    
            #    
            #        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            #        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)

    
    
    
    
            print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/Matched_Perms/%s .\n' % (jobid, args.analysis, args.sample, kvar)


##############################################################################################


#if args.plot == "Best_Perm_Combined_Cut":
#
#    # best perm categories pt, mass, eta, costh, mult
#
##        kinvar = ['Mass', 'Pt', 'Eta', 'Costh', 'Mult'] # kinematic variables
#    permcat = {'RIGHT' : 'black', 'MERGE_SWAP': 'red', 'MERGE' : 'green', 'WRONG' : 'blue'} # perm category
#    Cut = {'LCut' : ['< 2', 'black'], 'GCut' : ['#geq 2', 'red']}
#
#    for cut in Cut:
#        m_dir = 'Matched_Perms/Comb_Cut/Mass'
#        plotter.set_subdir(m_dir)
#        
#        titles = Cut[cut][0]
#        colors = Cut[cut][1]
#
#        bhad_to_draw = [] 
##        blep_to_draw = [] 
#            
##        for cat in permcat:
#        Invar_Mass = [
#            ('3J_BHadWJet_Mass_%s' % cut, '3J_BLepWJet_Mass_%s' % cut)
#        ]
#
#        for bhad, blep in Invar_Mass:
#            had_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+bhad)).Clone()
#            if had_hist.Integral() == 0:
#                print bhad 
#                continue
#        
#            plotter.set_histo_style(had_hist, title=titles, color=colors)
#            plotter.plot(had_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#            bhad_to_draw.append(had_hist)
#    
##            lep_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+blep)).Clone()
##            if lep_hist.Integral() == 0:
##                print blep 
##                continue
#        
##            plotter.set_histo_style(lep_hist, title=titles, color=colors)
##            plotter.plot(lep_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
##            blep_to_draw.append(lep_hist)
#    
#        had_stack, norm_had_stack = stack_plots(bhad_to_draw)
#        
#        plotter.plot(had_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='M_{b_{h}+wjet} [GeV]', drawstyle='hist')
##        box = plotter.make_text_box('Perm Value %s' % title, position='NE')
##        box.Draw()
#        plotter.save('3J_BHadWJet_Mass_Combined_Cut_Stack')
#
##        lep_stack, norm_lep_stack = stack_plots(blep_to_draw)
##        
##        plotter.plot(lep_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='M_{b_{l}+wjet} [GeV]', drawstyle='hist')
##        box = plotter.make_text_box('Perm Value %s' % Cut[cut], position='NE')
##        box.Draw()
##        plotter.save('3J_BLepWJet_Mass_%s_Stack' % cut)
#
#        #    
#        #        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
#        #        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)
#
#    
##    kinvar = {'Mass' : 'Mass [GeV]', 'Pt' : 'p_{T}', 'Eta' : '#eta', 'Costh' : 'Cos(#theta)', 'Mult' : 'Mult'} #kinematic variables
###    kinvar = {'Mass' : 'Mass [GeV]'} #kinematic variables
##    permcat = {'RIGHT' : 'black', 'MERGE_SWAP': 'red', 'MERGE' : 'green', 'WRONG' : 'blue'} # perm category
##    Cut = {'LCut' : '< 2', 'GCut' : '#geq 2'}
##
##    for cut in Cut:
##        for kvar in kinvar:
##            xlabel = kinvar[kvar]
##            m_dir = 'Matched_Perms/%s/%s' % (cut, kvar)
##            plotter.set_subdir(m_dir)
##            
##            bhad_to_draw = [] 
##            blep_to_draw = [] 
##            wjet_to_draw = [] 
##                
##            for cat in permcat:
##                VAR = [
##                    ('3J_BHad_%s_%s_%s' % (kvar, cat, cut), '3J_BLep_%s_%s_%s' % (kvar,cat,cut), '3J_WJet_%s_%s_%s' % (kvar, cat, cut))
##                ]
##    
##                for bhad, blep, wjet in VAR:
##                    had_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+bhad)).Clone()
##                    if had_hist.Integral() == 0:
##                        print bhad 
##                        continue
##                
##                    plotter.set_histo_style(had_hist, title=cat, color=permcat[cat])
##                    plotter.plot(had_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
##                    bhad_to_draw.append(had_hist)
##        
##                    lep_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+blep)).Clone()
##                    if lep_hist.Integral() == 0:
##                        print blep 
##                        continue
##                
##                    plotter.set_histo_style(lep_hist, title=cat, color=permcat[cat])
##                    plotter.plot(lep_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
##                    blep_to_draw.append(lep_hist)
##        
##                    w_hist = asrootpy(myfile.Get('Best_Perm_Plots/'+wjet)).Clone()
##                    if w_hist.Integral() == 0:
##                        print wjet 
##                        continue
##                
##                    plotter.set_histo_style(w_hist, title=cat, color=permcat[cat])
##                    plotter.plot(w_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
##                    wjet_to_draw.append(w_hist)
##        
##            had_stack, norm_had_stack = stack_plots(bhad_to_draw)
##            
##            plotter.plot(had_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='b_{h} %s' % xlabel , drawstyle='hist')
##            box = plotter.make_text_box('Perm Value %s' % Cut[cut], position='NE')
##            box.Draw()
##            plotter.save('3J_BHad_%s_%s_Stack' % (kvar,cut))
##    
##            lep_stack, norm_lep_stack = stack_plots(blep_to_draw)
##            
##            plotter.plot(lep_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='b_{l} %s' % xlabel, drawstyle='hist')
##            box = plotter.make_text_box('Perm Value %s' % Cut[cut], position='NE')
##            box.Draw()
##            plotter.save('3J_BLep_%s_%s_Stack' % (kvar, cut))
##    
##            w_stack, norm_w_stack = stack_plots(wjet_to_draw)
##            
##            plotter.plot(w_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='wjet %s' % xlabel, drawstyle='hist')
##            box = plotter.make_text_box('Perm Value %s' % Cut[cut], position='NE')
##            box.Draw()
##            plotter.save('3J_WJet_%s_%s_Stack' % (kvar, cut))
##    
##            #    
##            #        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
##            #        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)
##
##    
##    
##    
##    
##            print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s/Matched_Perms/%s .\n' % (jobid, args.analysis, args.sample, kvar)



##############################################
#print('cp ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/jet_perm_disc/%s/%s/%s .' % (jobid, args.analysis, args.sample))

