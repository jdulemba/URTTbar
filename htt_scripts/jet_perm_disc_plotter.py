'''
Jet_Perm_Disc Analyzer Plotter macro
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os, glob, sys
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
import rootpy.plotting as plotting
from rootpy.plotting import views, Graph, Hist
from argparse import ArgumentParser
import URAnalysis.Utilities.roottools as urtools
import ROOT
from URAnalysis.Utilities.tables import latex_table
import argparse
#import matplotlib.pyplot as plt
from collections import OrderedDict
import functions as fncts


analyzer = 'jet_perm_disc'

parser = argparse.ArgumentParser(description='Create plots using files from jet_perm_disc.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000).')
parser.add_argument('plot', help='Choose type of plots to generate:\n Discriminants(Discs, Matched_Best_Perm)\n Tagger_Vars(TP_Tagger, MP_Tagger)\n Everything')
#parser.add_argument('plot', help='Choose type of plots to generate (Had_Comp, Delta_Plots, Discs, NS_Chi, Combine_Discs, TP_Tagger, MP_Tagger, Best_Perm, Best_Perm_Combined).')
args = parser.parse_args()


##### check analysis type
if not (args.analysis == 'Test' or args.analysis == 'Full'):
    print "You chose %s as your analysis type.\n You must choose Full or Test!" % args.analysis
    sys.exit()


##### check sample type
results_files = []
if args.analysis == "Test":
    for f in os.listdir('../'):
        if not '%s.test.root' % analyzer in f:
            continue
        results_files.append(f.replace(".root", ""))

    if not '%s.%s.test' % (args.sample, analyzer) in results_files:
        print "You chose %s as your sample file.\nYou must choose from the %s files in URTTbar!" % (args.sample, analyzer)
        sys.exit()

if args.analysis == "Full":
    for f in os.listdir('../results/%s/%s' % (jobid, analyzer)):
        if not '.root' in f:
            continue
        results_files.append(f.replace(".root", ""))

    if not args.sample in results_files:
        print "You chose %s as your sample file.\nYou must choose from the files in results/%s/%sJ!" % (args.sample, jobid, analyzer)
        sys.exit()


##### check plot type
if not (args.plot == "Everything" or args.plot == "Discriminants" or args.plot == "Discs"\
        or args.plot == "Matched_Best_Perm" or args.plot == "Tagger_Vars" or args.plot == "TP_Tagger"\
        or args.plot == "MP_Tagger"):
    print "You chose %s as your plot type.\nYou must choose from the help list!"
    sys.exit()


print( 'Analysis: %s\nSample: %s\nPlot: %s' % (args.analysis, args.sample, args.plot) )


if args.analysis == "Test":
    myfile = root_open('../%s.%s.test.root' % (args.sample, analyzer), 'read')
    normfile = views.NormalizeView(root_open('../%s.%s.test.root' % (args.sample, analyzer), 'read'))#normalized file

elif args.analysis == "Full":
	myfile = root_open('../results/%s/%s/%s.root' % (jobid, analyzer, args.sample), 'read')
	normfile = views.NormalizeView(root_open('../results/%s/%s/%s.root' % (jobid, analyzer, args.sample), 'read'))

plotter = BasePlotter(
	'plots/%s/%s/%s/%s' % (analyzer, jobid, args.analysis, args.sample),
	defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}}
	#defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
)

#def stack_plots(lists):
#    lists.sort(key=lambda x: x.Integral())
#    total = 0
##    ratio_hists = []
#    Stack_hists = []
#
#    for i in lists:
#        total += i
#
#    for i in lists:
##        ratio_hists.append(i/total)
#       # i.SetFillStyle(1001)
#        Stack_hists.append(i)
#
#    stack = plotting.HistStack()
#    norm_stack = plotting.HistStack()
#    for i in Stack_hists:
#        stack.Add(i)
#        norm_stack.Add(i/total)
#
#    return stack, norm_stack
   

def plot_dir(plot):
    print '\ncp -r /uscms/home/jdulemba/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/%s/%s/%s/%s/%s .\n' % (analyzer, jobid, args.analysis, args.sample, plot)
 

##### Global Var. Definitions ####
defcol = 'black' # default color
defyax = 'A.U.' # yaxis title
if args.sample == "ttJetsM0":
	mass_min = 0
	mass_max = 2000
if args.sample == "ttJetsM700":
	mass_min = 700
	mass_max = 1000
if args.sample == "ttJetsM1000":
	mass_min = 1000
	mass_max = 2000


if 'ttJets' in args.sample:
    decay = 'SM t#bar t'

if 'AtoTT' in args.sample:
    decay = 'A->t#bar t'

if 'HtoTT' in args.sample:
    decay = 'H->t#bar t'


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

#if args.plot == "Had_Comp":
#### compares relative frequencies of:
#### 1.   three jets on hadronic side being resolved
#### 2.   only bhad and one w jet merging
#### 3.   only jets from w merging
#### 4.   all three hadronic jets merging
#### 5.   non-reconstructable events
#
#
#    ### Events with 3 jets
#	Hadronic_Events_3J = [
#		('3J_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
#		('3J_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
#		('3J_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
#		('3J_Merged_THad_vs_mttbar', 'All 3 Merged', 'green'),
#		('3J_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
#	]
#
#	for name, title, col in DRvals:
#
#		to_draw = []
#		Add_Ratio_Hists = []
#		Ratio_Hists = []
#	
#		for var, legends, colors in Hadronic_Events_3J:
#			hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
#			plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
#			to_draw.append(hist)
#	
#		total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]
#	
#		### plot log of event occurrences
#		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
#		plotter.save('3J_Hadronic_Events_vs_mttbar_log'+name)
#	
#	
#		### plot event ratios
#		NonReco_Ratio = to_draw[4]/total
#		Ratio_Hists.append(NonReco_Ratio)
#	
#		Resolved_Ratio = to_draw[0]/total
#		Ratio_Hists.append(Resolved_Ratio)
#	
#		M_BHadW_Ratio = to_draw[1]/total
#		Ratio_Hists.append(M_BHadW_Ratio)
#	
#		M_WJet_Ratio = to_draw[2]/total
#		Ratio_Hists.append(M_WJet_Ratio)
#	
#		M_THad_Ratio = to_draw[3]/total
#		Ratio_Hists.append(M_THad_Ratio)
#	
#		plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.15), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
#		plotter.save('3J_Hadronic_Events_vs_mttbar_ratio'+name)
#	
#	
#		### create stacked hists of ratios
#		stack_hists = []
#		for i in range(len(to_draw)):
#			stack_hists.append(to_draw[i].Integral())
#		stack_hists.sort()
#	
#		for i in range(len(stack_hists)):
#			for j in range(len(to_draw)):
#				if to_draw[j].Integral() == stack_hists[i]:
#					to_draw[j].SetFillStyle(1001)
#					Add_Ratio_Hists.append(to_draw[j])
#	
#		stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total, Add_Ratio_Hists[4]/total)
#		plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title)
#		plotter.save('3J_Hadronic_Events_vs_mttbar_stack'+name)
#
#
#
#   ### Events with 4 jets
#	Hadronic_Events_4J = [
#		('4J_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
#		('4J_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
#		('4J_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
#		('4J_Merged_THad_vs_mttbar', 'All 3 Merged', 'green'),
#		('4J_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
#	]
#
#	for name, title, col in DRvals:
#
#		to_draw = []
#		Add_Ratio_Hists = []
#		Ratio_Hists = []
#	
#		for var, legends, colors in Hadronic_Events_4J:
#			hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
#			plotter.set_histo_style(hist, title=legends, color=colors, drawstyle='hist', linestyle='solid')
#			to_draw.append(hist)
#	
#		total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]
#	
#		### plot log of event occurrences
#		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
#		plotter.save('4J_Hadronic_Events_vs_mttbar_log'+name)
#	
#	
#		### plot event ratios
#		NonReco_Ratio = to_draw[4]/total
#		Ratio_Hists.append(NonReco_Ratio)
#	
#		Resolved_Ratio = to_draw[0]/total
#		Ratio_Hists.append(Resolved_Ratio)
#	
#		M_BHadW_Ratio = to_draw[1]/total
#		Ratio_Hists.append(M_BHadW_Ratio)
#	
#		M_WJet_Ratio = to_draw[2]/total
#		Ratio_Hists.append(M_WJet_Ratio)
#	
#		M_THad_Ratio = to_draw[3]/total
#		Ratio_Hists.append(M_THad_Ratio)
#	
#		plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,0.85), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
#		plotter.save('4J_Hadronic_Events_vs_mttbar_ratio'+name)
#	
#	
#		### create stacked hists of ratios
#		stack_hists = []
#		for i in range(len(to_draw)):
#			stack_hists.append(to_draw[i].Integral())
#		stack_hists.sort()
#	
#		for i in range(len(stack_hists)):
#			for j in range(len(to_draw)):
#				if to_draw[j].Integral() == stack_hists[i]:
#					to_draw[j].SetFillStyle(1001)
#					Add_Ratio_Hists.append(to_draw[j])
#	
#		stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total, Add_Ratio_Hists[4]/total)
#		plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0, 1.3), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
#		plotter.save('4J_Hadronic_Events_vs_mttbar_stack'+name)
#
#
#
#    ##### Events with 5 or more jets
#	Hadronic_Events_5PJ = [
#		('5PJ_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
#		('5PJ_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
#		('5PJ_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
#		('5PJ_Merged_THad_vs_mttbar', 'All 3 Merged', 'green'),
#		('5PJ_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
#	]
#
#	for name, title, col in DRvals:
#
#		to_draw = []
#		Add_Ratio_Hists = []
#		Ratio_Hists = []
#
#		for var, legends, colors in Hadronic_Events_5PJ:
#			hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
#			plotter.set_histo_style(hist, title=legends, color=colors, drawstyle='hist', linestyle='solid')
#			to_draw.append(hist)
#
#		total = to_draw[0]+to_draw[1]+to_draw[2]+to_draw[3]+to_draw[4]
#
#		### plot log of event occurrences
#		plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
#		plotter.save('5PJ_Hadronic_Events_vs_mttbar_log'+name)
#
#
#		### plot event ratios
#		NonReco_Ratio = to_draw[4]/total
#		Ratio_Hists.append(NonReco_Ratio)
#
#		Resolved_Ratio = to_draw[0]/total
#		Ratio_Hists.append(Resolved_Ratio)
#
#		M_BHadW_Ratio = to_draw[1]/total
#		Ratio_Hists.append(M_BHadW_Ratio)
#
#		M_WJet_Ratio = to_draw[2]/total
#		Ratio_Hists.append(M_WJet_Ratio)
#
#		M_THad_Ratio = to_draw[3]/total
#		Ratio_Hists.append(M_THad_Ratio)
#
#		plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,0.75), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
#		plotter.save('5PJ_Hadronic_Events_vs_mttbar_ratio'+name)
#
#
#		### create stacked hists of ratios
#		stack_hists = []
#		for i in range(len(to_draw)):
#			stack_hists.append(to_draw[i].Integral())
#		stack_hists.sort()
#
#		for i in range(len(stack_hists)):
#			for j in range(len(to_draw)):
#				if to_draw[j].Integral() == stack_hists[i]:
#					to_draw[j].SetFillStyle(1001)
#					Add_Ratio_Hists.append(to_draw[j])
#
#		stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total, Add_Ratio_Hists[4]/total)
#		plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0, 1.3), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
#		plotter.save('5PJ_Hadronic_Events_vs_mttbar_stack'+name)
#
#
#
###############################################################################################
#
#if args.plot == "Delta_Plots":
#	# Delta (Gen - Reco) Plots
#
#	Sig_Background_mass = [
#		('3J_Preselection_Wrong_Perms_Delta_mass', 'Preselection Wrong', 'red', 'solid'),
#		('3J_Square_Cut_Wrong_Perms_Delta_mass', ' Signal Wrong', 'red', 'dashed'),
#		('3J_Preselection_Correct_Perms_Delta_mass', 'Preselection Correct', 'blue', 'solid'),
#		('3J_Square_Cut_Correct_Perms_Delta_mass', 'Signal Correct', 'blue', 'dashed')
#	]
#
#	to_draw = []
#
#	for var, legends, colors, linestyles in Sig_Background_mass:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=legends, color=colors, linestyle=linestyles)
#		to_draw.append(hist)
#
#	Sig_correct_mass = asrootpy(myfile.Get('Delta_Plots/3J_Square_Cut_Correct_Perms_Delta_mass')).Clone()
#	Presel_correct_mass = asrootpy(myfile.Get('Delta_Plots/3J_Preselection_Correct_Perms_Delta_mass')).Clone()
#	Sig_wrong_mass = asrootpy(myfile.Get('Delta_Plots/3J_Square_Cut_Wrong_Perms_Delta_mass')).Clone()
#	Presel_wrong_mass = asrootpy(myfile.Get('Delta_Plots/3J_Preselection_Wrong_Perms_Delta_mass')).Clone()
#	print 'Gen-Reco Mass:'
#	print 'Correct Preselection Purity = ', format(Presel_correct_mass.Integral()/(Presel_correct_mass.Integral()+Presel_wrong_mass.Integral()), '.4f')
#	print 'Wrong Preselection Purity = ', format(Presel_wrong_mass.Integral()/(Presel_correct_mass.Integral()+Presel_wrong_mass.Integral()), '.4f')
#	print 'Correct Signal Purity = ',  format(Sig_correct_mass.Integral()/(Sig_correct_mass.Integral()+Sig_wrong_mass.Integral()), '.4f')
#	print 'Wrong Signal Purity = ',  format(Sig_wrong_mass.Integral()/(Sig_correct_mass.Integral()+Sig_wrong_mass.Integral()), '.4f')
#	print 'Correct Efficiency = ', format(Sig_correct_mass.Integral()/Presel_correct_mass.Integral(), '.4f')
#	print 'Wrong Efficiency = ', format(Sig_wrong_mass.Integral()/Presel_wrong_mass.Integral(), '.4f')
#
#
#	#### plot log
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l',logy=True, x_range=(-1000, 200), y_range=(10**(-1),10**6), xtitle='#Delta m_{t_{h}} [GeV] (Gen - Reco)', drawstyle='hist')
#	plotter.save('3J_Sig_Background_Delta_mass_log')
#	
#	#### plot normal
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', x_range=(-1000, 200), xtitle='#Delta m_{t_{h}} [GeV] (Gen - Reco)', drawstyle='hist')
#	plotter.save('3J_Sig_Background_Delta_mass')
#	
#
#	Sig_Background_costh = [
#		('3J_Preselection_Wrong_Perms_Delta_costh', 'Preselection Wrong', 'red', 'solid'),
#		('3J_Square_Cut_Wrong_Perms_Delta_costh', 'Signal Wrong', 'red', 'dashed'),
#		('3J_Preselection_Correct_Perms_Delta_costh', 'Preselection Correct', 'blue', 'solid'),
#		('3J_Square_Cut_Correct_Perms_Delta_costh', 'Signal Correct', 'blue', 'dashed')
#	]
#
#	to_draw = []
#
#	for var, legends, colors, linestyles in Sig_Background_costh:
#		hist = asrootpy(myfile.Get('Delta_Plots/'+var)).Clone()
#		plotter.set_histo_style(hist, title=legends, color=colors, linestyle=linestyles)
#		to_draw.append(hist)
#
#	Sig_correct_costh = asrootpy(myfile.Get('Delta_Plots/3J_Square_Cut_Correct_Perms_Delta_costh')).Clone()
#	Presel_correct_costh = asrootpy(myfile.Get('Delta_Plots/3J_Preselection_Correct_Perms_Delta_costh')).Clone()
#	Sig_wrong_costh = asrootpy(myfile.Get('Delta_Plots/3J_Square_Cut_Wrong_Perms_Delta_costh')).Clone()
#	Presel_wrong_costh = asrootpy(myfile.Get('Delta_Plots/3J_Preselection_Wrong_Perms_Delta_costh')).Clone()
#	print 'Gen-Reco Cos(#theta):'
#	print 'Correct Presel/ Total Preselection = ', format(Presel_correct_costh.Integral()/(Presel_correct_costh.Integral()+Presel_wrong_costh.Integral()), '.4f')
#	print 'Wrong Presel/ Total Preselection = ', format(Presel_wrong_costh.Integral()/(Presel_correct_costh.Integral()+Presel_wrong_costh.Integral()), '.4f')
#	print 'Correct Signal/Total Signal = ',  format(Sig_correct_costh.Integral()/(Sig_correct_costh.Integral()+Sig_wrong_costh.Integral()), '.4f')
#	print 'Wrong Signal/Total Signal = ',  format(Sig_wrong_costh.Integral()/(Sig_correct_costh.Integral()+Sig_wrong_costh.Integral()), '.4f')
#
#	### plot log
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', logy=True, x_range=(-2, 2), y_range=(10**(-1),10**7), xtitle='#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', drawstyle='hist')
#	plotter.save('3J_Sig_Background_Delta_costh_log')
#	
#	### plot normal 
#	plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', x_range=(-2, 2), xtitle='#Delta Cos(#theta_{t_{h}}) (Gen - Reco)', drawstyle='hist')
#	plotter.save('3J_Sig_Background_Delta_costh')
	


##############################################################################################

def Discs():
    plots = 'Discriminants'

    disc_dict = {'Massdisc' : ['black', 'Mass_Disc', '#lambda_{mass}'], 'NSdisc' : ['black', 'NS_Disc', '#lambda_{NS}'], 'Totaldisc' : ['black', 'Total_Disc', '#lambda_{comb}'], 'NSchi' : ['black', 'NS_Chi', '#chi^{2}']}
    for disc in disc_dict:
        disc_col = disc_dict[disc][0]
        disc_dir = disc_dict[disc][1]
        disc_xaxis = disc_dict[disc][2]

        #sub_dir = disc_dir
        plotter.set_subdir('Discriminants/%s' % disc_dir)

        #### 3J Disc  values
        Likelihood_All_3J = [
            ('3J_Event_All_', 'All %s 3 jets' % disc_xaxis, '3J_Event_All_'),
            ('3J_Event_Lowest_', 'Lowest %s 3 jets' % disc_xaxis, '3J_Event_Lowest_'),
            ('3J_Event_Best_Perm_', 'Best Perm %s 3 jets' % disc_xaxis, '3J_Event_Best_Perm_')
        ]
       
        for var, xaxis, names in Likelihood_All_3J:
            hist = asrootpy(myfile.Get('Likelihood_Plots/'+var+disc)).Clone()
            if hist.Integral() == 0:
                print var
                continue
    
            mean = hist.GetMean()
            rms = hist.GetRMS()
    
            plotter.set_histo_style(hist, color=disc_col, xtitle=xaxis, ytitle=defyax)
            plotter.plot(hist, drawstyle='hist')
            box = plotter.make_text_box('Mean = %.3f\nRMS = %.3f' % (mean,rms), position='NE')
            box.Draw()
            plotter.save(names+disc)
    

    ### Plots for lowest likelihood values from event

        #m_dir = 'Matched_Perms/%s' % disc_dir
        plotter.set_subdir('Discriminants/%s/Best_Perm_Categories' % disc_dir)

            ### 3 jets

        Matched_Perm_Likelihood_3J = [
            ('3J_Best_Perm_RIGHT_', 'RIGHT', 'red', '%s 3 jets' % disc_xaxis, 3345),
            ('3J_Best_Perm_MERGE_SWAP_', 'MERGE_SWAP', 'orange', '%s 3 jets' % disc_xaxis, 0),
            ('3J_Best_Perm_MERGE_', 'MERGE', 'green', '%s 3 jets' % disc_xaxis, 3006),
            ('3J_Best_Perm_WRONG_', 'WRONG', 'blue', '%s 3 jets' % disc_xaxis, 3354)
        ]
       
        to_draw = [] 
        for var, legends, colors, xaxis, fill_style in Matched_Perm_Likelihood_3J:
            hist = asrootpy(myfile.Get('Likelihood_Plots/'+var+disc)).Clone()
            if hist.Integral() == 0:
                print var
                continue
    
            mean = hist.GetMean()
            rms = hist.GetRMS()
    
            plotter.set_histo_style(hist, title=legends, color=disc_col, xtitle=xaxis, ytitle=defyax)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            box = plotter.make_text_box(decay+'\nMean = %.3f\nRMS = %.3f' % (mean,rms), position='NE')
            box.Draw()
            plotter.save(var+disc+'_Categories')

            plotter.set_histo_style(hist, title=legends, color=colors, xtitle=xaxis, ytitle=defyax)
            hist.SetFillStyle(fill_style)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(hist)
    
    
        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        box = plotter.make_text_box(decay, position='NE')
        box.Draw()
        plotter.save('3J_Best_Perm_Categories_'+disc)
    
            
        stack, norm_stack, ratio = fncts.stack_plots(to_draw)
    
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), xtitle=xaxis, ytitle=defyax, legendstyle='l', drawstyle='hist')
        box = plotter.make_text_box(decay, position='NE')
        box.Draw()
        plotter.save('3J_Best_Perm_Categories_Stack_'+disc)
    
        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), xtitle=xaxis, ytitle=defyax, legendstyle='l', drawstyle='hist')
        box = plotter.make_text_box(decay, position='NE')
        box.Draw()
        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)

    
    plot_dir(plots)


##############################################################################################

def Combine_Discs():
    plots = 'Discriminants'
    plotter.set_subdir('Discriminants')

    disc_dict = {'Totaldisc' : ['black', '#lambda_{comb}'], 'Massdisc' : ['red', '#lambda_{mass}'], 'NSdisc' : ['blue', '#lambda_{NS}']}
    Evt_type = ['All', 'Lowest', 'Best_Perm']

    for evt in Evt_type:
        to_draw = []

        for disc in OrderedDict(sorted(disc_dict.items(), key=lambda t: t[1])): # makes sure total disc is first key
            disc_col = disc_dict[disc][0]
            disc_leg = disc_dict[disc][1]
            
            var = '3J_Event_%s_%s' % (evt, disc)
            hist = asrootpy(myfile.Get('Likelihood_Plots/'+var)).Clone()
            if hist.Integral() == 0:
                print var
                continue
    
            plotter.set_histo_style(hist, title=disc_leg, color=disc_col, xtitle='#lambda = -log(L)', ytitle=defyax)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(hist)

        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        plotter.save('3J_Event_%s_Disc_Comp' % evt)

    plot_dir(plots)


##############################################################################################

def TP_Tagger():
    plots = 'Tagger_Vars/Test_Perm'

    s_dir = 'Tagger_Vars/Test_Perm/Mass'
    plotter.set_subdir(s_dir)

    Test_Perm_Mass_3J = [
        ('3J_bj1wj1_Mass', 'm_{b_{1}+w_{1}} [GeV]', 'black'),
        ('3J_bj2wj1_Mass', 'm_{b_{2}+w_{1}} [GeV]', 'black')
    ]

    for var, xaxis, colors in Test_Perm_Mass_3J:
        hist = asrootpy(myfile.Get('Test_Perm_Plots/'+var)).Clone()
        if hist.Integral() == 0:
            print var
            continue

        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, color=colors, xtitle=xaxis, ytitle=defyax)
        plotter.plot(hist, drawstyle='hist')
        box = plotter.make_text_box('Mean = %.3f\nRMS = %.3f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(var+'_Test_Perm')

    tp_dict = {'Mass' : [' [GeV]', 'black'], 'MVA' : [' MVA', 'black'], 'CvsL' : [' CvsL', 'black'], 'Mult' : [' Multiplicity', 'black']}

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

            if key == 'Mass':
                plotter.set_histo_style(hist, color=dict_color, xtitle='m_{'+objects+'}'+dict_name, ytitle=defyax)
            else:
                plotter.set_histo_style(hist, color=dict_color, xtitle=objects+dict_name, ytitle=defyax)
            plotter.plot(hist, drawstyle='hist')
            box = plotter.make_text_box('Mean = %.3f\nRMS = %.3f' % (mean,rms), position='NE')
            box.Draw()
            plotter.save(var+key+'_Test_Perm')

            plotter.set_histo_style(hist, title=objects, color=colors, xtitle=dict_name, ytitle=defyax)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(hist)


        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        plotter.save('3J_'+key+'_Test_Perm_Comp')

    plot_dir(plots)


##############################################################################################

def MP_Tagger():

    hist_dir = [
        ('Partial_Merge_Matched_Perm_Plots/', '_Partial_Merge_Matched_Perm', 'Merged_Matched_Perm'),
        ('Unmerged_Matched_Perm_Plots/', '_Unmerged_Matched_Perm', 'Unmerged_Matched_Perm')
    ]

    for dir_name, hist_name, sdir_name in hist_dir:
        plots = 'Tagger_Vars/%s' % sdir_name

        s_dir = 'Tagger_Vars/%s/Mass' % sdir_name
        plotter.set_subdir(s_dir)

        Matched_Perm_Mass_3J = [
            ('3J_BHadWJet_Mass', 'm_{b_{h}+w_{j}} [GeV]', 'black'),
            ('3J_BLepWJet_Mass', 'm_{b_{l}+w_{j}} [GeV]', 'black')
        ]
    
        for var, xaxis, colors in Matched_Perm_Mass_3J:
            hist = asrootpy(myfile.Get(dir_name+var)).Clone()
            if hist.Integral() == 0:
                print var
                continue
    
            mean = hist.GetMean()
            rms = hist.GetRMS()
    
            plotter.set_histo_style(hist, color=colors, xtitle=xaxis, ytitle=defyax)
            plotter.plot(hist, drawstyle='hist')
            box = plotter.make_text_box('Mean = %.3f\nRMS = %.3f' % (mean,rms), position='NE')
            box.Draw()
            plotter.save(var+hist_name)
    
        mp_dict = {'Mass' : [' [GeV]', 'black'], 'MVA' : [' MVA', 'black'], 'CvsL' : [' CvsL', 'black'], 'Mult' : [' Multiplicity', 'black']}
    
        Matched_Perm_3J = [
            ('3J_BHad_', 'b_{h}', 'black'),
            ('3J_BLep_', 'b_{l}', 'blue'),
            ('3J_WJet_', 'w_{j}', 'red')
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
        
                if key == 'Mass':
                    plotter.set_histo_style(hist, color=dict_color, xtitle='m_{'+objects+'}'+dict_name, ytitle=defyax)
                else:
                    plotter.set_histo_style(hist, color=dict_color, xtitle=objects+dict_name, ytitle=defyax)
                plotter.plot(hist, drawstyle='hist')
                box = plotter.make_text_box('Mean = %.3f\nRMS = %.3f' % (mean,rms), position='NE')
                box.Draw()
                plotter.save(var+key+hist_name)
   
                if key == 'Mass': 
                    plotter.set_histo_style(hist, title=objects, color=colors, xtitle='Mass '+dict_name, ytitle=defyax)
                else:
                    plotter.set_histo_style(hist, title=objects, color=colors, xtitle=dict_name, ytitle=defyax)
                plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(hist)
    
    
            plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            plotter.save('3J_'+key+hist_name+'_Comp')

   
        plot_dir(plots) 


##############################################################################################

def Matched_Best_Perm():

    # best perm categories pt, mass, eta, costh, mult

    permcat = ['RIGHT', 'MERGE_SWAP', 'MERGE', 'WRONG'] # perm category


    for cat in permcat:
        m_dir = 'Discriminants/Matched_Best_Perms/Mass/%s' % cat
        plotter.set_subdir(m_dir)
    
        Invar_Mass = [
            ('3J_BHadWJet_Mass_%s' % cat, 'm_{b_{h}+w_{j}} [GeV]'),
            ('3J_BLepWJet_Mass_%s' % cat, 'm_{b_{l}+w_{j}} [GeV]')
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
    
            stack, norm_stack, ratio = fncts.stack_plots(to_draw)
        
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xaxis, ytitle=defyax, drawstyle='hist')
            box = plotter.make_text_box('%s' % cat, position='NE')
            box.Draw()
            plotter.save(var+'_Stack')
    #    
    #        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    #        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)
    
    kinvar = {'Mass' : ' ', 'Pt' : 'p_{T}', 'Eta' : '#eta', 'Costh' : 'Cos(#theta)', 'Mult' : 'Mult'}

    for kvar in kinvar:
        plots = 'Discriminants/Matched_Best_Perms/%s' % kvar

        xlabel = kinvar[kvar]

        for cat in permcat:
            m_dir = 'Discriminants/Matched_Best_Perms/%s/%s' % (kvar, cat)
            
            plotter.set_subdir(m_dir)
            
            VAR = [
                ('3J_BHad_%s_%s' % (kvar, cat), 'b_{h} %s' % xlabel ),
                ('3J_BLep_%s_%s' % (kvar, cat), 'b_{l} %s' % xlabel ),
                ('3J_WJet_%s_%s' % (kvar, cat), 'w_{j} %s' % xlabel )
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
        
                stack, norm_stack, ratio = fncts.stack_plots(to_draw)

                if kvar == 'Mass':
                    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{'+xaxis+'} [GeV]', ytitle=defyax, drawstyle='hist')
                else:
                    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xaxis, ytitle=defyax, drawstyle='hist')
                box = plotter.make_text_box('%s' % cat, position='NE')
                box.Draw()
                plotter.save(var+'_Stack')
        #    
        #        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        #        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)
    

        plot_dir(plots)    


##############################################################################################

def Matched_Best_Perm_Combined():

    # best perm categories pt, mass, eta, costh, mult

#        kinvar = ['Mass', 'Pt', 'Eta', 'Costh', 'Mult'] # kinematic variables
    permcat = {'RIGHT' : 'red', 'MERGE_SWAP': 'orange', 'MERGE' : 'green', 'WRONG' : 'blue'} # perm category
    Cut = {'LCut' : '< 2', 'GCut' : '#geq 2'}

    for cut in Cut:
        m_dir = 'Discriminants/Matched_Best_Perms/%s/Mass' % cut
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
    
        had_stack, norm_had_stack, ratio = fncts.stack_plots(bhad_to_draw)
        
        plotter.plot(had_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{b_{h}+w_{j}} [GeV]', ytitle=defyax, drawstyle='hist')
        box = plotter.make_text_box('#lambda_{comb} %s' % Cut[cut], position='NE')
        box.Draw()
        plotter.save('3J_BHadWJet_Mass_%s_Stack' % cut)

        lep_stack, norm_lep_stack, ratio = fncts.stack_plots(blep_to_draw)
        
        plotter.plot(lep_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{b_{l}+w_{j}} [GeV]', ytitle=defyax, drawstyle='hist')
        box = plotter.make_text_box('#lambda_{comb} %s' % Cut[cut], position='NE')
        box.Draw()
        plotter.save('3J_BLepWJet_Mass_%s_Stack' % cut)

        #    
        #        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        #        plotter.save('3J_Best_Perm_Categories_Stack_Norm_'+disc)

    
    kinvar = {'Mass' : '[GeV]', 'Pt' : 'p_{T}', 'Eta' : '#eta', 'Costh' : 'Cos(#theta)', 'Mult' : 'Mult'} #kinematic variables
#    kinvar = {'Mass' : 'Mass [GeV]'} #kinematic variables
    permcat = {'RIGHT' : 'red', 'MERGE_SWAP': 'orange', 'MERGE' : 'green', 'WRONG' : 'blue'} # perm category
    Cut = {'LCut' : '< 2', 'GCut' : '#geq 2'}

    for cut in Cut:
        for kvar in kinvar:
            plots = 'Discriminants/Matched_Best_Perms/%s' % kvar

            xlabel = kinvar[kvar]
            m_dir = 'Discriminants/Matched_Best_Perms/%s/%s' % (cut, kvar)
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
        
            had_stack, norm_had_stack, ratio = fncts.stack_plots(bhad_to_draw)

            if kvar == 'Mass':
                plotter.plot(had_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{b_{h}} %s' % xlabel, ytitle=defyax, drawstyle='hist')
            else:
                plotter.plot(had_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='b_{h} %s' % xlabel, ytitle=defyax, drawstyle='hist')
            box = plotter.make_text_box('#lambda_{comb} %s' % Cut[cut], position='NE')
            box.Draw()
            plotter.save('3J_BHad_%s_%s_Stack' % (kvar,cut))
    
            lep_stack, norm_lep_stack, ratio = fncts.stack_plots(blep_to_draw)

            if kvar == 'Mass':            
                plotter.plot(lep_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{b_{l}} %s' % xlabel, ytitle=defyax, drawstyle='hist')
            else:
                plotter.plot(lep_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='b_{l} %s' % xlabel, ytitle=defyax, drawstyle='hist')
            box = plotter.make_text_box('#lambda_{comb} %s' % Cut[cut], position='NE')
            box.Draw()
            plotter.save('3J_BLep_%s_%s_Stack' % (kvar, cut))
    
            w_stack, norm_w_stack, ratio = fncts.stack_plots(wjet_to_draw)
            
            if kvar == 'Mass':
                plotter.plot(w_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{w_{j}} %s' % xlabel, ytitle=defyax, drawstyle='hist')
            else:
                plotter.plot(w_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='w_{j} %s' % xlabel, ytitle=defyax, drawstyle='hist')
            box = plotter.make_text_box('#lambda_{comb} %s' % Cut[cut], position='NE')
            box.Draw()
            plotter.save('3J_WJet_%s_%s_Stack' % (kvar, cut))
    
   
            plot_dir(plots) 


###############################################################################################################################################################

if args.plot == "Discs":
    Discs()
    Combine_Discs()

if args.plot == "Matched_Best_Perm":
    Matched_Best_Perm()
    Matched_Best_Perm_Combined()

if args.plot == "Discriminants":
    Discs()
    Combine_Discs()
    Matched_Best_Perm()
    Matched_Best_Perm_Combined()

if args.plot == "TP_Tagger":
    TP_Tagger()

if args.plot == "MP_Tagger":
    MP_Tagger()

if args.plot == "Tagger_Vars":
    TP_Tagger()
    MP_Tagger()
   
if args.plot == "Everything":
    Discs()
    Combine_Discs()
    Matched_Best_Perm()
    Matched_Best_Perm_Combined()
    TP_Tagger()
    MP_Tagger()
 
