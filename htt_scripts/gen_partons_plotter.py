'''
Gen_Partons Analyzer Plotter macro
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os, glob, sys
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
from rootpy.plotting import views, Graph, Hist, Hist2D
from argparse import ArgumentParser
import URAnalysis.Utilities.roottools as urtools
import ROOT
from URAnalysis.Utilities.tables import latex_table
import argparse
import numpy as np
from URAnalysis.PlotTools.views.RebinView import RebinView
import rootpy.io as io


parser = argparse.ArgumentParser(description='Create plots using files from gen_partons.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, or Combined).')
parser.add_argument('plot', help='Choose type of plots to generate (DR, Kin_Vars, System, Merged_Fractions, 3J, 4J, 5PJ, AllJets, 3Partons, Had_Comp_DR).')
args = parser.parse_args()

if not (args.analysis == "Test" or args.analysis == "Full"):
    print "You chose %s as your analysis type.\n You must choose Full or Test!" % args.analysis
    sys.exit()

if not (args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000" or args.sample == "Combined"):
    print "You chose %s as your sample file.\n You must choose ttJetsM0, ttJetsM700, or ttJetsM1000!" % args.sample
    sys.exit()

if not (args.plot == "3J" or args.plot == "4J" or args.plot == "5PJ" or args.plot == "AllJets" or args.plot == "3Partons"\
    or args.plot == "Had_Comp_DR" or args.plot == "DR" or args.plot == "Kin_Vars"\
    or args.plot == "System" or args.plot == "Merged_Fractions"):
    print "You chose %s as your plot type.\nYou must choose from the help list!"
    sys.exit() 

if args.analysis == "Test":
    print( 'Your analysis type and sample file are %s and %s' % (args.analysis, args.sample) )

    if args.sample == "Combined":
        M0_file = root_open('../ttJetsM0.gen_partons.test.root', 'read')
        M0_normfile = views.NormalizeView(root_open('../ttJetsM0.gen_partons.test.root', 'read'))#normalized file
        M700_file = root_open('../ttJetsM700.gen_partons.test.root', 'read')
        M700_normfile = views.NormalizeView(root_open('../ttJetsM700.gen_partons.test.root', 'read'))#normalized file
        M1000_file = root_open('../ttJetsM1000.gen_partons.test.root', 'read')
        M1000_normfile = views.NormalizeView(root_open('../ttJetsM1000.gen_partons.test.root', 'read'))#normalized file

    else:
        myfile = root_open('../%s.gen_partons.test.root' % args.sample, 'read')
        normfile = views.NormalizeView(root_open('../%s.gen_partons.test.root' % args.sample, 'read'))#normalized file

    plotter = BasePlotter(
    	'plots/gen_partons/%s/%s/%s/%s' % (jobid, args.analysis, args.sample, args.plot),
    	defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}}
    )

elif args.analysis == "Full":
    print( 'Your analysis type and sample file are %s and %s' % (args.analysis, args.sample) )

    if args.sample == "Combined":
        ttJets_files = [
            glob.glob('../results/%s/gen_partons/ttJetsM0.root' % jobid),
            glob.glob('../results/%s/gen_partons/ttJetsM700.root' % jobid),
            glob.glob('../results/%s/gen_partons/ttJetsM1000.root' % jobid)
        ]
#        M0_file = root_open('../results/%s/gen_partons/ttJetsM0.root' % jobid, 'read')
#        M0_normfile = views.NormalizeView(root_open('../results/%s/gen_partons/ttJetsM0.root' % jobid, 'read'))
#        M700_file = root_open('../results/%s/gen_partons/ttJetsM700.root' % jobid, 'read')
#        M700_normfile = views.NormalizeView(root_open('../results/%s/gen_partons/ttJetsM700.root' % jobid, 'read'))
#        M1000_file = root_open('../results/%s/gen_partons/ttJetsM1000.root' % jobid, 'read')
#        M1000_normfile = views.NormalizeView(root_open('../results/%s/gen_partons/ttJetsM1000.root' % jobid, 'read'))

    else:
        myfile = root_open('../results/%s/gen_partons/%s.root' % (jobid, args.sample), 'read')
        normfile = views.NormalizeView(root_open('../results/%s/gen_partons/%s.root' % (jobid, args.sample), 'read'))

    plotter = BasePlotter(
    	'plots/gen_partons/%s/%s/%s/%s' % (jobid, args.analysis, args.sample, args.plot),
    	defaults = {'save' : {'png' : True, 'pdf' : False}}
    )

def stack_plots(lists):
    lists.sort()
    total = 0
    ratio_hists = []
    Stack_hists = []
 
    for i in lists:
        total += i
    for i in lists:
        ratio_hists.append(i/total)
        i.SetFillStyle(1001)
        Stack_hists.append(i)

    if len(Stack_hists) == 4:
        stack = plotter.create_stack(Stack_hists[0], Stack_hists[1], Stack_hists[2], Stack_hists[3])
        norm_stack = plotter.create_stack(Stack_hists[0]/total, Stack_hists[1]/total, Stack_hists[2]/total, Stack_hists[3]/total)
    elif len(Stack_hists) == 3:
        stack = plotter.create_stack(Stack_hists[0], Stack_hists[1], Stack_hists[2])
        norm_stack = plotter.create_stack(Stack_hists[0]/total, Stack_hists[1]/total, Stack_hists[2]/total)

    return stack, norm_stack, ratio_hists




##### Global Var. Definitions ####
defcol = 'black' # default color
defyax = 'a.u.' # yaxis title
if args.sample == "ttJetsM0" or args.sample == "Combined":
    mass_min = 250.
    mass_max = 2000.
if args.sample == "ttJetsM700":
	mass_min = 700.
	mass_max = 1000.
if args.sample == "ttJetsM1000":
	mass_min = 1000.
	mass_max = 2000.

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
		('DRmin_thad_vs_ptthad', 't_{h} p_{T}', 'min t_{h} #Delta R'),
		('DRmin_tlep_vs_pttlep', 't_{l} p_{T}', 'min t_{l} #Delta R'),
		('DRmax_thad_vs_ptthad', 't_{h} p_{T}', 'max t_{h} #Delta R')
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
if args.plot == "Kin_Vars":
	# Gen Object Plots
    Gen_Objs = { 'Lep' : 'l', 'BHad' : 'b_{h}', 'BLep' : 'b_{l}', 'WJa' : 'WJa', 'WJb' : 'WJb', 'THad' : 't_{h}', 'TLep' : 't_{l}' }

    for obj in Gen_Objs:
        label = Gen_Objs[obj]

#        kinvar = {'Mass' : ['m_{%s} [GeV]' % label, 'A.U.', 'Mass_Plots/']}
        kinvar = {'Mass' : ['m_{%s} [GeV]' % label, 'A.U.', 'Mass_Plots/'], 'Costh' : [' Cos #theta', 'A.U.', 'Costh_Plots/'], 'Eta' : [' #eta', 'A.U.', 'Eta_Plots/'], 'Pt' : [' p_{T} [GeV]', 'A.U.', 'Pt_Plots/']}

        for kvars in kinvar:
            plotter.set_subdir(kvars)
            xaxis = kinvar[kvars][0]
            yaxis = kinvar[kvars][1]
            kin_dir = kinvar[kvars][2]    

            var = kvars+'_'+obj

            hist = asrootpy(myfile.Get(kin_dir+var)).Clone()
            mean = hist.GetMean()
            rms = hist.GetRMS()

            plotter.set_histo_style(hist, color=defcol)
            if kvars == 'Mass':
                plotter.plot(hist, xtitle=xaxis, ytitle=defyax, drawstyle='hist')
            else:
                plotter.plot(hist, xtitle=label+xaxis, ytitle=defyax, drawstyle='hist')
            box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
            box.Draw()
            plotter.save(var)

            if kvars == 'Mass':
                continue
            
            var_2D = var+'_vs_Mtt'
            #print var_2D

            hist2D = asrootpy(myfile.Get(kin_dir+var_2D)).Clone()
            plotter.plot(hist2D)
            hist2D.Draw('colz')
            hist2D.xaxis.set_title('m_{t#bar t} [GeV]')
            hist2D.yaxis.set_title(label+xaxis)
            plotter.save(var_2D)

    
    print '\ncp -r '+os.getcwd()+'/plots/gen_partons/%s/%s/%s/%s . \n' % (jobid, args.analysis, args.sample, args.plot)

	
########################################################################################
if args.plot == "System":
	# Gen Object Plots
    kinvar = {'Mass_ttbar' : ['m_{t#bar t} [GeV]', 'A.U.', 'Mass_Plots/'], 'Costh_ttbar' : ['Cos #theta_{t#bar t}', 'A.U.', 'Costh_Plots/'],\
            'Eta_ttbar' : ['#eta_{t#bar t}', 'A.U.', 'Eta_Plots/'], 'Pt_ttbar' : [' p_{T_{t#bar t}} [GeV]', 'A.U.', 'Pt_Plots/'], 'nJets' : ['n_{jets}', 'A.U.', 'Sys_Plots/']}

    for kvar in kinvar:
        xaxis = kinvar[kvar][0]
        yaxis = kinvar[kvar][1]
        kin_dir = kinvar[kvar][2]

        hist = asrootpy(myfile.Get(kin_dir+kvar)).Clone()
        mean = hist.GetMean()
        rms = hist.GetRMS()

        plotter.set_histo_style(hist, color=defcol)
        plotter.plot(hist, xtitle=xaxis, ytitle=yaxis, drawstyle='hist')
        box = plotter.make_text_box('Mean = %f\nRMS = %f' % (mean,rms), position='NE')
        box.Draw()
        plotter.save(kvar)

#		efficiencies = []
#		total= hist.Integral()
#		for i in range(0, hist.GetXaxis().GetNbins()):
#			efficiencies.append(format(hist.GetBinContent(i+1)/total, '.4f'))
#		print(var, efficiencies)

    print '\ncp -r '+os.getcwd()+'/plots/gen_partons/%s/%s/%s/%s . \n' % (jobid, args.analysis, args.sample, args.plot)

##################################################################################################
if args.plot == "Merged_Fractions":
    other_pairs = { 'LepBLep' : ['(l, b_{l})', '#cc0066'], 'LepBHad' : ['(l, b_{h})', 'blue'], 'LepWJa' : ['(l, WJa)', '#ff9900'], 'LepWJb' : ['(l, WJb)', '#00ffff'],\
                'BHadBLep' : ['(b_{h}, b_{l})', '#ff66ff'], 'BLepWJa' : ['(b_{l}, WJa)', '#00ff00'], 'BLepWJb' : ['(b_{l}, WJb)', 'black']}
#                'BHadWJa' : ['(b_{h}, WJa)', 'red'], 'BHadWJb' : ['(b_{h}, WJb)', 'blue'], 'WJaWJb' : ['(WJa, WJb)', 'green']}

    had_pairs = {'BHadWJa' : ['(b_{h}, WJa)', 'red'], 'BHadWJb' : ['(b_{h}, WJb)', 'blue'], 'WJaWJb' : ['(WJa, WJb)', 'green']}

    xaxis = 'M_{t#bar t} [GeV]'
    yaxis = '#Delta R < 0.4 Fraction'

    hp_hist = []
    for objs in had_pairs:
        legends = had_pairs[objs][0]
        cols = had_pairs[objs][1]

        LP4 = asrootpy(myfile.Get('Gen_Plots/DR_'+objs+'_lDRP4_vs_mttbar')).Clone()
        GP4 = asrootpy(myfile.Get('Gen_Plots/DR_'+objs+'_gDRP4_vs_mttbar')).Clone()
        if LP4.Integral() == 0:
            continue
        Ratio = LP4/(LP4+GP4)
        plotter.set_histo_style(Ratio, title=legends, color=cols, markerstyle=20, linestyle='solid')
        hp_hist.append(Ratio)

    plotter.overlay(hp_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xaxis, ytitle=yaxis)
    #plotter.overlay(hp_hist, legend_def=LegendDefinition(position='NW'), logy=True, y_range=(10**(-2),10**0), legendstyle='l', xtitle=xaxis, ytitle=yaxis)
    plotter.save('Merged_Fractions_Had_Partons_DRP4')

    op_hist = []
    for objs in other_pairs:
        legends = other_pairs[objs][0]
        cols = other_pairs[objs][1]

        LP4 = asrootpy(myfile.Get('Gen_Plots/DR_'+objs+'_lDRP4_vs_mttbar')).Clone()
        GP4 = asrootpy(myfile.Get('Gen_Plots/DR_'+objs+'_gDRP4_vs_mttbar')).Clone()
        if LP4.Integral() == 0:
            continue
        Ratio = LP4/(LP4+GP4)
        plotter.set_histo_style(Ratio, title=legends, color=cols, markerstyle=20, linestyle='solid')
        op_hist.append(Ratio)

    plotter.overlay(op_hist, legend_def=LegendDefinition(position='NW'), logy=True, y_range=(10**(-5),10**2), legendstyle='l', xtitle=xaxis, ytitle=yaxis)
    plotter.save('Merged_Fractions_Other_Partons_DRP4')

#    DRmin_ratios = [
#    #		('DRmin_thad_lDRP4_vs_mttbar', 'DRmin_thad_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', 't_{h} min #Delta R', 'DRmin_thad_lp4_ratio_vs_mttbar'),
#    #		('DRmin_thad_lDRP4_vs_ptthad', 'DRmin_thad_gDRP4_vs_ptthad', 't_{h} p_{T}', '#Delta R < 0.4 Fraction', 't_{h} min #Delta R', 'DRmin_thad_lp4_ratio_vs_ptthad'),
#    #		('DRmin_tlep_lDRP4_vs_mttbar', 'DRmin_tlep_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', 't_{l} min #Delta R', 'DRmin_tlep_lp4_ratio_vs_mttbar'),
#    #		('DRmin_tlep_lDRP4_vs_pttlep', 'DRmin_tlep_gDRP4_vs_pttlep', 't_{l} p_{T}', '#Delta R < 0.4 Fraction', 't_{l} min #Delta R', 'DRmin_tlep_lp4_ratio_vs_pttlep'),
#    	
#    	('DR_LepBHad_lDRP4_vs_mttbar', 'DR_LepBHad_gDRP4_vs_mttbar', '#Delta R (l, b_{h})', 'DR_LepBHad_lp4_ratio_vs_mttbar'),
#        ('DR_LepBLep_lDRP4_vs_mttbar', 'DR_LepBLep_gDRP4_vs_mttbar', '#Delta R (l, b_{l})', 'DR_LepBLep_lp4_ratio_vs_mttbar'),
#        ('DR_LepWJa_lDRP4_vs_mttbar', 'DR_LepWJa_gDRP4_vs_mttbar', '#Delta R (l, WJa)', 'DR_LepWJa_lp4_ratio_vs_mttbar'),
#        ('DR_LepWJb_lDRP4_vs_mttbar', 'DR_LepWJb_gDRP4_vs_mttbar', '#Delta R (l, WJb)', 'DR_LepWJb_lp4_ratio_vs_mttbar'),
#        ('DR_BHadBLep_lDRP4_vs_mttbar', 'DR_BHadBLep_gDRP4_vs_mttbar', '#Delta R (b_{h}, b_{l})', 'DR_BHadBLep_lp4_ratio_vs_mttbar'),
#        #('DR_BHadWJa_lDRP4_vs_mttbar', 'DR_BHadWJa_gDRP4_vs_mttbar', '#Delta R (b_{h}, WJa)', 'DR_BHadWJa_lp4_ratio_vs_mttbar'),
#        #('DR_BHadWJb_lDRP4_vs_mttbar', 'DR_BHadWJb_gDRP4_vs_mttbar', '#Delta R (b_{h}, WJb)', 'DR_BHadWJb_lp4_ratio_vs_mttbar'),
#        ('DR_BLepWJa_lDRP4_vs_mttbar', 'DR_BLepWJa_gDRP4_vs_mttbar', '#Delta R (b_{l}, WJa)', 'DR_BLepWJa_lp4_ratio_vs_mttbar'),
#        ('DR_BLepWJb_lDRP4_vs_mttbar', 'DR_BLepWJb_gDRP4_vs_mttbar', '#Delta R (b_{l}, WJb)', 'DR_BLepWJb_lp4_ratio_vs_mttbar'),
#        #('DR_WJaWJb_lDRP4_vs_mttbar', 'DR_WJaWJb_gDRP4_vs_mttbar', '#Delta R (WJa, WJb)', 'DR_WJaWJb_lp4_ratio_vs_mttbar'),
#    #		('DRmax_thad_lDRP4_vs_mttbar', 'DRmax_thad_gDRP4_vs_mttbar', 'm_{t#bar t}', '#Delta R < 0.4 Fraction', 't_{h} max #Delta R', 'DRmax_thad_lp4_ratio_vs_mttbar'),
#    #		('DRmax_thad_lDRP4_vs_ptthad', 'DRmax_thad_gDRP4_vs_ptthad', 't_{h} p_{T}', '#Delta R < 0.4 Fraction', 't_{h} max #Delta R', 'DRmax_thad_lp4_ratio_vs_ptthad')
#    ]
#
#    Had_Side_Ratios = [
#        ('DR_BHadWJa_lDRP4_vs_mttbar', 'DR_BHadWJa_gDRP4_vs_mttbar', '#Delta R (b_{h}, WJa)', 'DR_BHadWJa_lp4_ratio_vs_mttbar'),
#        ('DR_BHadWJb_lDRP4_vs_mttbar', 'DR_BHadWJb_gDRP4_vs_mttbar', '#Delta R (b_{h}, WJb)', 'DR_BHadWJb_lp4_ratio_vs_mttbar'),
#        ('DR_WJaWJb_lDRP4_vs_mttbar', 'DR_WJaWJb_gDRP4_vs_mttbar', '#Delta R (WJa, WJb)', 'DR_WJaWJb_lp4_ratio_vs_mttbar'),
#    ]
#
#    for lp4, gp4, title, name in DRmin_ratios:# < 0.4, > 0.4, xaxis, yaxis, png name
#    	efficiencies = []
#    	
#    	LP4 = asrootpy(myfile.Get('Gen_Plots/'+lp4)).Clone()
#    	if LP4.Integral() != 0:
#    		GP4 = asrootpy(myfile.Get('Gen_Plots/'+gp4)).Clone()
#    		DR_Ratio = LP4/(LP4+GP4)
#    		plotter.set_histo_style(DR_Ratio, color=defcol, markerstyle=20)
#    		plotter.plot(DR_Ratio)
#    		DR_Ratio.Draw("P")
#    		DR_Ratio.xaxis.set_title(xaxis)
#    		DR_Ratio.yaxis.set_title(yaxis)
#    	#		DR_Ratio.yaxis.range_user = 0, 1
#    		box = plotter.make_text_box(title, position='NE')
#    		plotter.save(name)
#    	
#    	for i in range(0, DR_Ratio.GetXaxis().GetNbins()):
#    		efficiencies.append(format(DR_Ratio.GetBinContent(i+1), '.4f'))
#    	#	mass_range.append(Mu_Eff.GetXaxis().GetBinLowEdge(i+1))
#    	print(yaxis, efficiencies)
    
    print '\ncp -r '+os.getcwd()+'/plots/gen_partons/%s/%s/%s/%s . \n' % (jobid, args.analysis, args.sample, args.plot)


##################################################################################################
if args.plot == "Had_Comp_DR":
### compares relative frequencies of:
### 1.   three jets on hadronic side being resolved
### 2.   only bhad and one w jet merging
### 3.   only jets from w merging
### 4.   all three hadronic jets merging
### 5.   non-reconstructable events


    ### Events with all jets at different DR values
        Hadronic_Events_M = [
                ('Gen_Had_Resolved_vs_mttbar', 'Resolved', 'blue'),
                ('Gen_Merged_BHadWJet_vs_mttbar', 'Merged b_{h} and W jet', 'red'),
                ('Gen_Merged_WJets_vs_mttbar', 'Merged W jets', 'black'),
                ('Gen_Merged_THad_vs_mttbar', 'All 3 Merged', 'green')
#                ('Gen_Non_Reconstructable_vs_mttbar', 'Non Reconstructable', 'magenta')
        ]

        for name, title, col in DRvals:

                to_draw = []

                for var, legends, colors in Hadronic_Events_M:
                        hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
                        plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
                        to_draw.append(hist)


                ### plot log of event occurrences
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='m_{t#bar t} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
                plotter.save('AllJets_Hadronic_Events_vs_mttbar_log_'+name)


                stack, norm_stack, Ratio_Hists = stack_plots(to_draw)

                ### plot event ratios
                plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.15), xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
                plotter.save('AllJets_Hadronic_Events_vs_mttbar_ratio_'+name)


                ### create stacked hists of ratios
                plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title)
                plotter.save('AllJets_Hadronic_Events_vs_mttbar_stack_'+name)

                plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{t#bar t} [GeV]', ytitle='Event Fraction '+title)
                plotter.save('AllJets_Hadronic_Events_vs_mttbar_stack_Norm_'+name)


        Hadronic_Events_PT = [
                ('Gen_Had_Resolved_vs_thadpt', 'Resolved', 'blue'),
                ('Gen_Merged_BHadWJet_vs_thadpt', 'Merged b_{h} and W jet', 'red'),
                ('Gen_Merged_WJets_vs_thadpt', 'Merged W jets', 'black'),
                ('Gen_Merged_THad_vs_thadpt', 'All 3 Merged', 'green')
#                ('Gen_Non_Reconstructable_vs_thadpt', 'Non Reconstructable', 'magenta')
        ]

        for name, title, col in DRvals:

                to_draw = []

                for var, legends, colors in Hadronic_Events_PT:
                        hist = asrootpy(myfile.Get(name+'/'+var)).Clone()
                        plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
                        to_draw.append(hist)


                ### plot log of event occurrences
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='t_{h} p_{T} [GeV]', ytitle='A.U. '+title, drawstyle='hist')
                plotter.save('AllJets_Hadronic_Events_vs_thadpt_log_'+name)


                stack, norm_stack, Ratio_Hists = stack_plots(to_draw)


                ### plot event ratios
                plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.15), xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction '+title, drawstyle='hist')
                plotter.save('AllJets_Hadronic_Events_vs_thadpt_ratio_'+name)


                ### create stacked hists of ratios
                plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction '+title)
                plotter.save('AllJets_Hadronic_Events_vs_thadpt_stack_'+name)

                plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction '+title)
                plotter.save('AllJets_Hadronic_Events_vs_thadpt_stack_Norm_'+name)

#################################################################
if args.plot == "AllJets" or args.plot == "2LJ" or args.plot == "3J" or args.plot == "4J" or args.plot == "5PJ":

#    kinvar = {'THadpt_vs_Mttbar' : ['M_{t #bar{t}} [GeV]', 'Leading top p_{T} [GeV]']}
    kinvar = {'Mtt' : ['M_{t #bar{t}} [GeV]', 'A.U.'], 'THadpt' : ['t_{h} p_{T} [GeV]', 'A.U.'], 'THadpt_vs_Mttbar' : ['M_{t #bar{t}} [GeV]', 't_{h} p_{T} [GeV]']}
    event_type = {'Had_Resolved_DRP4' : ['Resolved #Delta R < 0.4', 'blue'], 'Merged_THad_DRP8' : ['All 3 Merged #Delta R < 0.8', 'green'], 'Partially_Merged_DRP8' : ['Partial Merge #Delta R < 0.8', 'black']}
    jets = {'AllJets' : '', '2LJ' : '_2LJ', '3J': '_3J', '4J': '_4J', '5PJ' : '_5PJ'}

    xbins = np.linspace(mass_min, mass_max, 21)
    low_pt = np.linspace(0, 500, 26)
    med_pt = np.linspace(550, 1000, 10)
    high_pt = [1200, 1400, 1700, 2000]
    ybins = np.concatenate((low_pt, med_pt, high_pt), axis=0)
    
    zmin, zmax = 0., 1.
    cmap_levels = np.linspace(zmin, zmax, 100)
    #log_cmap_levels = np.logspace(-1., 7., num=100)
    
    if args.sample == "Combined":
        for kin in kinvar:
            plotter.set_subdir(kin)
      
            xaxis = kinvar[kin][0]
            yaxis = kinvar[kin][1]  
        
            to_draw = []
            Pt_vs_M_Hists = []
        
            for evt in event_type:
                legends = event_type[evt][0]
                col = event_type[evt][1]
    
                var = 'Gen_%s_%s' % (evt, kin)

                comb_hist = 0
                for fname in ttJets_files:
                    for name in fname:
                        tfile = io.root_open(name)
                        hist = asrootpy(tfile.Get('Had_comp/'+var+jets[args.plot])).Clone()
                        comb_hist += hist
                if kin == 'THadpt_vs_Mttbar':
                    if 'Resolved' in var:
                        Res = comb_hist
                        Pt_vs_M_Hists.append(comb_hist)
                    elif 'Partially' in var:
                        Part = comb_hist
                        Pt_vs_M_Hists.append(comb_hist)
                    else:
                        Full = comb_hist
                        Pt_vs_M_Hists.append(comb_hist)

                    if len(Pt_vs_M_Hists) == 3:   
                        #rebin hists
                        Res = RebinView.newRebin2D(Res,xbins,ybins)
                        Full = RebinView.newRebin2D(Full,xbins,ybins)
                        Part = RebinView.newRebin2D(Part,xbins,ybins)
                        
                        tot = Res+Full+Part #combine all merging categories together
                        #set_trace()
                        plotter.plot(tot, drawstyle='hist')
                        tot.Draw('colz')
                        tot.xaxis.set_title(xaxis)
                        tot.yaxis.set_title(yaxis)
                        plotter.pad.SetLogz(True)
                        #tot.SetContour(len(log_cmap_levels), log_cmap_levels)
                        plotter.save('%s_Combined_Merge_Categories' % args.plot)

                        Res_frac = Res/tot # resolved plots
                        plotter.plot(Res_frac, drawstyle='hist')
                        Res_frac.Draw('colz')
                        Res_frac.zaxis.range_user = zmin, zmax
                        Res_frac.xaxis.set_title(xaxis)
                        Res_frac.yaxis.set_title(yaxis)
                        Res_frac.SetContour(len(cmap_levels), cmap_levels)
                        plotter.save('%s_Had_Resolved_DRP4_Ratio' % args.plot)
                        
                        Full_frac = Full/tot #fully merged plots
                        plotter.plot(Full_frac, drawstyle='hist')
                        Full_frac.Draw('colz')
                        Full_frac.zaxis.range_user = zmin, zmax
                        Full_frac.xaxis.set_title(xaxis)
                        Full_frac.yaxis.set_title(yaxis)
                        Full_frac.SetContour(len(cmap_levels), cmap_levels)
                        plotter.save('%s_Merged_THad_DRP8_Ratio' % args.plot)
                        
                        Part_frac = Part/tot #partially merged plots
                        plotter.plot(Part_frac, drawstyle='hist')
                        Part_frac.Draw('colz')
                        Part_frac.zaxis.range_user = zmin, zmax
                        Part_frac.xaxis.set_title(xaxis)
                        Part_frac.yaxis.set_title(yaxis)
                        Part_frac.SetContour(len(cmap_levels), cmap_levels)
                        plotter.save('%s_Partially_Merged_DRP8_Ratio' % args.plot)

                        Tot_Ratio = tot # total ratio plots
                        tot_xproj = tot.ProjectionX("",1,tot.GetNbinsY()) # get projection of total hist along x
                        for binsx in range(len(xbins)-1):
                            for binsy in range(len(ybins)-1):
                                #for each mtt bin (x), divide the number of events per pt bin (y)
                                #by the total number of events in the entire mtt bin
                                #the sum in each x bin is 1
                                Tot_Ratio.SetBinContent(binsx+1, binsy+1,tot.GetBinContent(binsx+1,binsy+1)/tot_xproj.GetBinContent(binsx+1))

                        plotter.plot(Tot_Ratio, drawstyle='hist')
                        Tot_Ratio.Draw('colz')
                        #Tot_Ratio.zaxis.range_user = zmin, zmax
                        Tot_Ratio.xaxis.set_title(xaxis)
                        Tot_Ratio.yaxis.set_title(yaxis)
                        plotter.save('%s_Combined_Merge_Categories_Ratio' % args.plot)
                        plotter.plot(Tot_Ratio, drawstyle='hist')
                        Tot_Ratio.Draw('colz')
                        Tot_Ratio.xaxis.set_title(xaxis)
                        Tot_Ratio.yaxis.set_title(yaxis)
                        plotter.pad.SetLogz(True)
                        #Tot_Ratio.SetContour(len(cmap_levels), cmap_levels)
                        plotter.save('%s_Combined_Merge_Categories_Ratio_Log' % args.plot)

                else:
                    plotter.set_histo_style(comb_hist, title=legends, drawstyle='hist', color=col, linestyle='solid')
                    to_draw.append(comb_hist)
    
                    if len(to_draw) == 3:
                        stack, norm_stack, Ratio_Hists = stack_plots(to_draw)
                        
                        ### plot event ratios
                        plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle=xaxis, ytitle='Event Fraction', drawstyle='hist')
                        plotter.save('%s_Hadronic_Events_%s_Comp_ratio' % (args.plot, kin))
                
                        ### create stacked hists of ratios
                        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xaxis, ytitle='Event Fraction')
                        plotter.save('%s_Hadronic_Events_%s_Comp_stack' % (args.plot, kin))
                
                        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xaxis, ytitle='Event Fraction')
                        plotter.save('%s_Hadronic_Events_%s_Comp_stack_Norm' % (args.plot, kin))
            
                        for i in to_draw:
                            i.SetFillStyle(0)
                        ### plot log of event occurrences
                        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle=xaxis, ytitle=yaxis, drawstyle='hist')
                        plotter.save('%s_Hadronic_Events_%s_Comp_log' % (args.plot, kin))


    else:
        for kin in kinvar:
            plotter.set_subdir(kin)
      
            xaxis = kinvar[kin][0]
            yaxis = kinvar[kin][1]  
        
            to_draw = []
        
            if kin == 'THadpt_vs_Mttbar':
                Res = asrootpy(myfile.Get(str('Had_comp/Gen_Had_Resolved_DRP4_%s'+jets[args.plot]) % kin)).Clone()
                Full = asrootpy(myfile.Get(str('Had_comp/Gen_Merged_THad_DRP8_%s'+jets[args.plot]) % kin)).Clone()
                Part = asrootpy(myfile.Get(str('Had_comp/Gen_Partially_Merged_DRP8_%s'+jets[args.plot]) % kin)).Clone()
    
                #rebin hists
                Res = RebinView.newRebin2D(Res,xbins,ybins)
                Full = RebinView.newRebin2D(Full,xbins,ybins)
                Part = RebinView.newRebin2D(Part,xbins,ybins)
    
                tot = Res+Full+Part #combine all merging categories together
                #set_trace()
                plotter.plot(tot, drawstyle='hist')
                tot.Draw('colz')
                tot.xaxis.set_title(xaxis)
                tot.yaxis.set_title(yaxis)
                plotter.pad.SetLogz(True)
                #tot.SetContour(len(log_cmap_levels), log_cmap_levels)
                plotter.save('%s_Combined_Merge_Categories' % args.plot)

                Res_frac = Res/tot #fully resolved plots
                plotter.plot(Res_frac, drawstyle='hist')
                Res_frac.Draw('colz')
                Res_frac.zaxis.range_user = zmin, zmax
                Res_frac.xaxis.set_title(xaxis)
                Res_frac.yaxis.set_title(yaxis)
                Res_frac.SetContour(len(cmap_levels), cmap_levels)
                plotter.save('%s_Had_Resolved_DRP4_Ratio' % args.plot)
    
                Full_frac = Full/tot #fully merged plots
                plotter.plot(Full_frac, drawstyle='hist')
                Full_frac.Draw('colz')
                Full_frac.zaxis.range_user = zmin, zmax
                Full_frac.xaxis.set_title(xaxis)
                Full_frac.yaxis.set_title(yaxis)
                Full_frac.SetContour(len(cmap_levels), cmap_levels)
                plotter.save('%s_Merged_THad_DRP8_Ratio' % args.plot)
    
                Part_frac = Part/tot #partially merged plots
                plotter.plot(Part_frac, drawstyle='hist')
                Part_frac.Draw('colz')
                Part_frac.zaxis.range_user = zmin, zmax
                Part_frac.xaxis.set_title(xaxis)
                Part_frac.yaxis.set_title(yaxis)
                Part_frac.SetContour(len(cmap_levels), cmap_levels)
                plotter.save('%s_Partially_Merged_DRP8_Ratio' % args.plot)
    
                Tot_Ratio = tot #total ratio plots
                tot_xproj = tot.ProjectionX("",1,tot.GetNbinsY()) # get projection of total hist along x
                for binsx in range(len(xbins)-1):
                    for binsy in range(len(ybins)-1):
                        #for each mtt bin (x), divide the number of events per pt bin (y)
                        #by the total number of events in the entire mtt bin
                        #the sum in each x bin is 1
                        Tot_Ratio.SetBinContent(binsx+1, binsy+1,tot.GetBinContent(binsx+1,binsy+1)/tot_xproj.GetBinContent(binsx+1))

                plotter.plot(Tot_Ratio, drawstyle='hist')
                Tot_Ratio.Draw('colz')
                #Tot_Ratio.zaxis.range_user = zmin, zmax
                Tot_Ratio.xaxis.set_title(xaxis)
                Tot_Ratio.yaxis.set_title(yaxis)
                plotter.save('%s_Combined_Merge_Categories_Ratio' % args.plot)
                plotter.plot(Tot_Ratio, drawstyle='hist')
                Tot_Ratio.Draw('colz')
                Tot_Ratio.xaxis.set_title(xaxis)
                Tot_Ratio.yaxis.set_title(yaxis)
                plotter.pad.SetLogz(True)
                #Tot_Ratio.SetContour(len(cmap_levels), cmap_levels)
                plotter.save('%s_Combined_Merge_Categories_Ratio_Log' % args.plot)

            else:    
                for evt in event_type:
                    legends = event_type[evt][0]
                    col = event_type[evt][1]
        
                    var = 'Gen_%s_%s' % (evt, kin)
                    hist = asrootpy(myfile.Get('Had_comp/'+var+jets[args.plot])).Clone()
                    plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=col, linestyle='solid')
                    to_draw.append(hist)
   
                stack, norm_stack, Ratio_Hists = stack_plots(to_draw)
    
                ### plot event ratios
                plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle=xaxis, ytitle='Event Fraction', drawstyle='hist')
                plotter.save('%s_Hadronic_Events_%s_Comp_ratio' % (args.plot, kin))
        
                ### create stacked hists of ratios
                plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xaxis, ytitle='Event Fraction')
                plotter.save('%s_Hadronic_Events_%s_Comp_stack' % (args.plot, kin))
        
                plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xaxis, ytitle='Event Fraction')
                plotter.save('%s_Hadronic_Events_%s_Comp_stack_Norm' % (args.plot, kin))
    
                for i in to_draw:
                    i.SetFillStyle(0)
                ### plot log of event occurrences
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle=xaxis, ytitle=yaxis, drawstyle='hist')
                plotter.save('%s_Hadronic_Events_%s_Comp_log' % (args.plot, kin))

    print '\ncp -r '+os.getcwd()+'/plots/gen_partons/%s/%s/%s/%s . \n' % (jobid, args.analysis, args.sample, args.plot)

#######################################################################

if args.plot == "3Partons":

#### one parton lost from acceptance

        Hadronic_Events_3Part = [
                ('Three_Parton_Gen_Had_Resolved_DRP4_vs_thadpt', 'Resolved #Delta R < 0.4', 'blue'),
                ('Three_Parton_Gen_Merged_BHadWJet_DRP8_vs_thadpt', 'Merged b_{h} and W jet #Delta R < 0.8', 'red'),
                ('Three_Parton_Gen_Merged_WJets_DRP8_vs_thadpt', 'Merged W jets #Delta R < 0.8', 'black'),
                ('Three_Parton_Gen_Merged_THad_DRP8_vs_thadpt', 'All 3 Merged #Delta R < 0.8', 'green')
#                ('Gen_Non_Reconstructable_vs_thadpt', 'Non Reconstructable', 'magenta')
        ]

        to_draw = []

        for var, legends, colors in Hadronic_Events_3Part:
                hist = asrootpy(myfile.Get('Had_comp/'+var)).Clone()
                plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
                to_draw.append(hist)

        ### plot log of event occurrences
        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='t_{h} p_{T} [GeV]', ytitle='A.U.', drawstyle='hist')
        plotter.save('Three_Partons_Hadronic_Events_Comp_vs_thadpt_log')


        stack, norm_stack, Ratio_Hists = stack_plots(to_draw)

        ### plot event ratios
        plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction', drawstyle='hist')
        plotter.save('Three_Partons_Hadronic_Events_Comp_vs_thadpt_ratio')

        ### create stacked hists of ratios
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
        plotter.save('Three_Partons_Hadronic_Events_Comp_vs_thadpt_stack')

        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
        plotter.save('Three_Partons_Hadronic_Events_Comp_vs_thadpt_stack_Norm')


    ### bhad missing
        Hadronic_Events_3Part_BHad_Missing = [
                ('Three_Parton_BHad_Missing_Gen_Had_Resolved_DRP4_vs_thadpt', 'Resolved #Delta R < 0.4', 'blue'),
                ('Three_Parton_BHad_Missing_Gen_Merged_BHadWJet_DRP8_vs_thadpt', 'Merged b_{h} and W jet #Delta R < 0.8', 'red'),
                ('Three_Parton_BHad_Missing_Gen_Merged_WJets_DRP8_vs_thadpt', 'Merged W jets #Delta R < 0.8', 'black'),
                ('Three_Parton_BHad_Missing_Gen_Merged_THad_DRP8_vs_thadpt', 'All 3 Merged #Delta R < 0.8', 'green')
#                ('Gen_Non_Reconstructable_vs_thadpt', 'Non Reconstructable', 'magenta')
        ]

        to_draw = []

        for var, legends, colors in Hadronic_Events_3Part_BHad_Missing:
                hist = asrootpy(myfile.Get('Had_comp/'+var)).Clone()
                plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
                to_draw.append(hist)

        ### plot log of event occurrences
        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='t_{h} p_{T} [GeV]', ytitle='A.U.', drawstyle='hist')
        box = plotter.make_text_box('b_{h} Lost', position='NW')
        plotter.save('Three_Partons_BHad_Missing_Hadronic_Events_Comp_vs_thadpt_log')


        stack, norm_stack, Ratio_Hists = stack_plots(to_draw)

        ### plot event ratios
        plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction', drawstyle='hist')
        box = plotter.make_text_box('b_{h} Lost', position='NW')
        plotter.save('Three_Partons_BHad_Missing_Hadronic_Events_Comp_vs_thadpt_ratio')

        ### create stacked hists of ratios
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
        box = plotter.make_text_box('b_{h} Lost', position='NW')
        plotter.save('Three_Partons_BHad_Missing_Hadronic_Events_Comp_vs_thadpt_stack')

        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
        box = plotter.make_text_box('b_{h} Lost', position='NW')
        plotter.save('Three_Partons_BHad_Missing_Hadronic_Events_Comp_vs_thadpt_stack_Norm')


    ### blep missing
        Hadronic_Events_3Part_BLep_Missing = [
                ('Three_Parton_BLep_Missing_Gen_Had_Resolved_DRP4_vs_thadpt', 'Resolved #Delta R < 0.4', 'blue'),
                ('Three_Parton_BLep_Missing_Gen_Merged_BHadWJet_DRP8_vs_thadpt', 'Merged b_{h} and W jet #Delta R < 0.8', 'red'),
                ('Three_Parton_BLep_Missing_Gen_Merged_WJets_DRP8_vs_thadpt', 'Merged W jets #Delta R < 0.8', 'black'),
                ('Three_Parton_BLep_Missing_Gen_Merged_THad_DRP8_vs_thadpt', 'All 3 Merged #Delta R < 0.8', 'green')
        ]

        to_draw = []

        for var, legends, colors in Hadronic_Events_3Part_BLep_Missing:
                hist = asrootpy(myfile.Get('Had_comp/'+var)).Clone()
                plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
                to_draw.append(hist)

        ### plot log of event occurrences
        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='t_{h} p_{T} [GeV]', ytitle='A.U.', drawstyle='hist')
        box = plotter.make_text_box('b_{l} Lost', position='NW')
        plotter.save('Three_Partons_BLep_Missing_Hadronic_Events_Comp_vs_thadpt_log')


        stack, norm_stack, Ratio_Hists = stack_plots(to_draw)

        ### plot event ratios
        plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction', drawstyle='hist')
        box = plotter.make_text_box('b_{l} Lost', position='NW')
        plotter.save('Three_Partons_BLep_Missing_Hadronic_Events_Comp_vs_thadpt_ratio')

        ### create stacked hists of ratios
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
        box = plotter.make_text_box('b_{l} Lost', position='NW')
        plotter.save('Three_Partons_BLep_Missing_Hadronic_Events_Comp_vs_thadpt_stack')

        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
        box = plotter.make_text_box('b_{l} Lost', position='NW')
        plotter.save('Three_Partons_BLep_Missing_Hadronic_Events_Comp_vs_thadpt_stack_Norm')



    ### wja missing
        Hadronic_Events_3Part_WJa_Missing = [
                ('Three_Parton_WJa_Missing_Gen_Had_Resolved_DRP4_vs_thadpt', 'Resolved #Delta R < 0.4', 'blue'),
                ('Three_Parton_WJa_Missing_Gen_Merged_BHadWJet_DRP8_vs_thadpt', 'Merged b_{h} and W jet #Delta R < 0.8', 'red'),
                ('Three_Parton_WJa_Missing_Gen_Merged_WJets_DRP8_vs_thadpt', 'Merged W jets #Delta R < 0.8', 'black'),
                ('Three_Parton_WJa_Missing_Gen_Merged_THad_DRP8_vs_thadpt', 'All 3 Merged #Delta R < 0.8', 'green')
        ]

        to_draw = []

        for var, legends, colors in Hadronic_Events_3Part_BLep_Missing:
                hist = asrootpy(myfile.Get('Had_comp/'+var)).Clone()
                plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
                to_draw.append(hist)

        ### plot log of event occurrences
        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='t_{h} p_{T} [GeV]', ytitle='A.U.', drawstyle='hist')
        box = plotter.make_text_box('WJa Lost', position='NW')
        plotter.save('Three_Partons_WJa_Missing_Hadronic_Events_Comp_vs_thadpt_log')


        stack, norm_stack, Ratio_Hists = stack_plots(to_draw)

        ### plot event ratios
        plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction', drawstyle='hist')
        box = plotter.make_text_box('WJa Lost', position='NW')
        plotter.save('Three_Partons_WJa_Missing_Hadronic_Events_Comp_vs_thadpt_ratio')

        ### create stacked hists of ratios
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
        box = plotter.make_text_box('WJa Lost', position='NW')
        plotter.save('Three_Partons_WJa_Missing_Hadronic_Events_Comp_vs_thadpt_stack')

        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
        box = plotter.make_text_box('WJa Lost', position='NW')
        plotter.save('Three_Partons_WJa_Missing_Hadronic_Events_Comp_vs_thadpt_stack_Norm')


    ### wjb missing
        Hadronic_Events_3Part_WJb_Missing = [
                ('Three_Parton_WJb_Missing_Gen_Had_Resolved_DRP4_vs_thadpt', 'Resolved #Delta R < 0.4', 'blue'),
                ('Three_Parton_WJb_Missing_Gen_Merged_BHadWJet_DRP8_vs_thadpt', 'Merged b_{h} and W jet #Delta R < 0.8', 'red'),
                ('Three_Parton_WJb_Missing_Gen_Merged_WJets_DRP8_vs_thadpt', 'Merged W jets #Delta R < 0.8', 'black'),
                ('Three_Parton_WJb_Missing_Gen_Merged_THad_DRP8_vs_thadpt', 'All 3 Merged #Delta R < 0.8', 'green')
        ]

        to_draw = []

        for var, legends, colors in Hadronic_Events_3Part_BLep_Missing:
                hist = asrootpy(myfile.Get('Had_comp/'+var)).Clone()
                plotter.set_histo_style(hist, title=legends, drawstyle='hist', color=colors, linestyle='solid')
                to_draw.append(hist)

        ### plot log of event occurrences
        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='t_{h} p_{T} [GeV]', ytitle='A.U.', drawstyle='hist')
        box = plotter.make_text_box('WJb Lost', position='NW')
        plotter.save('Three_Partons_WJb_Missing_Hadronic_Events_Comp_vs_thadpt_log')


        stack, norm_stack, Ratio_Hists = stack_plots(to_draw)

        ### plot event ratios
        plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction', drawstyle='hist')
        box = plotter.make_text_box('WJb Lost', position='NW')
        plotter.save('Three_Partons_WJb_Missing_Hadronic_Events_Comp_vs_thadpt_ratio')

        ### create stacked hists of ratios
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
        box = plotter.make_text_box('WJb Lost', position='NW')
        plotter.save('Three_Partons_WJb_Missing_Hadronic_Events_Comp_vs_thadpt_stack')

        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
        box = plotter.make_text_box('WJb Lost', position='NW')
        plotter.save('Three_Partons_WJb_Missing_Hadronic_Events_Comp_vs_thadpt_stack_Norm')

#####################################################################
#
#if args.plot == "5PJ":
#
##### events with 5+ jets

##### Jet event ratios
#
#        Hadronic_Resolved_Ratios = [
#                ('Gen_Had_Resolved_DRP4_vs_thadpt', 'Gen_Had_Resolved_DRP4_vs_thadpt_2LJ', 'Gen_Had_Resolved_DRP4_vs_thadpt_3J', 'Gen_Had_Resolved_DRP4_vs_thadpt_4J', 'Gen_Had_Resolved_DRP4_vs_thadpt_5PJ', 'Resolved #DeltaR < 0.4')
##                ('Gen_Merged_BHadWJet_DRP8_vs_thadpt', 'Gen_Merged_BHadWJet_DRP8_vs_thadpt_3J', 'Gen_Merged_BHadWJet_DRP8_vs_thadpt_4J', 'Gen_Merged_BHadWJet_DRP8_vs_thadpt_5PJ', 'Merged b_{h} and W jet #DeltaR < 0.8'),
##                 ('Gen_Merged_WJets_DRP8_vs_thadpt', 'Gen_Merged_WJets_DRP8_vs_thadpt_3J', 'Gen_Merged_WJets_DRP8_vs_thadpt_4J', 'Gen_Merged_WJets_DRP8_vs_thadpt_5PJ', 'Merged W jets #DeltaR < 0.8'),
##                 ('Gen_Merged_THad_DRP8_vs_thadpt', 'Gen_Merged_THad_DRP8_vs_thadpt_3J', 'Gen_Merged_THad_DRP8_vs_thadpt_4J', 'Gen_Merged_THad_DRP8_vs_thadpt_5PJ', 'All 3 Merged #DeltaR < 0.8')
##                ('Gen_Non_Reconstructable_vs_thadpt', 'Non Reconstructable', 'magenta')
#        ]
#
##        for name, title, col in DRvals:
#
#        to_draw = []
#        to_draw_all = []
#        Add_Ratio_Hists = []
#        Ratio_Hists = []
#
#        for All, Jet2L, Jet3, Jet4, Jet5P, legends in Hadronic_Resolved_Ratios:
#                total = asrootpy(myfile.Get('Had_comp/'+All)).Clone()
#                plotter.set_histo_style(total, title='All jets', drawstyle='hist', color='black', linestyle='solid')
#                to_draw_all.append(total)
#
#                J2 = asrootpy(myfile.Get('Had_comp/'+Jet2L)).Clone()
#                plotter.set_histo_style(J2, title=' < 3 jets', drawstyle='hist', color='cyan', linestyle='solid')
#                to_draw.append(J2)
#                to_draw_all.append(J2)
#
#                J3 = asrootpy(myfile.Get('Had_comp/'+Jet3)).Clone()
#                plotter.set_histo_style(J3, title='3 jets', drawstyle='hist', color='red', linestyle='solid')
#                to_draw.append(J3)
#                to_draw_all.append(J3)
#
#                J4 = asrootpy(myfile.Get('Had_comp/'+Jet4)).Clone()
#                plotter.set_histo_style(J4, title='4 jets', drawstyle='hist', color='blue', linestyle='solid')
#                to_draw.append(J4)
#                to_draw_all.append(J4)
#
#                J5P = asrootpy(myfile.Get('Had_comp/'+Jet5P)).Clone()
#                plotter.set_histo_style(J5P, title='5+ jets', drawstyle='hist', color='green', linestyle='solid')
#                to_draw.append(J5P)
#                to_draw_all.append(J5P)
#
#        to_draw.sort()
#
#        plotter.overlay(to_draw_all, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='t_{h} p_{T} [GeV]', ytitle='A.U.', drawstyle='hist')
#        plotter.save('Hadronic_Resolved_Ratio_log_vs_thadpt')
#
##        Unity = total/total
#        J2_Ratio = J2/total
#        J3_Ratio = J3/total
#        J4_Ratio = J4/total
#        J5P_Ratio = J5P/total
#
##        Ratio_Hists.append(Unity)
#        Ratio_Hists.append(J2_Ratio)
#        Ratio_Hists.append(J3_Ratio)
#        Ratio_Hists.append(J4_Ratio)
#        Ratio_Hists.append(J5P_Ratio)
#
#        Ratio_Hists.sort()
#
#        plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction', drawstyle='hist')
#        plotter.save('Hadronic_Resolved_Ratio_vs_thadpt')
#
#        for i in range(len(to_draw)):
#            to_draw[i].SetFillStyle(1001)
#            Add_Ratio_Hists.append(to_draw[i])
#
#        stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total)#, Add_Ratio_Hists[4]/total)
#        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
#        plotter.save('Hadronic_Resolved_Ratio_vs_thadpt_stack')
#
#
#### partial merge
#        Hadronic_Merge_BHadW_Ratios = [
#                ('Gen_Merged_BHadWJet_DRP8_vs_thadpt', 'Gen_Merged_BHadWJet_DRP8_vs_thadpt_2LJ', 'Gen_Merged_BHadWJet_DRP8_vs_thadpt_3J', 'Gen_Merged_BHadWJet_DRP8_vs_thadpt_4J', 'Gen_Merged_BHadWJet_DRP8_vs_thadpt_5PJ', 'Merged b_{h} and W jet #DeltaR < 0.8')
##                 ('Gen_Merged_WJets_DRP8_vs_thadpt', 'Gen_Merged_WJets_DRP8_vs_thadpt_3J', 'Gen_Merged_WJets_DRP8_vs_thadpt_4J', 'Gen_Merged_WJets_DRP8_vs_thadpt_5PJ', 'Merged W jets #DeltaR < 0.8')
##                 ('Gen_Merged_THad_DRP8_vs_thadpt', 'Gen_Merged_THad_DRP8_vs_thadpt_3J', 'Gen_Merged_THad_DRP8_vs_thadpt_4J', 'Gen_Merged_THad_DRP8_vs_thadpt_5PJ', 'All 3 Merged #DeltaR < 0.8')
#        ]
#
##        for name, title, col in DRvals:
#
#        to_draw = []
#        to_draw_all = []
#        Add_Ratio_Hists = []
#        Ratio_Hists = []
#
#        for All, Jet2, Jet3, Jet4, Jet5P, legends in Hadronic_Merge_BHadW_Ratios:
#            total = asrootpy(myfile.Get('Had_comp/'+All)).Clone()
#            plotter.set_histo_style(total, title='All jets', drawstyle='hist', color='black', linestyle='solid')
#            to_draw_all.append(total)
#
#            J2 = asrootpy(myfile.Get('Had_comp/'+Jet2)).Clone()
#            plotter.set_histo_style(J2, title=' < 3 jets', drawstyle='hist', color='cyan', linestyle='solid')
#            to_draw.append(J2)
#            to_draw_all.append(J2)
#
#            J3 = asrootpy(myfile.Get('Had_comp/'+Jet3)).Clone()
#            plotter.set_histo_style(J3, title='3 jets', drawstyle='hist', color='red', linestyle='solid')
#            to_draw.append(J3)
#            to_draw_all.append(J3)
#
#            J4 = asrootpy(myfile.Get('Had_comp/'+Jet4)).Clone()
#            plotter.set_histo_style(J4, title='4 jets', drawstyle='hist', color='blue', linestyle='solid')
#            to_draw.append(J4)
#            to_draw_all.append(J4)
#
#            J5P = asrootpy(myfile.Get('Had_comp/'+Jet5P)).Clone()
#            plotter.set_histo_style(J5P, title='5+ jets', drawstyle='hist', color='green', linestyle='solid')
#            to_draw.append(J5P)
#            to_draw_all.append(J5P)
#
#        to_draw.sort()
#
#        plotter.overlay(to_draw_all, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='t_{h} p_{T} [GeV]', ytitle='A.U.', drawstyle='hist')
#        plotter.save('Hadronic_Merged_BHadWjet_Ratio_log_vs_thadpt')
#
##        Unity = total/total
#        J2_Ratio = J2/total
#        J3_Ratio = J3/total
#        J4_Ratio = J4/total
#        J5P_Ratio = J5P/total
#
##        Ratio_Hists.append(Unity)
#        Ratio_Hists.append(J2_Ratio)
#        Ratio_Hists.append(J3_Ratio)
#        Ratio_Hists.append(J4_Ratio)
#        Ratio_Hists.append(J5P_Ratio)
#
#        Ratio_Hists.sort()
#
#
#        plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction', drawstyle='hist')
#        plotter.save('Hadronic_Merged_BHadWjet_Ratio_vs_thadpt')
#
#        for i in range(len(to_draw)):
#            to_draw[i].SetFillStyle(1001)
#            Add_Ratio_Hists.append(to_draw[i])
#
#        stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total)#, Add_Ratio_Hists[4]/total)
#        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
#        plotter.save('Hadronic_Merged_BHadWjet_Ratio_vs_thadpt_stack')
#
#
#
#### partial merge
#        Hadronic_Merge_Wjets_Ratios = [
#                 ('Gen_Merged_WJets_DRP8_vs_thadpt', 'Gen_Merged_WJets_DRP8_vs_thadpt_2LJ', 'Gen_Merged_WJets_DRP8_vs_thadpt_3J', 'Gen_Merged_WJets_DRP8_vs_thadpt_4J', 'Gen_Merged_WJets_DRP8_vs_thadpt_5PJ', 'Merged W jets #DeltaR < 0.8')
##                 ('Gen_Merged_THad_DRP8_vs_thadpt', 'Gen_Merged_THad_DRP8_vs_thadpt_3J', 'Gen_Merged_THad_DRP8_vs_thadpt_4J', 'Gen_Merged_THad_DRP8_vs_thadpt_5PJ', 'All 3 Merged #DeltaR < 0.8')
#        ]
#
#        to_draw = []
#        to_draw_all = []
#        Add_Ratio_Hists = []
#        Ratio_Hists = []
#
#        for All, Jet2, Jet3, Jet4, Jet5P, legends in Hadronic_Merge_Wjets_Ratios:
#                total = asrootpy(myfile.Get('Had_comp/'+All)).Clone()
#                plotter.set_histo_style(total, title='All jets', drawstyle='hist', color='black', linestyle='solid')
#                to_draw_all.append(total)
#
#                J2 = asrootpy(myfile.Get('Had_comp/'+Jet2)).Clone()
#                plotter.set_histo_style(J2, title=' < 3 jets', drawstyle='hist', color='cyan', linestyle='solid')
#                to_draw_all.append(J2)
#                to_draw.append(J2)
#
#                J3 = asrootpy(myfile.Get('Had_comp/'+Jet3)).Clone()
#                plotter.set_histo_style(J3, title='3 jets', drawstyle='hist', color='red', linestyle='solid')
#                to_draw.append(J3)
#                to_draw_all.append(J3)
#
#                J4 = asrootpy(myfile.Get('Had_comp/'+Jet4)).Clone()
#                plotter.set_histo_style(J4, title='4 jets', drawstyle='hist', color='blue', linestyle='solid')
#                to_draw.append(J4)
#                to_draw_all.append(J4)
#
#                J5P = asrootpy(myfile.Get('Had_comp/'+Jet5P)).Clone()
#                plotter.set_histo_style(J5P, title='5+ jets', drawstyle='hist', color='green', linestyle='solid')
#                to_draw.append(J5P)
#                to_draw_all.append(J5P)
#
#        to_draw.sort()
#
#        plotter.overlay(to_draw_all, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='t_{h} p_{T} [GeV]', ytitle='A.U.', drawstyle='hist')
#        plotter.save('Hadronic_Merged_Wjets_Ratio_log_vs_thadpt')
#
##        Unity = total/total
#        J2_Ratio = J2/total
#        J3_Ratio = J3/total
#        J4_Ratio = J4/total
#        J5P_Ratio = J5P/total
#
##        Ratio_Hists.append(Unity)
#        Ratio_Hists.append(J2_Ratio)
#        Ratio_Hists.append(J3_Ratio)
#        Ratio_Hists.append(J4_Ratio)
#        Ratio_Hists.append(J5P_Ratio)
#
#        Ratio_Hists.sort()
#
#
#        plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction', drawstyle='hist')
#        plotter.save('Hadronic_Merged_Wjets_Ratio_vs_thadpt')
#
#        for i in range(len(to_draw)):
#            to_draw[i].SetFillStyle(1001)
#            Add_Ratio_Hists.append(to_draw[i])
#
#        stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total)#, Add_Ratio_Hists[4]/total)
#        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
#        plotter.save('Hadronic_Merged_Wjets_Ratio_vs_thadpt_stack')
#
#
#
#### full merge
#        Hadronic_Full_Merge_Ratios = [
#                 ('Gen_Merged_THad_DRP8_vs_thadpt', 'Gen_Merged_THad_DRP8_vs_thadpt_2LJ', 'Gen_Merged_THad_DRP8_vs_thadpt_3J', 'Gen_Merged_THad_DRP8_vs_thadpt_4J', 'Gen_Merged_THad_DRP8_vs_thadpt_5PJ', 'All 3 Merged #DeltaR < 0.8')
#        ]
#
#        to_draw = []
#        to_draw_all = []
#        Add_Ratio_Hists = []
#        Ratio_Hists = []
#
#        for All, Jet2, Jet3, Jet4, Jet5P, legends in Hadronic_Full_Merge_Ratios:
#                total = asrootpy(myfile.Get('Had_comp/'+All)).Clone()
#                plotter.set_histo_style(total, title='All jets', drawstyle='hist', color='black', linestyle='solid')
#                to_draw_all.append(total)
#
#                J2 = asrootpy(myfile.Get('Had_comp/'+Jet2)).Clone()
#                plotter.set_histo_style(J2, title=' < 3 jets', drawstyle='hist', color='cyan', linestyle='solid')
#                to_draw.append(J2)
#                to_draw_all.append(J2)
#
#                J3 = asrootpy(myfile.Get('Had_comp/'+Jet3)).Clone()
#                plotter.set_histo_style(J3, title='3 jets', drawstyle='hist', color='red', linestyle='solid')
#                to_draw.append(J3)
#                to_draw_all.append(J3)
#
#                J4 = asrootpy(myfile.Get('Had_comp/'+Jet4)).Clone()
#                plotter.set_histo_style(J4, title='4 jets', drawstyle='hist', color='blue', linestyle='solid')
#                to_draw.append(J4)
#                to_draw_all.append(J4)
#
#                J5P = asrootpy(myfile.Get('Had_comp/'+Jet5P)).Clone()
#                plotter.set_histo_style(J5P, title='5+ jets', drawstyle='hist', color='green', linestyle='solid')
#                to_draw.append(J5P)
#                to_draw_all.append(J5P)
#
#        to_draw.sort()
#
#        plotter.overlay(to_draw_all, legend_def=LegendDefinition(position='NE'), legendstyle='l', logy=True, y_range=(10**(0),10**7), xtitle='t_{h} p_{T} [GeV]', ytitle='A.U.', drawstyle='hist')
#        plotter.save('Hadronic_Merged_THad_Ratio_log_vs_thadpt')
#
##        Unity = total/total
#        J2_Ratio = J2/total
#        J3_Ratio = J3/total
#        J4_Ratio = J4/total
#        J5P_Ratio = J5P/total
#
##        Ratio_Hists.append(Unity)
#        Ratio_Hists.append(J2_Ratio)
#        Ratio_Hists.append(J3_Ratio)
#        Ratio_Hists.append(J4_Ratio)
#        Ratio_Hists.append(J5P_Ratio)
#
#        Ratio_Hists.sort()
#
#
#        plotter.overlay(Ratio_Hists, legend_def=LegendDefinition(position='NW'), legendstyle='l', y_range=(0,1.4), xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction', drawstyle='hist')
#        plotter.save('Hadronic_Merged_THad_Ratio_vs_thadpt')
#
#        for i in range(len(to_draw)):
#            to_draw[i].SetFillStyle(1001)
#            Add_Ratio_Hists.append(to_draw[i])
#
#        stack = plotter.create_stack(Add_Ratio_Hists[0]/total, Add_Ratio_Hists[1]/total, Add_Ratio_Hists[2]/total, Add_Ratio_Hists[3]/total)#, Add_Ratio_Hists[4]/total)
#        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='t_{h} p_{T} [GeV]', ytitle='Event Fraction')
#        plotter.save('Hadronic_Merged_THad_Ratio_vs_thadpt_stack')


#############################################################

