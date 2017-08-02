'''
ttbar_reco_3J Analyzer Plotter macro
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

parser = argparse.ArgumentParser(description='Create plots using files from ttbar_reco_3J.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, AtoTT_M...).')
parser.add_argument('plot', help='Choose type of plots to generate (Gen_Plots, Reco_Plots, Resolution_Plots).')
args = parser.parse_args()

if args.analysis == "Test":
#    if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
    print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
    myfile = root_open('../%s.ttbar_reco_3J.test.root' % args.sample, 'read')
    normfile = views.NormalizeView(root_open('../%s.ttbar_reco_3J.test.root' % args.sample, 'read'))#normalized file
#	lumifile = open('../inputs/%s/%s.lumi' % (jobid, args.sample), 'read')
    plotter = BasePlotter(
    	'plots/ttbar_reco_3J/%s/%s/%s' % (jobid, args.analysis, args.sample),
    	defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}}
    	#defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
    )

elif args.analysis == "Full":
#	if args.sample == "ttJetsM0" or args.sample == "ttJetsM700" or args.sample == "ttJetsM1000":
    print( 'Your analysis: sample are %s: %s' % (args.analysis, args.sample) )
    myfile = root_open('../results/%s/ttbar_reco_3J/%s.root' % (jobid, args.sample), 'read')
    normfile = views.NormalizeView(root_open('../results/%s/ttbar_reco_3J/%s.root' % (jobid, args.sample), 'read'))
    lumifile = open('../inputs/%s/%s.lumi' % (jobid, args.sample), 'read')
    plotter = BasePlotter(
    	'plots/ttbar_reco_3J/%s/%s/%s' % (jobid, args.analysis, args.sample),
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


##############################################################################################

if args.plot == "Gen_Plots":

    kinvar = {'Mass' : ' Mass [GeV]', 'Pt' : ' p_{T}', 'Eta' : ' #eta', 'Costh' : ' Cos(#theta)'}

    for kin in kinvar:
        sdir = args.plot+'/'+kin
        plotter.set_subdir(sdir)

        xlabel = kinvar[kin]

    	Var = [
            ('thad_%s' % kin, 't_{h}'),
            ('tlep_%s' % kin, 't_{l}'),
            ('ttbar_%s' % kin, 't#bar{t}')
        ]

        for var, obj in Var:
            hist = asrootpy(myfile.Get(args.plot+'/'+var)).Clone()
            if hist.Integral() == 0:
                print var
                continue
    
            mean = hist.GetMean()
            rms = hist.GetRMS()
    
            plotter.set_histo_style(hist, color=defcol, xtitle='Gen '+obj+xlabel, ytitle=defyax)
            plotter.plot(hist, drawstyle='hist')
            box = plotter.make_text_box('Mean = %.2f\nRMS = %.2f' % (mean,rms), position='NE')
            box.Draw()
            plotter.save(var)
    
    print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/ttbar_reco_3J/%s/%s/%s/%s .\n' % (jobid, args.analysis, args.sample, args.plot)



##############################################################################################


if args.plot == "Reco_Plots":

    lumi = float(lumifile.readline()) #luminosity from lumi file

    Categories = ['RIGHT', 'MERGE_SWAP', 'MERGE', 'WRONG']
    Objects = {'THad' : 't_{h}', 'TLep' : 't_{l}', 'TTbar' : 't#bar{t}'}
    CUT = {'LCut' : ['< 2', 'black'], 'GCut' : ['#geq 2', 'red']}
    kinvar = {'Mass' : ' Mass [GeV]', 'Pt' : ' p_{T}', 'Eta' : ' #eta', 'Costh' : ' Cos(#theta)'}

        ### plots for separate categories (RIGHT, ...)
    for obj in Objects:
        for cat in Categories:
            for kin in kinvar:
                sdir = args.plot+'/'+cat+'/'+kin
                plotter.set_subdir(sdir)
       
                to_draw = [] 
                xlabel = kinvar[kin]
                for cut in CUT:
                    legends = CUT[cut][0]
                    col = CUT[cut][1]

#                    print '%s_%s_%s_%s' % (obj, kin, cat, cut)
                    Hists = [
                        ('%s_%s_%s_%s' % (obj, kin, cut, cat))
                    ]
        
                    for var in Hists:
                        hist = asrootpy(myfile.Get(args.plot+'/'+var)).Clone()
                        if hist.Integral() == 0:
#                            print var
                            continue
    
                        mean = hist.GetMean()
                        rms = hist.GetRMS()
                        
                        if kin == 'Mass' or kin == 'Pt':
                            plotter.set_histo_style(hist, color=col, title=legends+' Mean = %.0f, RMS = %.0f' % (mean,rms) 
)
                        if kin == 'Mass' and obj == 'THad':
                            plotter.set_histo_style(hist, color=col)
                            plotter.plot(hist, drawstyle='hist')
                            r = hist.Fit("gaus", "S")
                            Mean = r.Parameter(1)
                            sigma = r.Parameter(2)
                            hist.SetTitle(legends+' Fit Mean = %.0f, #sigma = %.0f' % (Mean, sigma))
#                            from pdb import set_trace; set_trace()

                        if kin == 'Eta' or kin == 'Costh':
                            plotter.set_histo_style(hist, color=col, title=legends+' Mean = %.2f, RMS = %.2f' % (mean,rms) )
    
                        plotter.plot(hist, drawstyle='hist')
                        to_draw.append(hist)
    
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='Reco '+Objects[obj]+xlabel, drawstyle='hist')
                box = plotter.make_text_box(cat, position='NE')
                box.Draw()
                plotter.save(obj+'_'+cat+'_'+kin+'_Combined_Cut')


        ### plots comparing LCut, GCut
    for obj in Objects:
        for kin in kinvar:
            sdir = args.plot+'/'+kin
            plotter.set_subdir(sdir)
       
            to_draw = []
            xlabel = kinvar[kin]

            Tot_hist = 0        
            for cut in CUT:
                legends = CUT[cut][0]
                col = CUT[cut][1]

                Cut_hist = 0        

                for cat in Categories:
                    Hists = [
                        ('%s_%s_%s_%s' % (obj, kin, cut, cat))
                    ]

                    for var in Hists:
                        hist = asrootpy(myfile.Get(args.plot+'/'+var)).Clone()

                        if hist.Integral() == 0:
                            continue

                        Cut_hist += hist
                Tot_hist += Cut_hist

                mean = Cut_hist.GetMean()
                rms = Cut_hist.GetRMS()

                if kin == 'Mass' or kin == 'Pt':
                    plotter.set_histo_style(Cut_hist, color=col, title=legends+' Mean = %.0f, RMS = %.0f' % (mean,rms) )

                if kin == 'Mass' and obj == 'THad':
                    plotter.set_histo_style(Cut_hist, color=col)
                    plotter.plot(Cut_hist, drawstyle='hist')
                    r = Cut_hist.Fit("gaus", "S")
                    Mean = r.Parameter(1)
                    sigma = r.Parameter(2)
                    Cut_hist.SetTitle(legends+' Fit Mean = %.0f, #sigma = %.0f' % (Mean, sigma))
#                            from pdb import set_trace; set_trace()

                if kin == 'Eta' or kin == 'Costh':
                    plotter.set_histo_style(Cut_hist, color=col, title=legends+' Mean = %.2f, RMS = %.2f' % (mean,rms) )
                
                plotter.plot(Cut_hist, drawstyle='hist')
                to_draw.append(Cut_hist)

            plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='Reco '+Objects[obj]+xlabel, drawstyle='hist')
            plotter.save(obj+'_'+kin+'_Combined_Cut')

            plotter.set_histo_style(Tot_hist*(50000/lumi), color='black', xtitle='Reco '+Objects[obj]+xlabel, drawstyle='hist')
            plotter.plot(Tot_hist*(50000/lumi))
#            plotter.plot(Tot_hist*50, drawstyle='hist')
            if kin == 'Mass':
                r = Cut_hist.Fit("gaus", "S")
                Mean = r.Parameter(1)
                sigma = r.Parameter(2)
            plotter.save(obj+'_'+kin+'_Norm')

    print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/ttbar_reco_3J/%s/%s/%s/%s .\n' % (jobid, args.analysis, args.sample, args.plot)
    


##############################################################################################


if args.plot == "Resolution_Plots":

    Categories = ['RIGHT', 'MERGE_SWAP', 'MERGE', 'WRONG']
    Objects = {'thad' : 't_{h}', 'tlep' : 't_{l}', 'ttbar' : 't#bar{t}'}
    CUT = {'LCut' : ['< 2', 'black'], 'GCut' : ['#geq 2', 'red']}
#    kinvar = {'Mass' : ' Mass [GeV]'}
    kinvar = {'Mass' : ' Mass [GeV]', 'Pt' : ' p_{T}', 'Eta' : ' #eta', 'Costh' : ' Cos(#theta)'}

        ### plots for separate categories (RIGHT, ...)
    for obj in Objects:
        for cat in Categories:
            for kin in kinvar:
                sdir = args.plot+'/'+cat+'/'+kin
                plotter.set_subdir(sdir)
       
                to_draw = [] 
                xlabel = kinvar[kin]
                for cut in CUT:
                    legends = CUT[cut][0]
                    col = CUT[cut][1]
    
                    Hists = [
                        ('%s_%s_%s_%s' % (obj, kin, cut, cat))
                    ]
        
                    for var in Hists:
                        hist = asrootpy(myfile.Get(args.plot+'/'+var)).Clone()
                        if hist.Integral() == 0:
                            print var
                            continue
    
                        mean = hist.GetMean()
                        rms = hist.GetRMS()

                        if kin == 'Mass':
                            plotter.set_histo_style(hist, color=col)
                            plotter.plot(hist, drawstyle='hist')
                            r = hist.Fit("gaus", "S")
                            Mean = r.Parameter(1)
                            sigma = r.Parameter(2)
                            hist.SetTitle(legends+' Fit Mean = %.0f, #sigma = %.0f' % (Mean, sigma))
#                            from pdb import set_trace; set_trace()


                        if kin == 'Pt':
#                        if kin == 'Mass' or kin == 'Pt':
                            plotter.set_histo_style(hist, color=col, title=legends+' Mean = %.0f, RMS = %.0f' % (mean,rms) )
                        if kin == 'Eta' or kin == 'Costh':
                            plotter.set_histo_style(hist, color=col, title=legends+' Mean = %.2f, RMS = %.2f' % (mean,rms) )
    
                        plotter.plot(hist, drawstyle='hist')
                        to_draw.append(hist)

    
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=Objects[obj]+' #Delta'+xlabel+' (Gen-Reco)', drawstyle='hist')
                box = plotter.make_text_box(cat, position='NE')
                box.Draw()
                plotter.save(obj+'_'+cat+'_'+kin+'_Combined_Cut')
    

        ### plots comparing LCut, GCut
    for obj in Objects:
        for kin in kinvar:
            sdir = args.plot+'/'+kin
            plotter.set_subdir(sdir)
       
            to_draw = []
            xlabel = kinvar[kin]

            for cut in CUT:
                legends = CUT[cut][0]
                col = CUT[cut][1]

                Cut_hist = 0        

                for cat in Categories:
                    Hists = [
                        ('%s_%s_%s_%s' % (obj, kin, cut, cat))
                    ]

                    for var in Hists:
                        hist = asrootpy(myfile.Get(args.plot+'/'+var)).Clone()
                        if hist.Integral() == 0:
                            continue

                        Cut_hist += hist

                if Cut_hist == 0:
                    continue

                mean = Cut_hist.GetMean()
                rms = Cut_hist.GetRMS()
                if kin == 'Mass':
                    plotter.set_histo_style(Cut_hist, color=col)
                    plotter.plot(Cut_hist, drawstyle='hist')
                    r = Cut_hist.Fit("gaus", "S")
                    Mean = r.Parameter(1)
                    sigma = r.Parameter(2)
                    Cut_hist.SetTitle(legends+' Fit Mean = %.0f, #sigma = %.0f' % (Mean, sigma))
#                            from pdb import set_trace; set_trace()

                if kin == 'Pt':
#                if kin == 'Mass' or kin == 'Pt':
                    plotter.set_histo_style(Cut_hist, color=col, title=legends+' Mean = %.0f, RMS = %.0f' % (mean,rms) )
                if kin == 'Eta' or kin == 'Costh':
                    plotter.set_histo_style(Cut_hist, color=col, title=legends+' Mean = %.2f, RMS = %.2f' % (mean,rms) )
                
                plotter.plot(Cut_hist, drawstyle='hist')
                to_draw.append(Cut_hist)

            plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=Objects[obj]+' #Delta'+xlabel+' (Gen-Reco)', drawstyle='hist')
            plotter.save(obj+'_'+kin+'_Combined_Cut')

    print '\ncp -r ~/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/ttbar_reco_3J/%s/%s/%s/%s .\n' % (jobid, args.analysis, args.sample, args.plot)
    

