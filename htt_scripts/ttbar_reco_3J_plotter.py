'''
ttbar_reco_3J Analyzer Plotter macro
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
from URAnalysis.PlotTools.views.RebinView import RebinView
import argparse
import matplotlib.pyplot as plt


analyzer = 'ttbar_reco_3J'

parser = argparse.ArgumentParser(description='Create plots using files from ttbar_reco_3J.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, AtoTT_M...).')
parser.add_argument('plot', help='Choose type of plots to generate (Gen_Plots, Reco_Plots, Resolution_Plots).')
args = parser.parse_args()


##### check analysis type
if not (args.analysis == "Test" or args.analysis == "Full"):
    print "You chose %s as your analysis type.\nYou must choose Full or Test!" % args.analysis
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
        print "You chose %s as your sample file.\nYou must choose from the files in results/%s/%s!" % (args.sample, jobid, analyzer)
        sys.exit()


if not ( args.plot == "Everything" or args.plot == 'Gen_Plots'\
        or args.plot == 'Reco_Plots' or args.plot == 'Resolution_Plots'):
    print "You chose %s as your plot type.\nYou must choose from the help list!" % args.plot
    sys.exit()


print( 'Analysis: %s\nSample: %s\nPlot: %s' % (args.analysis, args.sample, args.plot) )


if args.analysis == "Test":
    myfile = root_open('../%s.%s.test.root' % (args.sample, analyzer), 'read')
    normfile = views.NormalizeView(root_open('../%s.%s.test.root' % (args.sample, analzyer), 'read'))#normalized file
#	lumifile = open('../inputs/%s/%s.lumi' % (jobid, args.sample), 'read')

elif args.analysis == "Full":
    myfile = root_open('../results/%s/%s/%s.root' % (jobid, analyzer, args.sample), 'read')
    normfile = views.NormalizeView(root_open('../results/%s/%s/%s.root' % (jobid, analyzer, args.sample), 'read'))
    lumifile = open('../inputs/%s/%s.lumi' % (jobid, args.sample), 'read')

plotter = BasePlotter(
	'plots/%s/%s/%s/%s' % (analyzer, jobid, args.analysis, args.sample),
	defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}}
	#defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
)

def stack_plots(lists):
    lists.sort(key=lambda x: x.Integral())
    total = 0
#    ratio_hists = []
    Stack_hists = []

    for i in lists:
        total += i
    for i in lists:
#        ratio_hists.append(i/total)
        i.SetFillStyle(1001)
        Stack_hists.append(i)

    stack = plotting.HistStack()
    norm_stack = plotting.HistStack()
    for i in Stack_hists:
        stack.Add(i)
        norm_stack.Add(i/total)
        
    return stack, norm_stack
    


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


##############################################################################################

def Gen_Plots():
    plots = 'Gen_Plots'

    kinvar = {'Mass' : ' Mass [GeV]', 'Pt' : ' p_{T} [GeV]', 'Eta' : ' #eta', 'Costh' : ' cos(#theta^{*})'}

    for kin in kinvar:
        sdir = 'Gen_Plots/'+kin
        plotter.set_subdir(sdir)

        xlabel = kinvar[kin]

    	Var = [
            ('thad_%s' % kin, 't_{h}'),
            ('tlep_%s' % kin, 't_{l}'),
            ('ttbar_%s' % kin, 't#bar{t}')
        ]

        for var, obj in Var:
            hist = asrootpy(myfile.Get('Gen_Plots/'+var)).Clone()
            if hist.Integral() == 0:
                print var
                continue
    
            mean = hist.GetMean()
            rms = hist.GetRMS()

            if kin == 'Mass':    
                plotter.set_histo_style(hist, color=defcol, xtitle='Gen m_{'+obj+'} [GeV]', ytitle=defyax)
            else:    
                plotter.set_histo_style(hist, color=defcol, xtitle='Gen '+obj+xlabel, ytitle=defyax)
            plotter.plot(hist, drawstyle='hist')

            box = plotter.make_text_box('Mean = %.2f\nRMS = %.2f' % (mean,rms), position='NW')
            box.Draw()
            box2 = plotter.make_text_box(decay, position='NE')
            box2.Draw()

            plotter.save(var)
   
    plot_dir(plots) 



##############################################################################################


def Reco_Plots():
    plots = 'Reco_Plots'

    lumi = float(lumifile.readline()) #luminosity from lumi file

    Categories = ['RIGHT', 'MERGE_SWAP', 'MERGE', 'WRONG']
    Objects = {'THad' : 't_{h}', 'TLep' : 't_{l}', 'TTbar' : 't#bar{t}'}
    CUT = {'LCut' : ['< 2', 'black'], 'GCut' : ['#geq 2', 'red']}
    kinvar = {'Mass' : ' Mass [GeV]', 'Pt' : ' p_{T}', 'Eta' : ' #eta', 'Costh' : ' cos(#theta^{*})'}

        ### plots for separate categories (RIGHT, ...)
    for obj in Objects:
        for cat in Categories:
            for kin in kinvar:
                sdir = 'Reco_Plots/'+cat+'/'+kin
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
                        hist = asrootpy(myfile.Get('Reco_Plots/'+var)).Clone()
                        if hist.Integral() == 0:
#                            print var
                            continue

                        if 'AtoTT' in args.sample:
                            mttbins = [200., 291.13924051,382.27848101, 473.41772152, 564.55696203, 655.69620253, 746.83544304, 837.97468354, 929.11392405, 1020.25316456, 2000.]
                            if kin == 'Mass':
                                if obj == 'TTbar':
                                    hist = RebinView.rebin(hist, mttbins)
                            
    
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
    
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='Reco '+Objects[obj]+xlabel, ytitle=defyax, drawstyle='hist')
                box = plotter.make_text_box(decay+'\n'+cat, position='NE')
                box.Draw()
                #box2 = plotter.make_text_box(decay, position='NE')
                #box2.Draw()

                plotter.save(obj+'_'+cat+'_'+kin+'_Combined_Cut')


        ### plots comparing LCut, GCut
    for obj in Objects:
        for kin in kinvar:
            sdir = 'Reco_Plots/'+kin
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
                        hist = asrootpy(myfile.Get('Reco_Plots/'+var)).Clone()

                        if 'AtoTT' in args.sample:
                            if cut == 'GCut':
                                mttbins = [200., 291.13924051,382.27848101, 473.41772152, 564.55696203, 655.69620253, 746.83544304, 837.97468354, 929.11392405, 1020.25316456, 1293.67088608, 1567.08860759, 1840.50632911, 2000.]
                                mthadbins = [100., 130.3030303, 154.54545455, 178.78787879, 203.03030303, 227.27272727, 250.]
                                #costhbins = [-1., -0.8989899 , -0.7979798 , -0.6969697 , -0.5959596 , -0.49494949, -0.39393939, -0.29292929, -0.19191919, -0.09090909,\
                                #            0.01010101, 0.11111111, 0.21212121, 0.31313131, 0.41414141, 0.51515152, 0.61616162, 0.71717172, 0.81818182, 0.91919192, 1.]
                                costhbins = [-1., -0.7979798 , -0.5959596 , -0.39393939, -0.19191919, 0.01010101, 0.21212121, 0.41414141, 0.61616162, 0.81818182, 1.]

                            else:
                                mttbins = [200., 382.27848101, 473.41772152, 564.55696203, 610.12658228, 655.69620253, 701.26582278, 746.83544304, 769.62025316, 792.40506329, 837.97468354, 883.5443038, 929.11392405, 1020.25316456]
                                mthadbins = [100., 124.24242424, 142.42424242, 154.54545455, 166.66666667, 190.90909091, 233.33333333]
                                #costhbins = [-1., -0.8989899 , -0.7979798 , -0.6969697 , -0.5959596 , -0.49494949, -0.39393939, -0.29292929, -0.19191919, -0.09090909,\
                                #            0.01010101, 0.11111111, 0.21212121, 0.31313131, 0.41414141, 0.51515152, 0.61616162, 0.71717172, 0.81818182, 0.91919192, 1.]
                                costhbins = [-1., -0.7979798 , -0.5959596 , -0.39393939, -0.19191919, 0.01010101, 0.21212121, 0.41414141, 0.61616162, 0.81818182, 1.]

                            if kin == 'Mass':
                                if obj == 'TTbar':
                                    hist = RebinView.rebin(hist, mttbins)
                                    hist.xaxis.range_user = 200., 1000.
                                if obj == 'THad':
                                    hist = RebinView.rebin(hist, mthadbins)

                            if kin == 'Costh':
                                if obj == 'THad' or obj == 'TLep':
                                    hist = RebinView.rebin(hist, costhbins)    
                            
                        if hist.Integral() == 0:
                            continue

                        Cut_hist += hist
                Tot_hist += Cut_hist

                if Cut_hist == 0:
                    continue

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

            if kin == 'Mass':
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='Reco m_{'+Objects[obj]+'} [GeV]', ytitle=defyax, drawstyle='hist')
            else:
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='Reco '+Objects[obj]+xlabel, ytitle=defyax, drawstyle='hist')

            box2 = plotter.make_text_box(decay, position='NE')
            box2.Draw()

            plotter.save(obj+'_'+kin+'_Combined_Cut')

            if Tot_hist == 0:
                continue

            plotter.set_histo_style(Tot_hist*(50000/lumi), color='black', xtitle='Reco '+Objects[obj]+xlabel, ytitle=defyax, drawstyle='hist')
            plotter.plot(Tot_hist*(50000/lumi))
#            plotter.plot(Tot_hist*50, drawstyle='hist')
            if kin == 'Mass':
                r = Cut_hist.Fit("gaus", "S")
                Mean = r.Parameter(1)
                sigma = r.Parameter(2)
            plotter.save(obj+'_'+kin+'_Norm')

    plot_dir(plots) 
    


##############################################################################################


def Resolution_Plots():
    plots = 'Resolution_Plots'

    Categories = ['RIGHT', 'MERGE_SWAP', 'MERGE', 'WRONG']
    Objects = {'thad' : 't_{h}', 'tlep' : 't_{l}', 'ttbar' : 't#bar{t}'}
    CUT = {'LCut' : ['< 2', 'black'], 'GCut' : ['#geq 2', 'red']}
#    kinvar = {'Mass' : ' Mass [GeV]'}
    kinvar = {'Mass' : ' Mass [GeV]', 'Pt' : ' p_{T} [GeV]', 'Eta' : ' #eta', 'Costh' : ' cos(#theta^{*})'}

        ### plots for separate categories (RIGHT, ...)
    for obj in Objects:
        for cat in Categories:
            for kin in kinvar:
                sdir = 'Resolution_Plots/'+cat+'/'+kin
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
                        hist = asrootpy(myfile.Get('Resolution_Plots/'+var)).Clone()
                        if hist.Integral() == 0:
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

                if kin == 'Mass': 
                    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{'+Objects[obj]+'} Resolution [GeV]', ytitle=defyax, drawstyle='hist')

                else: 
                    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=Objects[obj]+' '+xlabel+' Resolution', ytitle=defyax, drawstyle='hist')
                box = plotter.make_text_box(decay+'\n'+cat, position='NE')
                box.Draw()
                #box2 = plotter.make_text_box(decay, position='NE')
                #box2.Draw()

                plotter.save(obj+'_'+cat+'_'+kin+'_Combined_Cut')
    

        ### plots comparing LCut, GCut
    for obj in Objects:
        for kin in kinvar:
            sdir = 'Resolution_Plots/'+kin
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
                        hist = asrootpy(myfile.Get('Resolution_Plots/'+var)).Clone()
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

            if kin == 'Mass':
                #print '\n'+'\n'+obj
                if obj == 'thad':
                    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), x_range=(-300., 200.), legendstyle='l', xtitle='m_{'+Objects[obj]+'} Resolution [GeV]', ytitle=defyax, drawstyle='hist')
                elif obj == 'ttbar':
                    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), x_range=(-500., 1500.), legendstyle='l', xtitle='m_{'+Objects[obj]+'} Resolution [GeV]', ytitle=defyax, drawstyle='hist')
                else:
                    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{'+Objects[obj]+'} Resolution [GeV]', ytitle=defyax, drawstyle='hist')
            else:
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=Objects[obj]+' '+xlabel+' Resolution', ytitle=defyax, drawstyle='hist')

            box2 = plotter.make_text_box(decay, position='NE')
            box2.Draw()
                
            plotter.save(obj+'_'+kin+'_Combined_Cut')

    plot_dir(plots) 
    

#####################################################################################################


if args.plot == 'Gen_Plots':
    Gen_Plots()

if args.plot == 'Reco_Plots':
    Reco_Plots()

if args.plot == 'Resolution_Plots':
    Resolution_Plots()

if args.plot == 'Everything':
    Gen_Plots()
    Reco_Plots()
    Resolution_Plots()

