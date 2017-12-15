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
import functions as fncts
import numpy as np

analyzer = 'ttbar_reco_3J'

parser = argparse.ArgumentParser(description='Create plots using files from ttbar_reco_3J.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, AtoTT_M...).')
parser.add_argument('--topology', default='', help='Choose between merged jet and lost jet events (MERGED, LOST) for plots besides Discriminant.')
parser.add_argument('plot', help='Choose type of plots to generate (Gen, Reconstruction, Resolution, Discriminant).')
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

#if not ( args.topology == 'MERGED' or args.topology == 'LOST' ):
#    print "You chose %s as your event topology.\nYou must choose from the help list!" % args.topology
#    sys.exit()


if not ( args.plot == "Everything" or args.plot == 'Gen'\
        or args.plot == 'Reconstruction' or args.plot == 'Resolution' or args.plot == 'Discriminant'):
    print "You chose %s as your plot type.\nYou must choose from the help list!" % args.plot
    sys.exit()


print( 'Analysis: %s\nSample: %s\nPlot: %s\nTopology: %s' % (args.analysis, args.sample, args.plot, args.topology) )


if args.analysis == "Test":
    myfile = root_open('../%s.%s.test.root' % (args.sample, analyzer), 'read')
    normfile = views.NormalizeView(root_open('../%s.%s.test.root' % (args.sample, analyzer), 'read'))#normalized file
    lumifile = open('../inputs/%s/%s.lumi' % (jobid, args.sample), 'read')

elif args.analysis == "Full":
    myfile = root_open('../results/%s/%s/%s.root' % (jobid, analyzer, args.sample), 'read')
    normfile = views.NormalizeView(root_open('../results/%s/%s/%s.root' % (jobid, analyzer, args.sample), 'read'))
    lumifile = open('../inputs/%s/%s.lumi' % (jobid, args.sample), 'read')

plotter = BasePlotter(
	'plots/%s/%s/%s/%s' % (analyzer, jobid, args.analysis, args.sample),
	defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}}
	#defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
)


def plot_dir(plot):
    print '\ncp -r /uscms/home/jdulemba/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/%s/%s/%s/%s/%s .\n' % (analyzer, jobid, args.analysis, args.sample, plot)

##### Global Var. Definitions ####
defcol = 'black' # default color
defyax = 'A.U.' # yaxis title
if args.sample == "ttJetsM0":
    mass_min = 0
    mass_max = 2000
    nbins = 41
    mtt_bins = np.linspace(mass_min, mass_max, nbins)
if args.sample == "ttJetsM700":
    mass_min = 700
    mass_max = 1000
    nbins = 31
    mtt_bins = np.linspace(mass_min, mass_max, nbins)
if args.sample == "ttJetsM1000":
    mass_min = 1000
    mass_max = 2000
    nbins = 11
    mtt_bins = np.linspace(mass_min, mass_max, nbins)


if 'ttJets' in args.sample:
    decay = 'SM t#bar t'

if 'AtoTT' in args.sample:
    decay = 'A->t#bar t'

if 'HtoTT' in args.sample:
    decay = 'H->t#bar t'


Categories = {'MERGED' : [('RIGHT', 'red', 3345), ('MERGED_SWAP', 'orange', 0), ('MERGED', 'green', 3006), ('WRONG', 'blue', 3354)],\
              'LOST' : [('RIGHT', 'red', 3345), ('LOST_SWAP', 'orange', 0), ('LOST', 'green', 3006), ('WRONG', 'blue', 3354)]}
Objects = {'THad' : 't_{h}', 'TLep' : 't_{l}', 'TTbar' : 't#bar{t}'}
CUT = {'LCut' : ['#lambda_{comb} < 2', 'black'], 'GCut' : ['#lambda_{comb} #geq 2', 'red']}
kinvar = {'Mass' : ' Mass [GeV]', 'Pt' : ' p_{T}', 'Eta' : ' #eta', 'Costh' : ' cos(#theta^{*})'}
DiscPlots = {'Massdisc' : ['#lambda_{mass}'], 'NSdisc' : ['#lambda_{NS}'], 'Totaldisc' : ['#lambda_{comb}']}

##############################################################################################

def Gen_Plots(topology):
    if topology == 'MERGED':
        directory = 'Merged_Plots'

    plots = 'Gen'

        ### create gen plots for thad, tlep, and ttbar objects for all likelihood values (LCUT and GCut) and event types (RIGHT, SWAP, MERGE, WRONG)
    for kvar in kinvar:
        xlabel = kinvar[kvar]

        sdir = directory+'/'+plots+'/'+kvar
        plotter.set_subdir(sdir)

        for obj in Objects:
            gen_hist = 0

            for i in range(len(Categories[topology])):
                evt_type = Categories[topology][i][0]
    
                for disc_cut in CUT:
                    legends = CUT[disc_cut][0]
                    col = CUT[disc_cut][1]
    
                    var = directory+'/'+plots+'/'+evt_type+'/'+disc_cut+'/'+obj+'_'+kvar
                    hist = asrootpy(myfile.Get(var)).Clone()

                    if hist.Integral() == 0:
                        continue
                    else:
                        gen_hist += hist

            mean = gen_hist.GetMean()
            rms = gen_hist.GetRMS()

            if kvar == 'Mass':
                plotter.set_histo_style(gen_hist, color=defcol, xtitle='Gen m_{'+Objects[obj]+'} [GeV]', ytitle=defyax)
                if obj == 'TTbar':
                    gen_hist = RebinView.rebin(gen_hist, mtt_bins)
                    gen_hist.xaxis.range_user = mass_min, mass_max
            elif kvar == 'Pt':
                plotter.set_histo_style(gen_hist, color=defcol, xtitle='Gen '+Objects[obj]+xlabel+' [GeV]', ytitle=defyax)
            else:
                plotter.set_histo_style(gen_hist, color=defcol, xtitle='Gen '+Objects[obj]+xlabel, ytitle=defyax)

            plotter.plot(gen_hist, drawstyle='hist')

            box = plotter.make_text_box('Mean = %.2f\nRMS = %.2f' % (mean,rms), position='NW')
            box.Draw()
            box2 = plotter.make_text_box(decay, position='NE')
            box2.Draw()

            hist_name = '%s_%s_%s_%s_%s' % (args.sample, topology, plots, obj, kvar)

            plotter.save(hist_name)
   
    plot_dir(directory+'/'+plots) 



##############################################################################################


def Reco_Plots(topology):
    if topology == 'MERGED':
        directory = 'Merged_Plots'

    plots = 'Reconstruction'

    lumi = float(lumifile.readline()) #luminosity from lumi file

    rebin_hist = {'Costh' : {'TTbar' : (-1., 1., 26), 'THad' : (-1., 1., 26), 'TLep' : (-1., 1., 26)},\
                  'Eta' : {'TTbar' : (-2.4, 2.4, 51), 'THad' : (-2.4, 2.4, 51), 'TLep' : (-2.4, 2.4, 51)},\
                  'Mass' : {'TTbar' : (200., 2000., 46), 'THad' : (100., 250., 46), 'TLep' : (100., 250., 46)},\
                  'Pt' : {'TTbar' : (0., 1000., 51), 'THad' : (0., 1000., 51), 'TLep' : (0., 1000., 51)}}

        ### Reco plots for thad, tlep, and ttbar objects looking at kinvar dists (mass, costh, ...) based on likelihood value
    for obj in Objects:

        for kvar in kinvar:
            xlabel = kinvar[kvar]

            to_draw = []

            sdir = directory+'/'+plots+'/'+kvar
            plotter.set_subdir(sdir)

            for disc_cut in CUT:
                legends = CUT[disc_cut][0]
                col = CUT[disc_cut][1]

                reco_hist = 0
                draw_types = []
    
                for i in range(len(Categories[topology])):
                    evt_type = Categories[topology][i][0]
                    evt_col = Categories[topology][i][1]
                    evt_fill = Categories[topology][i][2]
    
                    var = directory+'/'+plots+'/'+evt_type+'/'+disc_cut+'/'+obj+'_'+kvar
                    #print obj, kvar, disc_cut, evt_type

                    if (kvar == 'Mass' and obj == 'TLep') or (kvar == 'Costh' and obj == 'TTbar'):
                        continue

                    hist = asrootpy(myfile.Get(var)).Clone()
                    #print var, hist.Integral()

                    if hist.Integral() == 0:
                        continue

                    new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], rebin_hist[kvar][obj][2])
                    hist = RebinView.rebin(hist, new_bins)

#                    print var
                    reco_hist += hist

                    plotter.set_histo_style(hist, color=evt_col, title=evt_type, ytitle=defyax)
                    hist.SetFillStyle(evt_fill)
                    plotter.plot(hist, drawstyle='hist')
                    draw_types.append(hist)

                if reco_hist == 0:
                    continue

                stack, norm_stack, ratio =  fncts.stack_plots(draw_types)

                if kvar == 'Mass':
                    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist', xtitle='Reco m_{'+Objects[obj]+'} [GeV]', ytitle=defyax)
                elif kvar == 'Pt':
                    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist', xtitle='Reco '+Objects[obj]+xlabel+' [GeV]', ytitle=defyax)
                else:
                    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist', xtitle='Reco '+Objects[obj]+xlabel, ytitle=defyax)

                box1 = plotter.make_text_box(decay+'\n'+legends, position='NE')
                box1.Draw()

                hist_cut_name = '%s_%s_%s_%s_LikelihoodVal_%s' % (args.sample, plots, obj, kvar, disc_cut)
                plotter.save(hist_cut_name)
    
                mean = reco_hist.GetMean()
                rms = reco_hist.GetRMS()
                    
                if kvar == 'Mass' or kvar == 'Pt':
                    plotter.set_histo_style(reco_hist, color=col, title=legends+' Mean = %.0f, RMS = %.0f' % (mean,rms) )
                if kvar == 'Mass' and obj == 'THad':
                    plotter.set_histo_style(reco_hist, color=col)
                    plotter.plot(reco_hist, drawstyle='hist')
                    r = reco_hist.Fit("gaus", "S")
                    Mean = r.Parameter(1)
                    sigma = r.Parameter(2)
                    reco_hist.SetTitle(legends+' Fit Mean = %.0f, #sigma = %.0f' % (Mean, sigma))
#                    from pdb import set_trace; set_trace()

                if kvar == 'Eta' or kvar == 'Costh':
                    plotter.set_histo_style(reco_hist, color=col, title=legends+' Mean = %.2f, RMS = %.2f' % (mean,rms) )
    
                plotter.plot(reco_hist, drawstyle='hist')
                to_draw.append(reco_hist)

            if kvar == 'Mass':
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='Reco m_{'+Objects[obj]+'} [GeV]', ytitle=defyax, drawstyle='hist')
            elif kvar == 'Pt':
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='Reco '+Objects[obj]+xlabel+' [GeV]', ytitle=defyax, drawstyle='hist')
            else:
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='Reco '+Objects[obj]+xlabel, ytitle=defyax, drawstyle='hist')

            box2 = plotter.make_text_box(decay, position='NE')
            box2.Draw()

            hist_name = '%s_%s_%s_%s_LikelihoodVal_Combined' % (args.sample, plots, obj, kvar)
            plotter.save(hist_name)
 

    plot_dir(directory+'/'+plots) 



##############################################################################################


def Resolution_Plots(topology):
    if topology == 'MERGED':
        directory = 'Merged_Plots'

    plots = 'Resolution'

    lumi = float(lumifile.readline()) #luminosity from lumi file

    rebin_hist = {'Costh' : {'TTbar' : (-2., 2., 26), 'THad' : (-2., 2., 26), 'TLep' : (-2., 2., 26)},\
                  'Eta' : {'TTbar' : (-4.8, 4.8, 51), 'THad' : (-4.8, 4.8, 51), 'TLep' : (-4.8, 4.8, 51)},\
                  'Mass' : {'TTbar' : (-1000., 2000., 46), 'THad' : (-1000., 500., 46), 'TLep' : (-500., 500., 46)},\
                  'Pt' : {'TTbar' : (-1000., 500., 51), 'THad' : (-500., 500., 51), 'TLep' : (-1000., 1000., 51)}}

        ### Resolution plots for thad, tlep, and ttbar objects looking at kinvar dists (mass, costh, ...) based on likelihood value
    for obj in Objects:

        for kvar in kinvar:
            xlabel = kinvar[kvar]

            to_draw = []

            sdir = directory+'/'+plots+'/'+kvar
            plotter.set_subdir(sdir)

            for disc_cut in CUT:
                legends = CUT[disc_cut][0]
                col = CUT[disc_cut][1]

                reso_hist = 0
                draw_types = []
    
                for i in range(len(Categories[topology])):
                    evt_type = Categories[topology][i][0]
                    evt_col = Categories[topology][i][1]
                    evt_fill = Categories[topology][i][2]
    
                    var = directory+'/'+plots+'/'+evt_type+'/'+disc_cut+'/'+obj+'_'+kvar
                    #print obj, kvar, disc_cut, evt_type

                    if (kvar == 'Mass' and obj == 'TLep') or (kvar == 'Costh' and obj == 'TTbar'):
                        continue

                    hist = asrootpy(myfile.Get(var)).Clone()
                    #print var, hist.Integral()

                    if hist.Integral() == 0:
                        continue

                    new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], rebin_hist[kvar][obj][2])
                    hist = RebinView.rebin(hist, new_bins)

                    reso_hist += hist

                    plotter.set_histo_style(hist, color=evt_col, title=evt_type, ytitle=defyax)
                    hist.SetFillStyle(evt_fill)
                    plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                    draw_types.append(hist)

                if reso_hist == 0:
                    continue

                stack, norm_stack, ratio =  fncts.stack_plots(draw_types)

                if kvar == 'Mass':
                    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist', xtitle='m_{'+Objects[obj]+'} Resolution [GeV]', ytitle=defyax)
                elif kvar == 'Pt':
                    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist', xtitle=Objects[obj]+xlabel+' Resolution [GeV]', ytitle=defyax)
                else:
                    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist', xtitle=Objects[obj]+xlabel+' Resolution', ytitle=defyax)

                box1 = plotter.make_text_box(decay+'\n'+legends, position='NE')
                box1.Draw()

                hist_cut_name = '%s_%s_%s_%s_LikelihoodVal_%s' % (args.sample, plots, obj, kvar, disc_cut)
                plotter.save(hist_cut_name)
    
                mean = reso_hist.GetMean()
                rms = reso_hist.GetRMS()
                    
                if kvar == 'Mass' or kvar == 'Pt':
                    plotter.set_histo_style(reso_hist, color=col, title=legends+' Mean = %.0f, RMS = %.0f' % (mean,rms) )
                if kvar == 'Mass' and obj == 'THad':
                    plotter.set_histo_style(reso_hist, color=col)
                    plotter.plot(reso_hist, drawstyle='hist')
                    r = reso_hist.Fit("gaus", "S")
                    Mean = r.Parameter(1)
                    sigma = r.Parameter(2)
                    reso_hist.SetTitle(legends+' Fit Mean = %.0f, #sigma = %.0f' % (Mean, sigma))
#                    from pdb import set_trace; set_trace()

                if kvar == 'Eta' or kvar == 'Costh':
                    plotter.set_histo_style(reso_hist, color=col, title=legends+' Mean = %.2f, RMS = %.2f' % (mean,rms) )
    
                plotter.plot(reso_hist, drawstyle='hist')
                to_draw.append(reso_hist)
    
            if kvar == 'Mass':
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='m_{'+Objects[obj]+'} Resolution [GeV]', ytitle=defyax, drawstyle='hist')
            elif kvar == 'Pt':
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='Reco '+Objects[obj]+xlabel+' Resolution [GeV]', ytitle=defyax, drawstyle='hist')
            else:
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=Objects[obj]+xlabel+' Resolution', ytitle=defyax, drawstyle='hist')

            box2 = plotter.make_text_box(decay, position='NE')
            box2.Draw()

            hist_name = '%s_%s_%s_%s_LikelihoodVal_Combined' % (args.sample, plots, obj, kvar)
            plotter.save(hist_name)
 

    plot_dir(directory+'/'+plots) 

    
#####################################################################################################


def Discriminant_Plots_2D(topology, disc):
    directory = 'Best_Perm_Disc_Plots'

    if topology == 'MERGED':
        hist_save_name = '3J_Merged_Event'
    if topology == 'LOST':
        hist_save_name = '3J_Lost_Event'

    TotalDisc_2D = { 'LeadTopPt' : 'Leading top p_{T} [GeV]', 'SubleadTopPt' : 'Subleading top p_{T} [GeV]', 'AverageTopPt' : 'Average top p_{T} [GeV]'}
        ### Discriminant plots for Mass, NS, and Combined discs looking at evt_type (RIGHT, WRONG, ...)

    for toppt in TotalDisc_2D:
        hname = '3J_bp_%s_vs_%s' % (disc, toppt)  
        xlabel = TotalDisc_2D[toppt]
        ylabel = DiscPlots[disc][0]+' 3 jets'

        #set_trace()

        sdir_base = args.plot+'/'+topology+'/'+disc

        comb_disc_hist = 0

        for i in range(len(Categories[topology])):
            evt_type = Categories[topology][i][0]
            evt_col = Categories[topology][i][1]
            evt_fill = Categories[topology][i][2]

            plotter.set_subdir(sdir_base+'/'+evt_type)
            var = directory+'/'+topology+'/'+evt_type+'/'+hname
            #print var

            hist = asrootpy(myfile.Get(var)).Clone()
            #print var, hist.Integral()

            if hist.Integral() == 0:
                continue

#            set_trace()
            hist_xmean = hist.GetMean()
            hist_xrms = hist.GetRMS()

            plotter.plot(hist)
            hist.Draw('colz')
            hist.xaxis.set_title(xlabel)
            hist.yaxis.set_title(ylabel)

            box1 = plotter.make_text_box(decay+'\n'+evt_type, position='NE')
            box1.Draw()
            box2 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hist_xmean, hist_xrms), position='NW')
            box2.Draw()

            hist_name = '%s_%s_%s_vs_%s_%s' % (args.sample, hist_save_name, disc, toppt, evt_type)
            plotter.save(hist_name)

            comb_disc_hist += hist


        if comb_disc_hist == 0:
            continue

        comb_disc_hist_xmean = comb_disc_hist.GetMean()
        comb_disc_hist_xrms = comb_disc_hist.GetRMS()

        plotter.plot(comb_disc_hist)
        comb_disc_hist.Draw('colz')
        comb_disc_hist.xaxis.set_title(xlabel)
        comb_disc_hist.yaxis.set_title(ylabel)

        box1 = plotter.make_text_box(decay+'\nMean=%.2f\nRMS=%.2f' % (comb_disc_hist_xmean, comb_disc_hist_xrms), position='NE')
        box1.Draw()

        hist_name = '%s_%s_Combined_%s_vs_%s' % (args.sample, hist_save_name, disc, toppt)
        #set_trace()
        plotter.set_subdir(sdir_base)
        plotter.save(hist_name)
    

#####################################################################################################


def Discriminant_Plots(topology):
    directory = 'Best_Perm_Disc_Plots'

    if topology == 'MERGED':
        hist_save_name = '3J_Merged_Event'
    if topology == 'LOST':
        hist_save_name = '3J_Lost_Event'

    TotalDisc_2D = { 'LeadTopPt' : 'Leading top p_{T} [GeV]', 'SubleadTopPt' : 'Subleading top p_{T} [GeV]', 'AverageTopPt' : 'Average top p_{T} [GeV]'}
        ### Discriminant plots for Mass, NS, and Combined discs looking at evt_type (RIGHT, WRONG, ...)

    to_draw = []
    xaxes = []

    for disc in DiscPlots:

        hname = '3J_bp_%s' % disc    
        xlabel = DiscPlots[disc][0]
        xaxes.append(xlabel)

        if disc == 'Totaldisc':
            Discriminant_Plots_2D(topology, disc)
            #set_trace()

        sdir_base = args.plot+'/'+topology+'/'+disc

        draw_comb_discs = []
        comb_disc_hist = 0

        for i in range(len(Categories[topology])):
            evt_type = Categories[topology][i][0]
            evt_col = Categories[topology][i][1]
            evt_fill = Categories[topology][i][2]

            plotter.set_subdir(sdir_base+'/'+evt_type)
            var = directory+'/'+topology+'/'+evt_type+'/'+hname
            #print var

            hist = asrootpy(myfile.Get(var)).Clone()
            #print var, hist.Integral()

            if hist.Integral() == 0:
                continue
            comb_disc_hist += hist

            hist_mean = hist.GetMean()
            hist_rms = hist.GetRMS()
            plotter.set_histo_style(hist, color=evt_col, title=evt_type, xtitle=xlabel+' 3 jets', ytitle=defyax)
            hist.SetFillStyle(evt_fill)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            draw_comb_discs.append(hist)

            box1 = plotter.make_text_box(decay+'\nMean=%.2f\nRMS=%.2f' % (hist_mean, hist_rms), position='NE')
            box1.Draw()

            hist_type_name = '%s_%s_%s_%s' % (args.sample, hist_save_name, disc, evt_type)
            plotter.save(hist_type_name)

        if comb_disc_hist == 0:
            continue

        if topology == 'MERGED':
            plotter.set_histo_style(comb_disc_hist, color='red', title=topology)
            comb_disc_hist.SetFillStyle(3345)
        if topology == 'LOST':
            plotter.set_histo_style(comb_disc_hist, color='blue', title=topology)
            comb_disc_hist.SetFillStyle(3354)
        plotter.plot(comb_disc_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        to_draw.append(comb_disc_hist)

        comb_disc_mean = comb_disc_hist.GetMean()
        comb_disc_rms = comb_disc_hist.GetRMS()

        stack, norm_stack, ratio = fncts.stack_plots(draw_comb_discs)
        plotter.set_subdir(sdir_base)

        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        box2 = plotter.make_text_box(decay+'\nMean=%.2f\nRMS=%.2f' % (comb_disc_mean, comb_disc_rms), position='NE')
        box2.Draw()

        comb_disc_hist_name = '%s_%s_%s' % (args.sample, hist_save_name, disc)
        plotter.save(comb_disc_hist_name+'_Stack')

        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        box2.Draw()
        plotter.save(comb_disc_hist_name+'_Stack_Norm')

    plot_dir(args.plot+'/'+topology) 

    return to_draw, xaxes
    

#####################################################################################################

if args.plot == 'Gen':
    Gen_Plots('MERGED')

if args.plot == 'Reconstruction':
    Reco_Plots('MERGED')

if args.plot == 'Resolution':
    Resolution_Plots('MERGED')

if args.plot == 'Everything':
    Gen_Plots('MERGED')
    Reco_Plots('MERGED')
    Resolution_Plots('MERGED')


if args.plot == 'Discriminant':

    print '\nMERGED\n'
    merged_disc_hists, merged_xaxis = Discriminant_Plots('MERGED')
    #set_trace()
    print '\nLost\n'
    lost_disc_hists, lost_xaxis = Discriminant_Plots('LOST')

    sdir = args.plot
    plotter.set_subdir(sdir)
   
    to_draw = []
    for i in range(len(merged_disc_hists)):
        to_draw.append(merged_disc_hists[i])
        to_draw.append(lost_disc_hists[i])
        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=merged_xaxis[i]+' 3 jets', ytitle=defyax, drawstyle='hist')
        box = plotter.make_text_box(decay, position='NE')
        box.Draw()

        plotter.save(args.sample+'_3J_Combined_'+DiscPlots.keys()[i])
        to_draw = []

    plot_dir(args.plot) 
