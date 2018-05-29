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
from rootpy.plotting import views, Graph, Hist, Hist2D
from argparse import ArgumentParser
import URAnalysis.Utilities.roottools as urtools
import ROOT
from URAnalysis.Utilities.tables import latex_table
from URAnalysis.PlotTools.views.RebinView import RebinView
import argparse
import matplotlib.pyplot as plt
import functions as fncts
import numpy as np
from styles import styles
import URAnalysis.Utilities.prettyjson as prettyjson


out_analyzer = 'ttbar_reco_3J'

parser = argparse.ArgumentParser(description='Create plots using files from ttbar_reco_3J.')

jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']

parser.add_argument('--analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('--sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, AtoTT_M...).')
parser.add_argument('--perm', help='Choose between matched_perms and 3-jet event perms (Matched, Event).')
parser.add_argument('--plot', help='Choose type of plots to generate (Gen, Reconstruction, Resolution, Discriminant, All).')
parser.add_argument('--everything', default=False, help='Run Gen, Reco, Resolution, and Discr plots for matched and event perms.')
args = parser.parse_args()


##### check analysis type
if not (args.analysis == "Test" or args.analysis == "Full"):
    print "You chose %s as your analysis type.\nYou must choose Full or Test!" % args.analysis
    sys.exit()

### check perm type
if not ( args.perm == "Matched" or args.perm == "Best_Perm" or args.perm == "Final_Class" ):
    print "You chose %s as your perm type.\nYou must choose from the help list!" % args.perm
    sys.exit()

analyzer = ''
if args.perm == 'Matched':
    analyzer = 'ttbar_matched_perms'
elif args.perm == 'Best_Perm':
    analyzer = 'ttbar_bestperm_solutions'
else: # args.perm == 'Final_Class'
    analyzer = 'ttbar_final_reco_3J'

##### check sample type
results_files = []
results_file = '' 
if args.analysis == "Test":
    results_files = filter(lambda x: '%s.%s.test.root' % (args.sample, analyzer) in x, os.listdir(project))
    results_file = '%s/%s' % (project, results_files[0])    

    plotter = BasePlotter(
    	'%s/htt_scripts/plots/%s/%s/%s/%s' % (project, out_analyzer, jobid, args.analysis, args.sample),
    	#defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
    	defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}},
        styles = {
            'RIGHT' : styles['*RIGHT'],
            'MERGED_SWAP' : styles['*SWAP'],
            'MERGED' : styles['*OTHER'],
            'WRONG' : styles['*WRONG'],
            'LOST_SWAP' : styles['*SWAP'],
            'LOST' : styles['*OTHER'],
            'sample' : styles[args.sample]
        }
    )
    #set_trace()

if args.analysis == "Full":
    results_files = filter(lambda x: '%s.root' % args.sample in x, os.listdir('%s/results/%s/%s' % (project, jobid, analyzer)))
    results_file  = '%s/results/%s/%s/%s' % (project, jobid, analyzer, results_files[0])

    plotter = BasePlotter(
    	'%s/plots/%s/%s/%s' % (project, jobid, out_analyzer, args.sample),
    	#defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
    	defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}},
        styles = {
            'RIGHT' : styles['*RIGHT'],
            'MERGED_SWAP' : styles['*SWAP'],
            'MERGED' : styles['*OTHER'],
            'WRONG' : styles['*WRONG'],
            'LOST_SWAP' : styles['*SWAP'],
            'LOST' : styles['*OTHER'],
            'sample' : styles[args.sample]
        }
    )


allowed_plot_types = ['Gen', 'Reconstruction', 'Resolution', 'Discriminant', 'Alpha_Correction', 'Post_Alpha', 'All']

if not args.plot in allowed_plot_types:
    print "You chose %s as your plot type.\nYou must choose from the help list!" % args.plot
    sys.exit()



myfile = root_open(results_file, 'read')
normfile = views.NormalizeView( root_open(results_file, 'read') )
lumifile = open('%s/inputs/%s/%s.lumi' % (project, jobid, args.sample), 'read')


#set_trace()

def plot_dir(plot):
    print '\ncp -r %s/htt_scripts/plots/%s/%s/%s/%s/%s .\n' % (project, out_analyzer, jobid, args.analysis, args.sample, plot)

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


m_range = plotter.styles['sample']['name']
decay = plotter.styles['sample']['decay']

##########
'''
This section creates dictionaries used for making plots based on different classifications
1. Categories- categorization of (MERGED or LOST) best perms based on their objects compared to the matched perm's
2. Objects- objects used in reconstruction of system variables and their labels
3. CUT- the label and color of each region from a cut made using the combined discriminant
4. kinvar- kinematic variables used to make plots and their axis labels
5. DiscPlots- the 3 different types of discriminants and their axis labels
6. var_types- the objects used for each kinematic variable and their axis labels
'''

Categories = {'MERGED' : ['RIGHT', 'MERGED_SWAP', 'MERGED', 'WRONG'],\
              'LOST' : ['RIGHT', 'LOST_SWAP', 'LOST', 'WRONG']}
#Categories = {'MERGED' : [('RIGHT', 'red', 3345), ('MERGED_SWAP', 'orange', 0), ('MERGED', 'green', 3006), ('WRONG', 'blue', 3354)],\
#              'LOST' : [('RIGHT', 'red', 3345), ('LOST_SWAP', 'orange', 0), ('LOST', 'green', 3006), ('WRONG', 'blue', 3354)]}
Objects = {'THad' : 't_{h}', 'TLep' : 't_{l}', 'TTbar' : 't#bar{t}'}
CUT = {'LCut' : ['#lambda_{comb} < 2', 'black'], 'GCut' : ['#lambda_{comb} #geq 2', 'red']}
kinvar = {'Mass' : ' Mass [GeV]', 'Pt' : ' p_{T}', 'Eta' : ' #eta', 'Costh' : ' cos(#theta^{*})'}
DiscPlots = {'Massdisc' : ['#lambda_{mass}'], 'NSdisc' : ['#lambda_{NS}'], 'Totaldisc' : ['#lambda_{comb}']}

var_types = {}
if args.perm == 'Matched':
    var_types = {'Costh' : {'THad' : 't_{h} cos(#theta^{*})', 'TLep' : 't_{l} cos(#theta^{*})'},\
            'Mass' : {'THad' : 'M(t_{h}) [GeV]', 'TTbar' : 'M(t#bar{t}) [GeV]'},\
            'Pt' : {'THad' : 'p_{T}(t_{h}) [GeV]', 'THad_PT_P' : 't_{h} sin(#theta)', 'THad_PZ_P' : 't_{h} cos(#theta)', 'TLep' : 'p_{T}(t_{l}) [GeV]', 'TTbar' : 'p_{T}(t#bar{t}) [GeV]'},\
            'Eta' : {'THad' : '#eta(t_{h})', 'TLep' : '#eta(t_{l})'}}#, 'TTbar' : 't#bar{t}'} }
else: #elif args.perm == 'Event':
    var_types = {'Costh' : {'THad' : 't_{h} cos(#theta^{*})', 'TLep' : 't_{l} cos(#theta^{*})'},\
            'Mass' : {'THad' : 'M(t_{h}) [GeV]', 'TTbar' : 'M(t#bar{t}) [GeV]'},\
            'Pt' : {'THad' : 'p_{T}(t_{h}) [GeV]', 'TLep' : 'p_{T}(t_{l}) [GeV]', 'TTbar' : 'p_{T}(t#bar{t}) [GeV]'},\
            'Eta' : {'THad' : '#eta(t_{h})', 'TLep' : '#eta(t_{l})'}}#, 'TTbar' : 't#bar{t}'} }
##########

## beginning of plot function definitions


##############################################################################################
'''
This is the function made to create the alpha(mthad) plots for lost-jet events
so a correction can be made to the reconstructed ttbar objects for these events.
'''
def post_alpha_corrections(directory, subdir, topology):

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]


############# plots for reconstructed vars

    reco_hists = {
        'Mass' : {
                    'THad' : { 'Uncorrected' : ('', 50., 200., 'b'), 'Corrected' : ('_Corrected', 100., 200., 'r') },
                    'TTbar' : { 'Uncorrected' : ('', 200., 2000., 'b'), 'Corrected' : ('_Corrected', 200., 2000., 'r') },
                    'Reco_vs_Gen_TTbar' : { 'Uncorrected' : ('', 'Gen M(t#bar{t}) [GeV]', 'Uncorrected Reco M(t#bar{t}) [GeV]') },
                    'Reco_vs_Gen_Corrected_TTbar' : { 'Corrected' : ('', 'Gen M(t#bar{t}) [GeV]', 'Corrected Reco M(t#bar{t}) [GeV]') }
         },
        'Costh' : {
                    'THad' : { 'Uncorrected' : ('', -1., 1., 'b'), 'Corrected' : ('_Corrected', -1., 1., 'r') },
                    'THad_Labframe' : { 'Uncorrected' : ('', -1., 1., 'b'), 'Corrected' : ('_Corrected', -1., 1., 'r') },
                    'Reco_vs_Gen_THad' : { 'Uncorrected' : ('', 'Gen cos(#theta^{*})', 'Reco cos(#theta^{*})') },
                    'Reco_vs_Gen_Corrected_THad' : { 'Corrected' : ('', 'Gen cos(#theta^{*})', 'Corrected Reco cos(#theta^{*})') }
        }
    }


    #for kvar in reco_hists.keys():
    #    plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Reconstruction', kvar]))
   
    #    print '\n\n%s\n' % kvar
 
    #    for obj in reco_hists[kvar].keys():

    #        if var_types[kvar].has_key(obj):
    #            xlabel = 'Reco %s' % var_types[kvar][obj]
    #        elif obj == 'THad_Labframe':
    #            xlabel = 'Reco cos(#theta)(t_{h}) LF'
    #        else:
    #            xlabel = ''

    #        to_draw = []
    #        for corr_type in reco_hists[kvar][obj].keys():
    #            hvar_extension = reco_hists[kvar][obj][corr_type][0]
    #            hname = '/'.join([directory, 'Post_Alpha_Correction', 'Reconstruction', kvar, obj+hvar_extension])

    #            hist = asrootpy(myfile.Get(hname)).Clone()

    #            if hist.Integral() == 0:
    #                continue
  
    #            if hist.DIM == 2:
    #                plotter.plot(hist)
    #                hist.Draw('colz')
    #                hist.set_x_title(reco_hists[kvar][obj][corr_type][1])
    #                hist.set_y_title(reco_hists[kvar][obj][corr_type][2])

    #                plotter.save('%s_%s_%s' % (obj, kvar, corr_type))
    #                set_trace()
    #                continue


    #                ## set hist range, get hist mean and rms 
    #            hvar_xmin = reco_hists[kvar][obj][corr_type][1]
    #            hvar_xmax = reco_hists[kvar][obj][corr_type][2]
    #            hvar_col = reco_hists[kvar][obj][corr_type][3]

    #            hist.xaxis.range_user = hvar_xmin, hvar_xmax
    #            hmean = hist.GetMean()
    #            hrms = hist.GetRMS()

    #                ## set hist style, plot, and save
    #            plotter.set_histo_style(hist, xtitle='%s %s' % (corr_type, xlabel), ytitle=defyax)
    #            plotter.plot(hist, drawstyle='hist')
    #            box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
    #            box1.Draw()
    #            plotter.save('Reco_%s_%s_%s' % (kvar, obj, corr_type))

    #                ## set hist style to be compared with other type
    #            plotter.set_histo_style(hist, color=hvar_col, title='%s Mean=%.2f, RMS=%.2f' % (corr_type, hmean, hrms) )
    #            to_draw.append(hist)


    #            ## compare uncorreced/corrected hists
    #        if not to_draw:
    #            continue
    #        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
    #        plotter.save('Reco_%s_%s_Comparison' % (kvar, obj) )
    #        set_trace()



###### plots for resolution

    #reso_hists = {
    #    'Mass' : {
    #                'THad' : { 'Uncorrected' : ('', -200., 200., 'b'), 'Corrected' : ('_Corrected', -100., 100., 'r') },
    #                'TTbar' : { 'Uncorrected' : ('', -400., 500., 'b'), 'Corrected' : ('_Corrected', -400., 400., 'r') },
    #                'Frac_THad' : { 'Uncorrected' : ('', -1., 1., 'b'), 'Corrected' : ('_Corrected', -0.5, 0.5, 'r') },
    #                'Frac_TTbar' : { 'Uncorrected' : ('', -1., 1., 'b'), 'Corrected' : ('_Corrected', -1., 1., 'r') },
    #                'Reso_MTTbar_vs_Gen_MTTbar' : { 'Uncorrected' : ('', 'Gen M(t#bar{t}) [GeV]', 'Uncorrected M(t#bar{t}) Resolution [GeV]') },
    #                'Reso_MTTbar_Corrected_vs_Gen_MTTbar' : { 'Corrected' : ('', 'Gen M(t#bar{t}) [GeV]', 'Corrected Reco M(t#bar{t}) Resolution [GeV]') }
    #     },
    #    'Costh' : {
    #                'THad' : { 'Uncorrected' : ('', -1., 1., 'b'), 'Corrected' : ('_Corrected', -1., 1., 'r') }
    #    }
    #}

    #for kvar in reso_hists.keys():
    #    plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', kvar]))
   
    #    print '\n\n%s\n' % kvar
 
    #    for obj in reso_hists[kvar].keys():

    #        if var_types[kvar].has_key(obj):
    #            xlabel = '%s Resolution' % var_types[kvar][obj]
    #        elif obj == 'Frac_THad':
    #            xlabel = 'M(t_{h}) Fractional Resolution (G-R/G)'
    #        elif obj == 'Frac_TTbar':
    #            xlabel = 'M(t#bar{t}) Fractional Resolution (G-R/G)'
    #        else:
    #            xlabel = ''

    #        to_draw = []
    #        for corr_type in reso_hists[kvar][obj].keys():
    #            hvar_extension = reso_hists[kvar][obj][corr_type][0]
    #            hname = '/'.join([directory, 'Post_Alpha_Correction', 'Resolution', kvar, obj+hvar_extension])

    #            hist = asrootpy(myfile.Get(hname)).Clone()

    #            if hist.Integral() == 0:
    #                continue
  
    #            if hist.DIM == 2:
    #                plotter.plot(hist)
    #                hist.Draw('colz')
    #                hist.set_x_title(reso_hists[kvar][obj][corr_type][1])
    #                hist.set_y_title(reso_hists[kvar][obj][corr_type][2])

    #                plotter.save('%s_%s_%s' % (obj, kvar, corr_type))
    #                continue


    #                ## set hist range, get hist mean and rms 
    #            hvar_xmin = reso_hists[kvar][obj][corr_type][1]
    #            hvar_xmax = reso_hists[kvar][obj][corr_type][2]
    #            hvar_col = reso_hists[kvar][obj][corr_type][3]

    #            hist.xaxis.range_user = hvar_xmin, hvar_xmax
    #            hmean = hist.GetMean()
    #            hrms = hist.GetRMS()

    #                ## set hist style, plot, and save
    #            plotter.set_histo_style(hist, xtitle='%s %s' % (corr_type, xlabel), ytitle=defyax)
    #            plotter.plot(hist, drawstyle='hist')
    #            box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
    #            box1.Draw()
    #            plotter.save('Reso_%s_%s_%s' % (kvar, obj, corr_type))

    #                ## set hist style to be compared with other type
    #            plotter.set_histo_style(hist, color=hvar_col, title='%s Mean=%.2f, RMS=%.2f' % (corr_type, hmean, hrms) )
    #            to_draw.append(hist)


    #            ## compare uncorreced/corrected hists
    #        if not to_draw:
    #            continue
    #        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
    #        plotter.save('Reso_%s_%s_Comparison' % (kvar, obj) )


    parton_acceptances = {
                '2Partons' :{
                            '1D' : {
                                    'xlabel' : 'M(t#bar{t}) Resolution [GeV], 2 partons',
                                    'title' : 'Reso_MTTbar_2Partons',
                                    'hists' : {
                                                'Uncorrected' : ('Reso_MTTbar_2Partons', 'b'),
                                                'Corrected' : ('Reso_MTTbar_Corrected_2Partons', 'r')
                                    }
                            },
                            '2D' : {
                                    'xlabel' : 'Gen M(t#bar{t}) [GeV], 2 partons',
                                    'ylabel' : 'M(t#bar{t}) Resolution [GeV]',
                                    'hists' : {
                                                'Uncorrected' : 'Reso_MTTbar_vs_Gen_MTTbar_2Partons',
                                                'Corrected' : 'Reso_MTTbar_Corrected_vs_Gen_MTTbar_2Partons'
                                    }
                            }
                },
                '3Partons' :{
                            '1D' : {
                                    'xlabel' : 'M(t#bar{t}) Resolution [GeV], 3 partons',
                                    'title' : 'Reso_MTTbar_3Partons',
                                    'hists' : {
                                                'Uncorrected' : ('Reso_MTTbar_3Partons', 'b'),
                                                'Corrected' : ('Reso_MTTbar_Corrected_3Partons', 'r')
                                    }
                            },
                            '2D' : {
                                    'xlabel' : 'Gen M(t#bar{t}) [GeV], 3 partons',
                                    'ylabel' : 'M(t#bar{t}) Resolution [GeV]',
                                    'hists' : {
                                                'Uncorrected' : 'Reso_MTTbar_vs_Gen_MTTbar_3Partons',
                                                'Corrected' : 'Reso_MTTbar_Corrected_vs_Gen_MTTbar_3Partons'
                                    }
                            }
                },
                '4Partons' :{
                            '1D' : {
                                    'xlabel' : 'M(t#bar{t}) Resolution [GeV], 4 partons',
                                    'title' : 'Reso_MTTbar_4Partons',
                                    'hists' : {
                                                'Uncorrected' : ('Reso_MTTbar_4Partons', 'b'),
                                                'Corrected' : ('Reso_MTTbar_Corrected_4Partons', 'r')
                                    }
                            },
                            '2D' : {
                                    'xlabel' : 'Gen M(t#bar{t}) [GeV], 4 partons',
                                    'ylabel' : 'M(t#bar{t}) Resolution [GeV]',
                                    'hists' : {
                                                'Uncorrected' : 'Reso_MTTbar_vs_Gen_MTTbar_4Partons',
                                                'Corrected' : 'Reso_MTTbar_Corrected_vs_Gen_MTTbar_4Partons'
                                    }
                            }
                },
                'All' : {
                            '2D' : {
                                    'xlabel' : 'Gen M(t#bar{t}) [GeV]',
                                    'ylabel' : 'M(t#bar{t}) Resolution [GeV]',
                                    'hists' : {
                                                'Uncorrected' : 'Reso_MTTbar_vs_Gen_MTTbar',
                                                'Corrected' : 'Reso_MTTbar_Corrected_vs_Gen_MTTbar'
                                    }
                            }
                }
    }

    for nparts in parton_acceptances.keys():
        if nparts == 'All' :
            plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances']))
        else:
            plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances', nparts]))
        #plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances', nparts]))
   
        print '\n\n%s\n' % nparts

        for ndim in parton_acceptances[nparts].keys():

            if ndim == '1D':

                xlabel = parton_acceptances[nparts][ndim]['xlabel']
                title = parton_acceptances[nparts][ndim]['title']
                to_draw = []

                for corr_type in parton_acceptances[nparts][ndim]['hists'].keys():
                    varname = parton_acceptances[nparts][ndim]['hists'][corr_type][0]
                    hname = '/'.join([directory, 'Post_Alpha_Correction', 'Resolution', 'Parton_Acceptance', varname])
 
                    hist = asrootpy(myfile.Get(hname)).Clone()

                    if hist.Integral() == 0:
                        continue

                    hvar_col = parton_acceptances[nparts][ndim]['hists'][corr_type][1]

            #        hist.xaxis.range_user = hvar_xmin, hvar_xmax
                    hmean = hist.GetMean()
                    hrms = hist.GetRMS()

                        ## set hist style, plot, and save
                    plotter.set_histo_style(hist, xtitle='%s %s' % (corr_type, xlabel), ytitle=defyax)
                    plotter.plot(hist, drawstyle='hist')
                    box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
                    box1.Draw()
                    plotter.save(varname)

                        ## set hist style to be compared with other type
                    plotter.set_histo_style(hist, color=hvar_col, title='%s Mean=%.2f, RMS=%.2f' % (corr_type, hmean, hrms) )
                    to_draw.append(hist)


                    ## compare uncorreced/corrected hists
                if not to_draw:
                    continue
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                plotter.save('%s_Comparison' % title )

            elif ndim == '2D':
                
                xlabel = parton_acceptances[nparts][ndim]['xlabel']
                ylabel = parton_acceptances[nparts][ndim]['ylabel']

                for corr_type in parton_acceptances[nparts][ndim]['hists'].keys():
                    varname = parton_acceptances[nparts][ndim]['hists'][corr_type]

                    if nparts == 'All':
                        hname = '/'.join([directory, 'Post_Alpha_Correction', 'Resolution', 'Mass', varname])
                    else:
                        hname = '/'.join([directory, 'Post_Alpha_Correction', 'Resolution', 'Parton_Acceptance', varname])
 
                    hist = asrootpy(myfile.Get(hname)).Clone()

                    if hist.Integral() == 0:
                        continue

                    plotter.plot(hist)
                    hist.Draw('colz')
                    hist.set_x_title(xlabel)
                    hist.set_y_title('%s %s' % (corr_type, ylabel) )

                    plotter.save( varname )

                #    #### split 2D plots into bins of gen M(ttbar), overlaying uncorrected/corrected hists
                #plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances', nparts, 'yproj']))
                #uncorr_varname = parton_acceptances[nparts][ndim]['hists']['Uncorrected']
                #corr_varname = parton_acceptances[nparts][ndim]['hists']['Corrected']

                #if nparts == 'All':
                #    continue

                #uncorr_hname = '/'.join([directory, 'Post_Alpha_Correction', 'Resolution', 'Parton_Acceptance', uncorr_varname])
                #corr_hname = '/'.join([directory, 'Post_Alpha_Correction', 'Resolution', 'Parton_Acceptance', corr_varname])

                #uncorr_hist = asrootpy(myfile.Get(uncorr_hname)).Clone()
                #corr_hist = asrootpy(myfile.Get(corr_hname)).Clone()

                #if uncorr_hist.Integral() == 0 or corr_hist.Integral() == 0:
                #    print 'One of the hists is empty!'
                #    sys.exit()

                #xbins = np.linspace(uncorr_hist.GetXaxis().GetBinLowEdge(1), uncorr_hist.GetXaxis().GetBinUpEdge(uncorr_hist.GetXaxis().GetNbins()), uncorr_hist.GetNbinsX()/2+1) # combine 2 xbins
                #ybins = np.linspace(uncorr_hist.GetYaxis().GetBinLowEdge(1), uncorr_hist.GetYaxis().GetBinUpEdge(uncorr_hist.GetYaxis().GetNbins()), uncorr_hist.GetYaxis().GetNbins()+1 ) # ybinning remains unchanged

                #uncorr_hist = RebinView.newRebin2D(uncorr_hist, xbins, ybins)
                #corr_hist = RebinView.newRebin2D(corr_hist, xbins, ybins)
                #for xbin in range(uncorr_hist.GetNbinsX() + 1):
                #    uncorr_hist_yproj = uncorr_hist.ProjectionY("", xbin, xbin)
                #    corr_hist_yproj = corr_hist.ProjectionY("", xbin, xbin)

                #    if uncorr_hist_yproj.Integral() == 0 and corr_hist_yproj.Integral() == 0:
                #        continue

                #    set_trace()
                #    uncorr_hist_yproj.Draw('hist')
                #    uncorr_hist_yproj.SetXTitle('Uncorrected')
                #    uncorr_hist_yproj.SetYTitle(defyax)
                #    #uncorr_hist_yproj.SetLineColor('b')
                #    #plotter.set_histo_style(uncorr_hist_yproj, xtitle='Uncorrected', ytitle=defyax)
                #    #plotter.plot(hist, drawstyle='hist')

                #    corr_hist_yproj.Draw('hist same')
                #    corr_hist_yproj.SetXTitle('Corrected')
                #    corr_hist_yproj.SetYTitle(defyax)
                #    #corr_hist_yproj.SetLineColor('r')

                #    #mean = hist_yproj.GetMean()
                #    #stddev= hist_yproj.GetStdDev()
                #    genmttbar_range = '%.1f #leq Gen M(t#bar{t}) #leq %.1f' % (uncorr_hist.GetXaxis().GetBinLowEdge(xbin), uncorr_hist.GetXaxis().GetBinUpEdge(xbin))

                #    box1 = plotter.make_text_box('%s' % genmttbar_range, position='NE')
                #    box1.Draw()

                #    #plotter.overlay([uncorr_hist_yproj, corr_hist_yproj], legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=ylabel, ytitle=defyax, drawstyle='hist')
                    
    #            box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
    #            box1.Draw()
    #        plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')

                    #plotter.save('%s_%s' % (varname, hist.GetXaxis().GetBinLowEdge(xbin)) )
                    #set_trace()




def alpha_corrections(directory, subdir, topology):

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]


    hvars = [ # variables for hists
        ('THad_E/Alpha_THad_E', '173.1/Reco M(t_{h})', '#alpha_{E} = Gen E(t_{h})/Reco E(t_{h})', ''),
        ('THad_P/Alpha_THad_P', '173.1/Reco M(t_{h})', '#alpha_{P} = Gen P(t_{h})/Reco P(t_{h})', ''),
    ]

    xranges = np.linspace(0.9, 2.5, 9) ## rebin xaxis so [0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
    for num in range(len(xranges)-1):
        txt_box = "%.1f < 173.1/Reco M(t_{h}) #leq %.1f" % (xranges[num], xranges[num+1])
        Evar = ('THad_E/Gen_vs_Reco_THadE_%.1fto%.1f' % (xranges[num], xranges[num+1]), 'Reco E(t_{h})', 'Gen E(t_{h})', txt_box )
        Pvar = ('THad_P/Gen_vs_Reco_THadP_%.1fto%.1f' % (xranges[num], xranges[num+1]), 'Reco E(t_{h})', 'Gen E(t_{h})', txt_box )
        hvars.append(Evar)
        hvars.append(Pvar)

    fits = {
        'Alpha_THad_E' : {},
        'Alpha_THad_P' : {},
    }

    for hvar, xlabel, ylabel, txt_box_label in hvars:
        hname = '/'.join([directory, 'Alpha_Correction', hvar])
        hist = asrootpy(myfile.Get(hname)).Clone()

        if hist.Integral() == 0:
            continue

        plotter.plot(hist)
        #set_trace()
        if hvar == 'THad_E/Alpha_THad_E' or hvar == 'THad_P/Alpha_THad_P':

            xbins = np.linspace(0.9, 2.5, 9) ## rebin xaxis so [0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
            ybins = np.linspace( hist.GetYaxis().GetBinLowEdge(1), hist.GetYaxis().GetBinUpEdge(hist.GetYaxis().GetNbins()), hist.GetYaxis().GetNbins()+1 ) # ybinning remains unchanged

            hist = RebinView.newRebin2D(hist, xbins, ybins)
            hist.xaxis.range_user = 0.9, 2.5
            #hist.yaxis.range_user = 0.0, 4.0
            #set_trace()

            medians = []
            median_weights = []
            median_errors = []
            maxvals = []
            for xbin in range(hist.GetNbinsX() + 1):
                if hist.GetXaxis().GetBinLowEdge(xbin) >= hist.xaxis.range_user[0] and hist.GetXaxis().GetBinUpEdge(xbin) <= hist.xaxis.range_user[1]:
                    hist_yproj = hist.ProjectionY("", xbin, xbin)
                    if hist_yproj.Integral() == 0:
                        continue

                    hist_yproj.Draw()
                    hist_yproj.GetXaxis().SetRangeUser(0.0,15.0)
                    hist_yproj.GetXaxis().SetTitle(ylabel)

                    probs = np.zeros(len(ybins))
                    binvals = np.zeros(len(ybins))

                    for bins in range(len(ybins)):
                        probs[bins] = hist_yproj.Integral(1, bins+1)/hist_yproj.Integral()
                        binvals[bins] = hist_yproj.GetBinContent(bins)

                    median = hist_yproj.GetBinCenter(np.where(probs>= 0.5)[0][0])
                    median_error = 1.2533*hist_yproj.GetMeanError() #standard error of median
                    if median_error == 0:
                        #set_trace()
                        median_weight = 1000 # 1/(standard error of median) to be used in fit
                    else:
                        median_weight = 1/(1.2533*hist_yproj.GetMeanError()) # 1/(standard error of median) to be used in fit
                    maxval = round(hist_yproj.GetBinCenter(np.where(binvals==binvals.max())[0][0]), 2)

                    medians.append(median)
                    median_weights.append(median_weight)
                    median_errors.append(median_error)
                    maxvals.append(maxval)

                    box1 = plotter.make_text_box('%s #leq %s #leq %s' % (hist.GetXaxis().GetBinLowEdge(xbin), xlabel, hist.GetXaxis().GetBinUpEdge(xbin)), position='NE')
                    box1.Draw()
                    plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', hvar.split('/')[0], 'y_proj']))
                    plotter.save('%s_%s' % (hvar.split('/')[1], hist.GetXaxis().GetBinLowEdge(xbin)) )
                    #set_trace()

            xbin_ints = [hist.Integral(xbinid+1, xbinid+1, 1, hist.GetNbinsY()+1) for xbinid in range(hist.GetNbinsX())]
            xbin_ratios = [xbin_ints[i]/sum(xbin_ints) for i in range(len(xbin_ints))]

            hist.yaxis.range_user = 0.0, 4.0

            #### section to fit polynomials to median values from xaxis bins
            #xfit_bins = np.linspace(hist.GetXaxis().GetXmin(),hist.GetXaxis().GetXmax(), 200) #makes points to be used for fit results
            xfit_bins = np.linspace(0, 5.0, 200) #makes points to be used for fit results
            hist.GetXaxis().GetCenter(xbins) # change array of xbin values from edges to centers

            fit_1d_groom = np.polyfit(xbins[:-1], medians, 1, w=median_weights) # fit medians to 1 degree polynomial for 0.9 < xbin < 2.5
            fit_1d_1p1 = np.polyfit(xbins[1:-1], medians[1:], 1, w=median_weights[1:]) # weighted fit medians to 1 degree polynomial for 1.1 < xbin < 2.5
            fit_2d_groom = np.polyfit(xbins[:-1], medians, 2, w=median_weights) # fit medians to 2 degree polynomial for 0.9 < xbin < 2.5
            fit_2d_1p1 = np.polyfit(xbins[1:-1], medians[1:], 2, w=median_weights[1:]) # fit medians to 2 degree polynomial for 1.1 < xbin < 2.5

            p1_groom = np.poly1d(fit_1d_groom) #gets parameters for 1 degree fit getting rid of first 2 points
            p1_1p1 = np.poly1d(fit_1d_1p1) #gets parameters for weighted 1 degree fit getting rid of first 3 points
            p2_groom = np.poly1d(fit_2d_groom) #gets parameters for 2 degree fit getting rid of first 2 points
            p2_1p1 = np.poly1d(fit_2d_1p1) #gets parameters for weighted 2 degree fit getting rid of first 3 points

            fig = plt.figure()
            fit_medians = plt.errorbar(xbins[:-1], medians, yerr=median_errors, label='Medians', linestyle='None', color='black', marker='.') #plots median values with errors
                # 1 degree fits
            p1_1p1_weighted = plt.plot(xfit_bins, p1_1p1(xfit_bins), label='weighted deg=1 bins > 1.1', linestyle='-', color='red') #plots p1 for weighted bins > 1.1
            p1_groom_fit = plt.plot(xfit_bins, p1_groom(xfit_bins), label='weighted deg=1 bins > 0.9', linestyle='--', color='red') #plots p1 for bins > 0.9
                # 2 degree fits
            p2_1p1_weighted = plt.plot(xfit_bins, p2_1p1(xfit_bins), label='weighted deg=2 bins > 1.1', linestyle='-', color='green') #plots p2 for weighted bins > 1.1
            p2_groom_fit = plt.plot(xfit_bins, p2_groom(xfit_bins), label='weighted deg=2 bins > 0.9', linestyle='--', color='green') #plots p2 for bins > 0.9

            fits[hvar.split('/')[1]]['1D_g1.1'] = fit_1d_1p1.tolist()
            fits[hvar.split('/')[1]]['2D_g1.1'] = fit_2d_1p1.tolist()

            plt.xlabel('173.1/Reco M($t_{h}$)')
            if hvar == 'THad_E/Alpha_THad_E':
                plt.ylabel('$\\alpha_{E}$ = Gen E($t_{h}$)/Reco E($t_{h}$)')
            if hvar == 'THad_P/Alpha_THad_P':
               plt.ylabel('$\\alpha_{P}$ = Gen P($t_{h}$)/Reco P($t_{h}$)')
            plt.ylim(0.0, 10.0)
            plt.grid()
            plt.legend(loc='upper right',fontsize=8, numpoints=1)
            plotter.set_subdir('/'.join([subdir, 'Alpha_Correction']))
            plotting_dir = plotter.outputdir
            fig.savefig('%s/%s_fits' % (plotting_dir, hvar.split('/')[1]))
            #set_trace()

        hist.Draw('colz')
        hist.xaxis.set_title(xlabel)
        hist.yaxis.set_title(ylabel)

        if 'Gen_vs_Reco' in hvar:
            box1 = plotter.make_text_box( txt_box_label, position='NE')
            box1.Draw()
            if 'THadE' in hvar:
                plotter.set_subdir('/'.join([subdir, 'Alpha_Correction/THad_E/Gen_vs_Reco']))
            if 'THadP' in hvar:
                plotter.set_subdir('/'.join([subdir, 'Alpha_Correction/THad_P/Gen_vs_Reco']))
        else:
            plotter.set_subdir('/'.join([subdir, 'Alpha_Correction']))

        plotter.save(hvar.split('/')[1])
        
        #set_trace()

        #### create hists based on perm category (RIGHT/WRONG, etc...)
        #if topology:
        #    for cat in Categories[topology]:
        #        plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', cat]))
        #        #evt_col=plotter.styles[cat]['fillcolor']
        #        evt_type = plotter.styles[cat]['name']
        #        cat_hname = '/'.join([directory, 'Alpha_Correction', cat, hvar])
        #        cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

        #        if cat_hist.Integral() == 0:
        #            continue

        #        plotter.plot(cat_hist)
    
        #        cat_hist.Draw('colz')
        #        cat_hist.xaxis.set_title(xlabel)
        #        cat_hist.yaxis.set_title(ylabel)

        #        #r = cat_hist.Fit("pol1", "Q") # fit hist with polynomial of degree 1 with no output

        #        box1 = plotter.make_text_box( cat, position='NW')
        #        box1.Draw()
    
        #        plotter.save('%s_%s' % (hvar, cat))

    with open('%s/fit_parameters.json' % '/'.join([plotter.base_out_dir, subdir, 'Alpha_Correction']), 'w') as f:
        #set_trace()
        f.write(prettyjson.dumps(fits))



##############################################################################################
'''
These are plots whose permutations are made from 3-jet events (at least 2 pass b-tagging reqs)
and which have at least one solution from using the merged- and lost-jet assumption discriminants.

solutions_dict-either only one solution exists (Only_Merged(Lost)) from the discriminants,
    both solutions exist but give the same classification (Both_BP/Class_Merged(Lost)),
    or all distributions in which a discriminant has a solution (Merged(Lost)_BP)

both_sols_dict- classifications for events in which both discriminants have solutions (some of the possibilities)
    1. disctributions for both classifications for a discriminant assumption (Merged(Lost)_BP/Class_Merged(Lost))
    2. discr solutions give opposite classifications (Opposite_Class/Merged(Lost)_BP/Class...)
'''
def Event_Plots( plot ):
    directory = '3J_Event_Plots'

    #set_trace()
    solutions_dict = {'Only_Merged' : 'MERGED', 'Only_Lost' : 'LOST', 'Merged_BP' : 'MERGED', 'Lost_BP' : 'LOST',\
                    'Both_BP/Class_Merged' : 'MERGED', 'Both_BP/Class_Lost' : 'LOST'}
    for sol in solutions_dict.keys():
        if plot== 'Discriminant':
            print '\nMaking discriminant plots for %s 3-jet events\n\n' % sol
            Discriminant_Plots('%s/%s' % (directory, sol), '%s/%s' % (directory, sol), solutions_dict[sol])
        if plot== 'Reconstruction':
            print '\nMaking reco plots for %s 3-jet events\n\n' % sol
            Reco_Plots('%s/%s' % (directory, sol), '%s/%s' % (directory, sol), solutions_dict[sol])
        if plot== 'Resolution':
            print '\nMaking resolution plots for %s 3-jet events\n\n' % sol
            Resolution_Plots('%s/%s' % (directory, sol), '%s/%s' % (directory, sol), solutions_dict[sol])
        if plot== 'Gen':
            print '\nMaking gen plots for %s 3-jet events\n\n' % sol
            Gen_Plots('%s/%s' % (directory, sol), '%s/%s' % (directory, sol), solutions_dict[sol])


    both_sols_dict = {'Merged_BP' : {'Class_Merged' : 'MERGED', 'Class_Lost' : 'MERGED'}, 'Lost_BP' : {'Class_Merged' : 'LOST', 'Class_Lost' : 'LOST'},\
                    'Opposite_Class/Merged_BP' : {'Class_Merged' : 'MERGED', 'Class_Lost' : 'MERGED'}, 'Opposite_Class/Lost_BP' : {'Class_Merged' : 'LOST', 'Class_Lost' : 'LOST'}}
    for bp in both_sols_dict.keys():
        for classification in both_sols_dict[bp]:
            subdir = 'Both_BP/'+'/'.join([bp, classification])

            if plot== 'Discriminant':
                print '\nMaking discriminant plots for %s 3-jet events\n\n' % subdir
                Discriminant_Plots('%s/%s' % (directory, subdir), '%s/%s' % (directory, subdir), both_sols_dict[bp][classification])
            if plot== 'Reconstruction':
                print '\nMaking reco plots for %s 3-jet events\n\n' % subdir
                Reco_Plots('%s/%s' % (directory, subdir), '%s/%s' % (directory, subdir), both_sols_dict[bp][classification])
            if plot== 'Resolution':
                print '\nMaking resolution plots for %s 3-jet events\n\n' % subdir
                Resolution_Plots('%s/%s' % (directory, subdir), '%s/%s' % (directory, subdir), both_sols_dict[bp][classification])
            if plot== 'Gen':
                print '\nMaking gen plots for %s 3-jet events\n\n' % subdir
                Gen_Plots('%s/%s' % (directory, subdir),'%s/%s' % (directory, subdir), both_sols_dict[bp][classification])


def Final_Reco_Plots( plot ):
    directory = '3J_Event_Plots'

    plot_types = {'Reconstruction' : 'Reco', 'Resolution' : 'Resolution (Gen-Reco)', 'Gen' : 'Gen'}#, 'Discr' : 'Discr'}

    classifications = {'Class_Merged' : 'MERGED', 'Class_Lost' : 'LOST'}
    reco_dirs = {'Clear_Classification' : ['r', 'Clear Class'], 'Final_Classification' : ['b', 'Classes Incl Cut'] }

    reco_subdir_names = ['Clear_Classification', 'Clear_and_MassCut_Classes']

    dirs = [ '/'.join([directory, reco, classes]) for reco in reco_dirs.keys() for classes in classifications.keys()]
    subdirs = ['/'.join([directory, 'Final_Reco', reco, classes]) for reco in reco_subdir_names for classes in classifications.keys() ]

    inputs = [ (dirs[i], subdirs[i]) for i in range(len(dirs)) ]

    #set_trace()

    for reco_dir, reco_subdir in inputs:
        if plot == 'Discriminant':
            print '\nMaking Discr plots for %s 3-jet events\n\n' % reco_dir
            [Discriminant_Plots( reco_dir , reco_subdir, classifications[classes] ) for classes in classifications.keys() if classes in reco_dir]

        if plot == 'Gen':
            print '\nMaking Gen plots for %s 3-jet events\n\n' % reco_dir
            [Gen_Plots( reco_dir , reco_subdir, classifications[classes] ) for classes in classifications.keys() if classes in reco_dir]

        if plot == 'Reconstruction':
            print '\nMaking reco plots for %s 3-jet events\n\n' % reco_dir
            [Reco_Plots( reco_dir , reco_subdir, classifications[classes] ) for classes in classifications.keys() if classes in reco_dir]

        if plot == 'Resolution':
            print '\nMaking resolution plots for %s 3-jet events\n\n' % reco_dir
            [Resolution_Plots( reco_dir , reco_subdir, classifications[classes] ) for classes in classifications.keys() if classes in reco_dir]


        if plot == 'Alpha_Correction':
            ### create alpha correction plots for lost-jet events
            if reco_dir == '3J_Event_Plots/Clear_Classification/Class_Lost' or reco_dir == '3J_Event_Plots/Final_Classification/Class_Lost':
                [alpha_corrections( reco_dir , reco_subdir, classifications[classes] ) for classes in classifications.keys() if classes in reco_dir]

        if plot == 'Post_Alpha':
            ### create post alpha correction plots for lost-jet events
            if reco_dir == '3J_Event_Plots/Clear_Classification/Class_Lost' or reco_dir == '3J_Event_Plots/Final_Classification/Class_Lost':
                [post_alpha_corrections( reco_dir , reco_subdir, classifications[classes] ) for classes in classifications.keys() if classes in reco_dir]


    #plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]
    #gen_rebin_hist = {'Costh' : {'THad' : (-1., 1., 26), 'TLep' : (-1., 1., 26)},\
    #              'Eta' : {'TTbar' : (-2.5, 2.5, 26), 'THad' : (-2.5, 2.5, 26), 'TLep' : (-2.5, 2.5, 26)},\
    #              'Mass' : {'TTbar' : (200., 2000., 46), 'THad' : (100., 250., 46)},\
    #              'Pt' : {'TTbar' : (0., 1000., 51), 'THad' : (0., 1000., 51), 'TLep' : (0., 1000., 51)}}
    #reco_rebin_hist = {'Costh' : {'THad' : (-1., 1., 26), 'TLep' : (-1., 1., 26)},\
    #              'Eta' : {'TTbar' : (-2.5, 2.5, 26), 'THad' : (-2.5, 2.5, 26), 'TLep' : (-2.5, 2.5, 26)},\
    #              'Mass' : {'TTbar' : (200., 2000., 46), 'THad' : (50., 250., 46)},\
    #              'Pt' : {'TTbar' : (0., 1000., 51), 'THad' : (0., 1000., 51), 'TLep' : (0., 1000., 51)}}
    #reso_rebin_hist = {'Costh' : {'TTbar' : (-2., 2., 26), 'THad' : (-2., 2., 26), 'TLep' : (-2., 2., 26)},\
    #              'Eta' : {'TTbar' : (-5.0, 50, 26), 'THad' : (-5.0, 5.0, 26), 'TLep' : (-5.0, 5.0, 26)},\
    #              'Mass' : {'TTbar' : (-1000., 2000., 46), 'THad' : (-1000., 500., 46), 'TLep' : (-500., 500., 46)},\
    #              'Pt' : {'TTbar' : (-1000., 500., 51), 'THad' : (-500., 500., 51), 'TLep' : (-1000., 1000., 51)}}

    #    ### Reco and Res plots for thad, tlep, and ttbar objects looking at kinvar dists (mass, costh, ...)
    #for plot_type in plot_types.keys():
    #    for kvar in var_types.keys():
    #        for obj in var_types[kvar].keys():
    #            if plot_type == 'Reconstruction' or plot_type == 'Gen':
    #                xlabel = '%s %s' % (plot_types[plot_type], var_types[kvar][obj])
    #            elif plot_type == 'Resolution':
    #                xlabel = '%s %s' % (var_types[kvar][obj], plot_types[plot_type])
    #    
    #            for classes in classifications:
    #                plotter.set_subdir('/'.join([directory, 'Final_Reco/Comparison', classes, plot_type, kvar]))
    #    
    #                to_draw = []
    #    
    #                for reco in reco_dirs.keys():
    #                    evt_col=reco_dirs[reco][0]
    #            #        evt_type = plotter.styles[cat]['name']
    #                    hname = '/'.join([directory, reco, classes, plot_type, kvar, obj])
    #                    hist = asrootpy(myfile.Get(hname)).Clone()

    #                    if plot_type == 'Gen':    
    #                        new_bins = np.linspace(gen_rebin_hist[kvar][obj][0], gen_rebin_hist[kvar][obj][1], hist.nbins()/2+1)
    #                    if plot_type == 'Reconstruction':    
    #                        new_bins = np.linspace(reco_rebin_hist[kvar][obj][0], reco_rebin_hist[kvar][obj][1], hist.nbins()/2+1)
    #                    if plot_type == 'Resolution':    
    #                        new_bins = np.linspace(reso_rebin_hist[kvar][obj][0], reso_rebin_hist[kvar][obj][1], hist.nbins()/2+1)
    #                    hist = RebinView.rebin(hist, new_bins)
    #                    if hist.Integral() == 0:
    #                        continue
    #                    mean = hist.GetMean()
    #                    rms = hist.GetRMS()
    #    
    #                    plotter.set_histo_style(hist, color=evt_col, title=reco_dirs[reco][1]+' Mean = %.2f, RMS = %.2f' % (mean,rms) )
    #                    to_draw.append(hist)
    #    
    #                if not to_draw:
    #                    continue
    #                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
    #    
    #                plotter.save('Final_Reco_%s_%s_%s_Comparison' % (classes, obj, kvar) )

##############################################################################################
'''
These are plots for 3-jet events in which a matched permutation was created and whose distributions
can be compared based on its event type and discriminant assumption

event_types- matched perm is considered Merged or Lost based on its objects
treat_types- lost and merged-jet discriminant assumptions are applied
'''


def Matched_Perm_Plots( plot_type ):
    event_types = {'Merged' : 'r', 'Lost' : 'b'}
    treat_types = ['TreatMerged', 'TreatLost']

    directory = 'Matched_Perm_Plots'

    if plot_type == 'Discriminant':
        for treatment in treat_types:
            plotter.set_subdir('/'.join([directory, treatment, 'Discr']))
            for disc in DiscPlots:
                xlabel = DiscPlots[disc][0]

                to_draw = []
                for evt_type in event_types.keys():
                    var = '/'.join([directory, evt_type, treatment, 'Discr', disc])

                    hist = asrootpy(myfile.Get(var)).Clone()
                    if hist.Integral() == 0:
                        continue
                    mean = hist.GetMean()
                    rms = hist.GetRMS()

                    plotter.set_histo_style(hist, color=event_types[evt_type], title=evt_type+' Mean = %.2f, RMS = %.2f' % (mean,rms) )
                    to_draw.append(hist)
    
                plotter.defaults['watermark'] = ['%s Comparison (13 TeV, 25ns)' % treatment, False]
                plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
    
                plotter.save('%s_%s_Comparison' % (disc, treatment))
                #set_trace()

    else:
        for treatment in treat_types:
            for kvar in var_types.keys():
                plotter.set_subdir('/'.join([directory, treatment, plot_type, kvar]))
                for obj in var_types[kvar].keys():
                    if plot_type == 'Gen':
                        xlabel = 'Gen %s' % var_types[kvar][obj]
                    elif plot_type == 'Reconstruction':
                        xlabel = 'Reco %s' % var_types[kvar][obj]
                    elif plot_type == 'Resolution':
                        xlabel = '%s Resolution (Gen-Reco)' % var_types[kvar][obj]
                    to_draw = []
                    #set_trace()
                    for evt_type in event_types.keys():
                        var = '/'.join([directory, evt_type, treatment, plot_type, kvar, obj])
                        hist = asrootpy(myfile.Get(var)).Clone()
    
                        #if obj == 'TTbar' and plot_type == 'Gen':
                        #    hist = RebinView.rebin(hist, mtt_bins)
                        #    hist.xaxis.range_user = mass_min, mass_max
    
                        mean = hist.GetMean()
                        rms = hist.GetRMS()
    
                        #set_trace()
                        plotter.set_histo_style(hist, color=event_types[evt_type], title=evt_type+' Mean = %.2f, RMS = %.2f' % (mean,rms) )
                        to_draw.append(hist)
    
                    #stack, norm_stack, ratio =  fncts.stack_plots(to_draw)
                    plotter.defaults['watermark'] = ['%s Comparison (13 TeV, 25ns)' % treatment, False]
                    plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
    
                    plotter.save('%s_%s_%s_%s_Comparison' % (plot_type, kvar, obj, treatment))
                    #set_trace()


##############################################################################################

def Gen_Plots(directory, subdir, topology):

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]
    rebin_hist = {'Costh' : {'THad' : (-1., 1., 26), 'TLep' : (-1., 1., 26)},\
                  'Eta' : {'TTbar' : (-2.5, 2.5, 26), 'THad' : (-2.5, 2.5, 26), 'TLep' : (-2.5, 2.5, 26)},\
                  'Mass' : {'TTbar' : (200., 2000., 46), 'THad' : (100., 250., 46)},\
                  'Pt' : {'TTbar' : (0., 1000., 51), 'THad' : (0., 1000., 51), 'TLep' : (0., 1000., 51)}}

        ### Gen plots for thad, tlep, and ttbar objects looking at kinvar dists (mass, costh, ...) based on likelihood value
    for kvar in var_types.keys():
        plotter.set_subdir('/'.join([subdir, 'Gen', kvar]))
        #plotter.set_subdir('/'.join([directory, 'Gen', kvar]))
        for obj in var_types[kvar].keys():
            xlabel = 'Gen %s' % var_types[kvar][obj]
            hname = '/'.join([directory, 'Gen', kvar, obj])

            hist = asrootpy(myfile.Get(hname)).Clone()

            if hist.Integral() == 0:
                continue

            #new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], rebin_hist[kvar][obj][2])
            new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], hist.nbins()/2+1)
            hist = RebinView.rebin(hist, new_bins)
            #hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

            hist_mean = hist.GetMean()
            hist_rms = hist.GetRMS()

            plotter.set_histo_style(hist, xtitle=xlabel, ytitle=defyax)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hist_mean, hist_rms), position='NE')
            box1.Draw()
            plotter.save('Gen_%s_%s' % (obj, kvar))

                ### create hists based on perm category (RIGHT/WRONG, etc...)
            to_draw = []
            for cat in Categories[topology]:
                evt_col=plotter.styles[cat]['fillcolor']
                evt_type = plotter.styles[cat]['name']
                cat_hname = '/'.join([directory, 'Gen', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

                if cat_hist.Integral() == 0:
                    continue
                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(plotter.styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)

            if not to_draw:
                continue
            stack, norm_stack, ratio = fncts.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Gen_%s_%s_Stack' % (obj, kvar))

            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('gen_%s_%s_Stack_Norm' % (obj, kvar))


##############################################################################################


def Reco_Plots(directory, subdir, topology):

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]
    rebin_hist = {'Costh' : {'THad' : (-1., 1., 26), 'TLep' : (-1., 1., 26)},\
                  'Eta' : {'TTbar' : (-2.5, 2.5, 26), 'THad' : (-2.5, 2.5, 26), 'TLep' : (-2.5, 2.5, 26)},\
                  'Mass' : {'TTbar' : (200., 2000., 46), 'THad' : (50., 250., 46)},\
                  'Pt' : {'TTbar' : (0., 1000., 51), 'THad' : (0., 1000., 51), 'TLep' : (0., 1000., 51)}}

        ### Reco plots for thad, tlep, and ttbar objects looking at kinvar dists (mass, costh, ...) based on likelihood value
    for kvar in var_types.keys():

        plotter.set_subdir('/'.join([subdir, 'Reconstruction', kvar]))
        #plotter.set_subdir('/'.join([directory, 'Reconstruction', kvar]))
        for obj in var_types[kvar].keys():
            xlabel = 'Reco %s' % var_types[kvar][obj]
            hname = '/'.join([directory, 'Reconstruction', kvar, obj])

            hist = asrootpy(myfile.Get(hname)).Clone()

            if hist.Integral() == 0:
                continue

            #if kvar == 'Mass':
            #    pre_rebin = [hist.Integral(i+1,i+1) for i in range(hist.nbins())]
            #    #new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], rebin_hist[kvar][obj][2])
            #    post_rebin = [hist.Integral(i+1,i+1) for i in range(hist.nbins())]
            #set_trace()
            new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], hist.nbins()/2+1)
            hist = RebinView.rebin(hist, new_bins)
            #hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

            hist_mean = hist.GetMean()
            hist_rms = hist.GetRMS()

            plotter.set_histo_style(hist, xtitle=xlabel, ytitle=defyax)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hist_mean, hist_rms), position='NE')
            box1.Draw()
            plotter.save('Reco_%s_%s' % (obj, kvar))


                ### create hists based on perm category (RIGHT/WRONG, etc...)

            if topology:
                to_draw = []
                for cat in Categories[topology]:
                    evt_col=plotter.styles[cat]['fillcolor']
                    evt_type = plotter.styles[cat]['name']
                    cat_hname = '/'.join([directory, 'Reconstruction', cat, kvar, obj])
                    cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                    cat_hist = RebinView.rebin(cat_hist, new_bins)
                    #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                    if cat_hist.Integral() == 0:
                        continue
    
                    plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                    cat_hist.SetFillStyle(plotter.styles[cat]['fillstyle'])
                    plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                    to_draw.append(cat_hist)
    
                if not to_draw:
                    continue
                stack, norm_stack, ratio = fncts.stack_plots(to_draw)
                plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                box1.Draw()
                plotter.save('Reco_%s_%s_Stack' % (obj, kvar))
    
                plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                box1.Draw()
                plotter.save('Reco_%s_%s_Stack_Norm' % (obj, kvar))
    #set_trace()
    #            if kvar == 'Mass' and obj == 'THad':
    #                plotter.set_histo_style(reco_hist, color=col)
    #                plotter.plot(reco_hist, drawstyle='hist')
    #                r = reco_hist.Fit("gaus", "S")
    #                Mean = r.Parameter(1)
    #                sigma = r.Parameter(2)
    #                reco_hist.SetTitle(legends+' Fit Mean = %.0f, #sigma = %.0f' % (Mean, sigma))

    #plot_dir(directory+'/'+plots) 



##############################################################################################


def Resolution_Plots(directory, subdir, topology):

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]
    rebin_hist = {'Costh' : {'TTbar' : (-2., 2., 26), 'THad' : (-2., 2., 26), 'TLep' : (-2., 2., 26)},\
                  'Eta' : {'TTbar' : (-5.0, 50, 26), 'THad' : (-5.0, 5.0, 26), 'TLep' : (-5.0, 5.0, 26)},\
                  'Mass' : {'TTbar' : (-1000., 2000., 46), 'THad' : (-1000., 500., 46), 'TLep' : (-500., 500., 46)},\
                  'Pt' : {'TTbar' : (-1000., 500., 51), 'THad' : (-500., 500., 51), 'TLep' : (-1000., 1000., 51)}}

        ### Resolution plots for thad, tlep, and ttbar objects looking at kinvar dists (mass, costh, ...) based on likelihood value
    for kvar in var_types.keys():
        plotter.set_subdir('/'.join([subdir, 'Resolution', kvar]))
        #plotter.set_subdir('/'.join([directory, 'Resolution', kvar]))
        for obj in var_types[kvar].keys():
            xlabel = '%s Resolution (Gen-Reco)' % var_types[kvar][obj]
            hname = '/'.join([directory, 'Resolution', kvar, obj])

            hist = asrootpy(myfile.Get(hname)).Clone()

            if hist.Integral() == 0:
                continue

            #new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], rebin_hist[kvar][obj][2])
            new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], hist.nbins()/2+1)
            hist = RebinView.rebin(hist, new_bins)
            #hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

            hist_mean = hist.GetMean()
            hist_rms = hist.GetRMS()

            plotter.set_histo_style(hist, xtitle=xlabel, ytitle=defyax)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hist_mean, hist_rms), position='NE')
            box1.Draw()
            plotter.save('Reso_%s_%s' % (obj, kvar))


                    ### create hists based on perm category (RIGHT/WRONG, etc...)
            if topology:

                to_draw = []
                for cat in Categories[topology]:
                    evt_col=plotter.styles[cat]['fillcolor']
                    evt_type = plotter.styles[cat]['name']
                    cat_hname = '/'.join([directory, 'Resolution', cat, kvar, obj])
                    cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                    if cat_hist.Integral() == 0:
                        continue
                    cat_hist = RebinView.rebin(cat_hist, new_bins)
                    #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                    plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                    cat_hist.SetFillStyle(plotter.styles[cat]['fillstyle'])
                    plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                    to_draw.append(cat_hist)
    
                if not to_draw:
                    continue
                stack, norm_stack, ratio = fncts.stack_plots(to_draw)
                plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                box1.Draw()
                plotter.save('Reso_%s_%s_Stack' % (obj, kvar))
    
                plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                box1.Draw()
                plotter.save('Reso_%s_%s_Stack_Norm' % (obj, kvar))

    
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

            box1 = plotter.make_text_box(evt_type, position='NE')
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


def Discriminant_Plots(directory, subdir, topology):

        ### Discriminant plots for Mass, NS, and Combined discs looking at evt_type (RIGHT, WRONG, ...)

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]
    plotter.set_subdir('/'.join([subdir, 'Discr']))
    #plotter.set_subdir('/'.join([directory, 'Discr']))
    for disc in DiscPlots:
        xlabel = DiscPlots[disc][0]
        hname = '/'.join([directory, 'Discr', disc])

        hist = asrootpy(myfile.Get(hname)).Clone()
        if hist.Integral() == 0:
            continue
        hist_mean = hist.GetMean()
        hist_rms = hist.GetRMS()

        plotter.set_histo_style(hist, xtitle=xlabel+' 3 jets', ytitle=defyax)
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hist_mean, hist_rms), position='NE')
        box1.Draw()
        plotter.save(disc)

            ### create hists based on perm category (RIGHT/WRONG, etc...)
        if topology:

            to_draw = []
            for cat in Categories[topology]:
                cat_hname = '/'.join([directory, 'Discr', cat, disc])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
                evt_col=plotter.styles[cat]['fillcolor']
                evt_type = plotter.styles[cat]['name']
    
                if cat_hist.Integral() == 0:
                    continue
                cat_hist_mean = cat_hist.GetMean()
                cat_hist_rms = cat_hist.GetRMS()
    
                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel+' 3 jets', ytitle=defyax)
                cat_hist.SetFillStyle(plotter.styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)
    
            if not to_draw:
                continue
            stack, norm_stack, ratio = fncts.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save(disc+'_Stack')
    
            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save(disc+'_Stack_Norm')
    #set_trace()


#    for disc in DiscPlots:
#
#        hname = '3J_bp_%s' % disc    
#        xlabel = DiscPlots[disc][0]
#        xaxes.append(xlabel)
#
#        if disc == 'Totaldisc':
#            Discriminant_Plots_2D(topology, disc)
#            #set_trace()
#
#        sdir_base = args.plot+'/'+topology+'/'+disc
#
#        draw_comb_discs = []
#        comb_disc_hist = 0
#
    #    for i in range(len(Categories[topology])):
    #        evt_type = Categories[topology][i][0]
    #        evt_col = Categories[topology][i][1]
    #        evt_fill = Categories[topology][i][2]

    #        plotter.set_subdir(sdir_base+'/'+evt_type)
    #        var = directory+'/'+topology+'/'+evt_type+'/'+hname
    #        #print var

    #        hist = asrootpy(myfile.Get(var)).Clone()
    #        #print var, hist.Integral()

    #        if hist.Integral() == 0:
    #            continue
    #        comb_disc_hist += hist

    #        hist_mean = hist.GetMean()
    #        hist_rms = hist.GetRMS()
    #        plotter.set_histo_style(hist, color=evt_col, title=evt_type, xtitle=xlabel+' 3 jets', ytitle=defyax)
    #        hist.SetFillStyle(evt_fill)
    #        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    #        draw_comb_discs.append(hist)

    #        box1 = plotter.make_text_box(decay+'\nMean=%.2f\nRMS=%.2f' % (hist_mean, hist_rms), position='NE')
    #        box1.Draw()

    #        hist_type_name = '%s_%s_%s_%s' % (args.sample, hist_save_name, disc, evt_type)
    #        plotter.save(hist_type_name)

    #    if comb_disc_hist == 0:
    #        continue

    #    if topology == 'MERGED':
    #        plotter.set_histo_style(comb_disc_hist, color='red', title=topology)
    #        comb_disc_hist.SetFillStyle(3345)
    #    if topology == 'LOST':
    #        plotter.set_histo_style(comb_disc_hist, color='blue', title=topology)
    #        comb_disc_hist.SetFillStyle(3354)
    #    plotter.plot(comb_disc_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
    #    to_draw.append(comb_disc_hist)

    #    comb_disc_mean = comb_disc_hist.GetMean()
    #    comb_disc_rms = comb_disc_hist.GetRMS()

    #    stack, norm_stack, ratio = fncts.stack_plots(draw_comb_discs)
    #    plotter.set_subdir(sdir_base)

    #    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
    #    box2 = plotter.make_text_box(decay+'\nMean=%.2f\nRMS=%.2f' % (comb_disc_mean, comb_disc_rms), position='NE')
    #    box2.Draw()

    #    comb_disc_hist_name = '%s_%s_%s' % (args.sample, hist_save_name, disc)
    #    plotter.save(comb_disc_hist_name+'_Stack')

    #    plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
    #    box2.Draw()
    #    plotter.save(comb_disc_hist_name+'_Stack_Norm')

    #plot_dir(args.plot+'/'+topology) 

    #return to_draw, xaxes
    

#####################################################################################################

#if args.perm == 'Event':
#
#    if args.plot == 'Final_Reco':
#        Final_Reco_Plots()
#
#    else:
#        if args.plot == 'All':
#            Event_Plots('Discriminant')
#            Event_Plots('Gen')
#            Event_Plots('Reconstruction')
#            Event_Plots('Resolution')
#        else:
#            Event_Plots(args.plot)

if args.perm == 'Final_Class':
    if args.plot == 'All':
        Final_Reco_Plots('Discriminant')
        Final_Reco_Plots('Gen')
        Final_Reco_Plots('Reconstruction')
        Final_Reco_Plots('Resolution')
    else:
        Final_Reco_Plots(args.plot)

if args.perm == 'Best_Perm':
    if args.plot == 'All':
        Event_Plots('Discriminant')
        Event_Plots('Gen')
        Event_Plots('Reconstruction')
        Event_Plots('Resolution')
    else:
        Event_Plots(args.plot)

if args.perm == 'Matched':
    if args.plot == 'All':
        Matched_Perm_Plots('Discriminant')
        Matched_Perm_Plots('Gen')
        Matched_Perm_Plots('Reconstruction')
        Matched_Perm_Plots('Resolution')
    else:
        Matched_Perm_Plots(args.plot)

