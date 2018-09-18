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


out_analyzer = 'ttbar_post_alpha_reco'

parser = argparse.ArgumentParser(description='Create plots using files from ttbar_reco_3J.')

jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']

parser.add_argument('--analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('--sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, AtoTT_M...).')
parser.add_argument('--perm', help='Choose type of perm to take information from (Event->best perm, Matched->matched perm).')
parser.add_argument('--plot', help='Choose type of plots to generate (Reconstruction, Resolution, Partons_Acceptances, All).')
parser.add_argument('--everything', default=False, help='Run Gen, Reco, Resolution, and Discr plots for matched and event perms.')
args = parser.parse_args()


##### check analysis type
if not (args.analysis == "Test" or args.analysis == "Full"):
    print "You chose %s as your analysis type.\nYou must choose Full or Test!" % args.analysis
    sys.exit()

##### check perm type
if not (args.perm == "Event" or args.perm == "Matched"):
    print "You chose %s as your permutation type.\nYou must choose Event or Matched!" % args.perm
    sys.exit()


analyzer = 'ttbar_post_alpha_reco'


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
            'sample' : styles['ttJetsM0'] if args.sample=='ttJets' else styles[args.sample]
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
            'sample' : styles['ttJetsM0'] if args.sample=='ttJets' else styles[args.sample]
        }
    )


allowed_plot_types = ['Reconstruction', 'Resolution', 'Partons_Acceptances', 'All']

if not args.plot in allowed_plot_types:
    print "You chose %s as your plot type.\nYou must choose from the help list!" % args.plot
    sys.exit()



myfile = root_open(results_file, 'read')
normfile = views.NormalizeView( root_open(results_file, 'read') )
lumifile = open('%s/inputs/%s/%s.lumi' % (project, jobid, args.sample), 'read')


#set_trace()

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

plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]

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

#Perm_Categories = ['NO_MP', 'CORRECT_WJET_CORRECT_Bs', 'CORRECT_WJET_SWAPPED_Bs', 'CORRECT_WJET_CORRECT_BHAD', 'CORRECT_WJET_CORRECT_BLEP', 'CORRECT_WJET_WRONG_Bs',
#                    'WRONG_WJET_CORRECT_Bs', 'WRONG_WJET_SWAPPED_Bs', 'WRONG_WJET_CORRECT_BHAD', 'WRONG_WJET_CORRECT_BLEP', 'WRONG_WJET_WRONG_Bs'
#                ]
Correct_WJet_Categories = ['CORRECT_WJET_CORRECT_Bs', 'CORRECT_WJET_SWAPPED_Bs', 'CORRECT_WJET_CORRECT_BHAD', 'CORRECT_WJET_CORRECT_BLEP', 'CORRECT_WJET_WRONG_Bs']
Wrong_WJet_Categories = ['WRONG_WJET_CORRECT_Bs', 'WRONG_WJET_SWAPPED_Bs', 'WRONG_WJET_CORRECT_BHAD', 'WRONG_WJET_CORRECT_BLEP', 'WRONG_WJET_WRONG_Bs']


#Categories = {'MERGED' : ['RIGHT', 'MERGED_SWAP', 'MERGED', 'WRONG'],\
#              'LOST' : ['RIGHT', 'LOST_SWAP', 'LOST', 'WRONG']}
#Categories = {'MERGED' : [('RIGHT', 'red', 3345), ('MERGED_SWAP', 'orange', 0), ('MERGED', 'green', 3006), ('WRONG', 'blue', 3354)],\
#              'LOST' : [('RIGHT', 'red', 3345), ('LOST_SWAP', 'orange', 0), ('LOST', 'green', 3006), ('WRONG', 'blue', 3354)]}
Objects = {'THad' : 't_{h}', 'TLep' : 't_{l}', 'TTbar' : 't#bar{t}'}
CUT = {'LCut' : ['#lambda_{comb} < 2', 'black'], 'GCut' : ['#lambda_{comb} #geq 2', 'red']}
kinvar = {'Mass' : ' Mass [GeV]', 'Pt' : ' p_{T}', 'Eta' : ' #eta', 'Costh' : ' cos(#theta^{*})'}
DiscPlots = {'Massdisc' : ['#lambda_{mass}'], 'NSdisc' : ['#lambda_{NS}'], 'Totaldisc' : ['#lambda_{comb}']}

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
def Reconstruction_Plots(directory, subdir):

#    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]


############# plots for reconstructed vars

    reco_hists = {
        'Mass' : {
                    'THad' : {  'Uncorrected' : ('', 50., 200., 'k', 'Uncorrected', 'solid'),
                                'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 50., 200., 'r', '1d #alpha_{E} All', 'solid'), 'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 50., 200., 'g', '2d #alpha_{E} All', 'solid'),
                                'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 50., 200., 'b', '1d #alpha_{P} All', 'solid'), 'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 50., 200., 'm', '2d #alpha_{P} All', 'solid'),
                                'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 50., 200., 'r', '1d #alpha_{E} Mtt', 'dashed'), 'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 50., 200., 'g', '2d #alpha_{E} Mtt', 'dashed'),
                                'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 50., 200., 'b', '1d #alpha_{P} Mtt', 'dashed'), 'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 50., 200., 'm', '2d #alpha_{P} Mtt', 'dashed')
                    },
                    'TTbar' : { 'Uncorrected' : ('', 200., 2000., 'k', 'Uncorrected', 'solid'),
                                'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 200., 2000., 'r', '1d #alpha_{E} All', 'solid'), 'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 200., 2000., 'g', '2d #alpha_{E} All', 'solid'),
                                'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 200., 2000., 'b', '1d #alpha_{P} All', 'solid'), 'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 200., 2000., 'm', '2d #alpha_{P} All', 'solid'),
                                'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 200., 2000., 'r', '1d #alpha_{E} Mtt', 'dashed'), 'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 200., 2000., 'g', '2d #alpha_{E} Mtt', 'dashed'),
                                'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 200., 2000., 'b', '1d #alpha_{P} Mtt', 'dashed'), 'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 200., 2000., 'm', '2d #alpha_{P} Mtt', 'dashed')
                    },
                    'Reco_vs_Gen_TTbar' : { 'Uncorrected' : ('', 'Gen M(t#bar{t}) (GeV)', 'Uncorrected Reco M(t#bar{t}) (GeV)'),
                                            'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 'Gen M(t#bar{t}) (GeV)', '1d #alpha_{E} All Reco M(t#bar{t}) (GeV)'),
                                            'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 'Gen M(t#bar{t}) (GeV)', '2d #alpha_{E} All Reco M(t#bar{t}) (GeV)'),
                                            'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 'Gen M(t#bar{t}) (GeV)', '1d #alpha_{P} All Reco M(t#bar{t}) (GeV)'),
                                            'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 'Gen M(t#bar{t}) (GeV)', '2d #alpha_{P} All Reco M(t#bar{t}) (GeV)'),
                                            'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 'Gen M(t#bar{t}) (GeV)', '1d #alpha_{E} Mtt Reco M(t#bar{t}) (GeV)'),
                                            'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 'Gen M(t#bar{t}) (GeV)', '2d #alpha_{E} Mtt Reco M(t#bar{t}) (GeV)'),
                                            'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 'Gen M(t#bar{t}) (GeV)', '1d #alpha_{P} Mtt Reco M(t#bar{t}) (GeV)'),
                                            'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 'Gen M(t#bar{t}) (GeV)', '2d #alpha_{P} Mtt Reco M(t#bar{t}) (GeV)')
                    },
         },
        'Costh' : {
                    'THad' : { 'Uncorrected' : ('', -1., 1., 'k', 'Uncorrected', 'solid'),
                               'Corrected_1D_E_All' : ('_Corrected_1D_E_All', -1., 1., 'r', '1d #alpha_{E} All', 'solid'), 'Corrected_2D_E_All' : ('_Corrected_2D_E_All', -1., 1., 'g', '2d #alpha_{E} All', 'solid'),
                               'Corrected_1D_P_All' : ('_Corrected_1D_P_All', -1., 1., 'b', '1d #alpha_{P} All', 'solid'), 'Corrected_2D_P_All' : ('_Corrected_2D_P_All', -1., 1., 'm', '2d #alpha_{P} All', 'solid'),
                               'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', -1., 1., 'r', '1d #alpha_{E} Mtt', 'dashed'), 'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', -1., 1., 'g', '2d #alpha_{E} Mtt', 'dashed'),
                               'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', -1., 1., 'b', '1d #alpha_{P} Mtt', 'dashed'), 'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', -1., 1., 'm', '2d #alpha_{P} Mtt', 'dashed')
                    },
                    #'THad_Labframe' : { 'Uncorrected' : ('', -1., 1., 'k', 'Uncorrected', 'solid'),
                    #                    'Corrected_1D_E' : ('_Corrected_1D_E', -1., 1., 'r', 'Corrected 1D #alpha_{E}'), 'Corrected_2D_E' : ('_Corrected_2D_E', -1., 1., 'g', 'Corrected 2D #alpha_{E}'),
                    #                    'Corrected_1D_P' : ('_Corrected_1D_P', -1., 1., 'b', 'Corrected 1D #alpha_{P}'), 'Corrected_2D_P' : ('_Corrected_2D_P', -1., 1., 'm', 'Corrected 2D #alpha_{P}')
                    #},
                    'Reco_vs_Gen_THad' : { 'Uncorrected' : ('', 'Gen cos(#theta^{*})', 'Reco cos(#theta^{*})'),
                                           'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 'Gen cos(#theta^{*})', '1d #alpha_{E} All Reco cos(#theta^{*})'),
                                           'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 'Gen cos(#theta^{*})', '2d #alpha_{E} All Reco cos(#theta^{*})'), 
                                           'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 'Gen cos(#theta^{*})', '1d #alpha_{P} All Reco cos(#theta^{*})'),
                                           'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 'Gen cos(#theta^{*})', '2d #alpha_{P} All Reco cos(#theta^{*})'),
                                           'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 'Gen cos(#theta^{*})', '1d #alpha_{E} Mtt Reco cos(#theta^{*})'),
                                           'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 'Gen cos(#theta^{*})', '2d #alpha_{E} Mtt Reco cos(#theta^{*})'), 
                                           'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 'Gen cos(#theta^{*})', '1d #alpha_{P} Mtt Reco cos(#theta^{*})'),
                                           'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 'Gen cos(#theta^{*})', '2d #alpha_{P} Mtt Reco cos(#theta^{*})') 
                    },
        }
    }

    #set_trace()

    for kvar in reco_hists.keys():
        plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Reconstruction', kvar]))
   
        for obj in reco_hists[kvar].keys():

            if var_types[kvar].has_key(obj):
                xlabel = 'Reco %s' % var_types[kvar][obj]
            elif obj == 'THad_Labframe':
                xlabel = 'Reco cos(#theta)(t_{h}) LF'
            else:
                xlabel = ''

            to_draw = []
            for corr_type in reco_hists[kvar][obj].keys():
                hvar_extension = reco_hists[kvar][obj][corr_type][0]
                hname = '/'.join([directory, 'Post_Alpha_Correction', 'Reconstruction', kvar, obj+hvar_extension])

                hist = asrootpy(myfile.Get(hname)).Clone()

                if hist.Integral() == 0:
                    continue
  
                if hist.DIM == 2:
                    plotter.set_histo_style(hist, xtitle=reco_hists[kvar][obj][corr_type][1], ytitle=reco_hists[kvar][obj][corr_type][2])
                    plotter.plot(hist, drawstyle='colz')

                    plotter.save('%s_%s_%s' % (obj, kvar, corr_type))
                    #set_trace()
                    continue
                #set_trace()

                    ## set hist range, get hist mean and rms 
                hvar_xmin = reco_hists[kvar][obj][corr_type][1] #xmin value
                hvar_xmax = reco_hists[kvar][obj][corr_type][2] #xmax value
                hvar_col = reco_hists[kvar][obj][corr_type][3] #color
                hvar_title = reco_hists[kvar][obj][corr_type][4] #title
                hvar_lstyle = reco_hists[kvar][obj][corr_type][5] #linestyle

                hist.xaxis.range_user = hvar_xmin, hvar_xmax
                hmean = hist.GetMean()
                hrms = hist.GetRMS()

                    ## set hist style, plot, and save
                plotter.set_histo_style(hist, xtitle='%s %s' % (hvar_title, xlabel), ytitle=defyax)
                plotter.plot(hist, drawstyle='hist')
                box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
                box1.Draw()
                plotter.save('Reco_%s_%s_%s' % (kvar, obj, corr_type))

                    ## set hist style to be compared with other type
                plotter.set_histo_style(hist, color=hvar_col, title='%s Mean=%.2f, RMS=%.2f' % (hvar_title, hmean, hrms), linestyle=hvar_lstyle )
                to_draw.append(hist)


                ## compare uncorreced/corrected hists
            if not to_draw:
                continue
            plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            plotter.save('Reco_%s_%s_Comparison' % (kvar, obj) )
            #set_trace()



def Reso_Plots(directory, subdir):
##### plots for resolution

    reso_hists = {
        'Mass' : {
                    'THad' : { 'Uncorrected' : ('', -40., 160., 'k', 'Uncorrected', 'solid'),
                               'Corrected_1D_E_All' : ('_Corrected_1D_E_All', -40., 160., 'r', '1d #alpha_{E} All', 'solid'), 'Corrected_2D_E_All' : ('_Corrected_2D_E_All', -40., 160., 'g', '2d #alpha_{E} All', 'solid'),
                               'Corrected_1D_P_All' : ('_Corrected_1D_P_All', -40., 160., 'b', '1d #alpha_{P} All', 'solid'), 'Corrected_2D_P_All' : ('_Corrected_2D_P_All', -40., 160., 'm', '2d #alpha_{P} All', 'solid'),
                               'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', -40., 160., 'r', '1d #alpha_{E} Mtt', 'dashed'), 'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', -40., 160., 'g', '2d #alpha_{E} Mtt', 'dashed'),
                               'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', -40., 160., 'b', '1d #alpha_{P} Mtt', 'dashed'), 'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', -40., 160., 'm', '2d #alpha_{P} Mtt', 'dashed')
                    },
                    'TTbar' : { 'Uncorrected' : ('', -300., 300., 'k', 'Uncorrected', 'solid'),
                                'Corrected_1D_E_All' : ('_Corrected_1D_E_All', -300., 300., 'r', '1d #alpha_{E} All', 'solid'), 'Corrected_2D_E_All' : ('_Corrected_2D_E_All', -300., 300., 'g', '2d #alpha_{E} All', 'solid'),
                                'Corrected_1D_P_All' : ('_Corrected_1D_P_All', -300., 300., 'b', '1d #alpha_{P} All', 'solid'), 'Corrected_2D_P_All' : ('_Corrected_2D_P_All', -300., 300., 'm', '2d #alpha_{P} All', 'solid'),
                                'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', -300., 300., 'r', '1d #alpha_{E} Mtt', 'dashed'), 'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', -300., 300., 'g', '2d #alpha_{E} Mtt', 'dashed'),
                                'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', -300., 300., 'b', '1d #alpha_{P} Mtt', 'dashed'), 'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', -300., 300., 'm', '2d #alpha_{P} Mtt', 'dashed')
                    },
                    'TLep' : { 'Uncorrected' : ('', -200., 200., 'k', 'Uncorrected', 'solid'),
                               #'Corrected_1D_E_All' : ('_Corrected_1D_E_All', -200., 200., 'r', '1d #alpha_{E} All', 'solid'), 'Corrected_2D_E_All' : ('_Corrected_2D_E_All', -200., 200., 'g', '2d #alpha_{E} All', 'solid'),
                               #'Corrected_1D_P_All' : ('_Corrected_1D_P_All', -200., 200., 'b', '1d #alpha_{P} All', 'solid'), 'Corrected_2D_P_All' : ('_Corrected_2D_P_All', -200., 200., 'm', '2d #alpha_{P} All', 'solid'),
                               #'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', -200., 200., 'r', '1d #alpha_{E} Mtt', 'dashed'), 'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', -200., 200., 'g', '2d #alpha_{E} Mtt', 'dashed'),
                               #'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', -200., 200., 'b', '1d #alpha_{P} Mtt', 'dashed'), 'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', -200., 200., 'm', '2d #alpha_{P} Mtt', 'dashed')
                    },
                    'Frac_THad' : { 'Uncorrected' : ('', -1., 1., 'k', 'Uncorrected', 'solid'),
                                    'Corrected_1D_E_All' : ('_Corrected_1D_E_All', -1.0, 1.0, 'r', '1d #alpha_{E} All', 'solid'), 'Corrected_2D_E_All' : ('_Corrected_2D_E_All', -1.0, 1.0, 'g', '2d #alpha_{E} All', 'solid'),
                                    'Corrected_1D_P_All' : ('_Corrected_1D_P_All', -1.0, 1.0, 'b', '1d #alpha_{P} All', 'solid'), 'Corrected_2D_P_All' : ('_Corrected_2D_P_All', -1.0, 1.0, 'm', '2d #alpha_{P} All', 'solid'),
                                    'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', -1.0, 1.0, 'r', '1d #alpha_{E} Mtt', 'dashed'), 'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', -1.0, 1.0, 'g', '2d #alpha_{E} Mtt', 'dashed'),
                                    'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', -1.0, 1.0, 'b', '1d #alpha_{P} Mtt', 'dashed'), 'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', -1.0, 1.0, 'm', '2d #alpha_{P} Mtt', 'dashed')
                    },
                    'Frac_TTbar' : { 'Uncorrected' : ('', -1., 1., 'k', 'Uncorrected', 'solid'),
                                     'Corrected_1D_E_All' : ('_Corrected_1D_E_All', -1., 1., 'r', '1d #alpha_{E} All', 'solid'), 'Corrected_2D_E_All' : ('_Corrected_2D_E_All', -1., 1., 'g', '2d #alpha_{E} All', 'solid'),
                                     'Corrected_1D_P_All' : ('_Corrected_1D_P_All', -1., 1., 'b', '1d #alpha_{P} All', 'solid'), 'Corrected_2D_P_All' : ('_Corrected_2D_P_All', -1., 1., 'm', '2d #alpha_{P} All', 'solid'),
                                     'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', -1., 1., 'r', '1d #alpha_{E} Mtt', 'dashed'), 'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', -1., 1., 'g', '2d #alpha_{E} Mtt', 'dashed'),
                                     'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', -1., 1., 'b', '1d #alpha_{P} Mtt', 'dashed'), 'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', -1., 1., 'm', '2d #alpha_{P} Mtt', 'dashed')
                    },
                    'Reso_MTTbar_vs_Gen_MTTbar' : { 'Uncorrected' : ('', 'Gen M(t#bar{t}) [GeV]', 'Uncorrected M(t#bar{t}) Resolution [GeV]'),
                                                    'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 'Gen M(t#bar{t}) [GeV]', '1d #alpha_{E} All M(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 'Gen M(t#bar{t}) [GeV]', '2d #alpha_{E} All M(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 'Gen M(t#bar{t}) [GeV]', '1d #alpha_{P} All M(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 'Gen M(t#bar{t}) [GeV]', '2d #alpha_{P} All M(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 'Gen M(t#bar{t}) [GeV]', '1d #alpha_{E} Mtt M(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 'Gen M(t#bar{t}) [GeV]', '2d #alpha_{E} Mtt M(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 'Gen M(t#bar{t}) [GeV]', '1d #alpha_{P} Mtt M(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 'Gen M(t#bar{t}) [GeV]', '2d #alpha_{P} Mtt M(t#bar{t}) Res. [GeV]')
                    },
                    'Reso_MTHad_vs_Gen_MTTbar' : { 'Uncorrected' : ('', 'Gen M(t#bar{t}) [GeV]', 'Uncorrected M(t_{h}) Resolution [GeV]'),
                                                    'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 'Gen M(t#bar{t}) [GeV]', '1d #alpha_{E} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 'Gen M(t#bar{t}) [GeV]', '2d #alpha_{E} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 'Gen M(t#bar{t}) [GeV]', '1d #alpha_{P} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 'Gen M(t#bar{t}) [GeV]', '2d #alpha_{P} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 'Gen M(t#bar{t}) [GeV]', '1d #alpha_{E} Mtt M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 'Gen M(t#bar{t}) [GeV]', '2d #alpha_{E} Mtt M(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 'Gen M(t#bar{t}) [GeV]', '1d #alpha_{P} Mtt M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 'Gen M(t#bar{t}) [GeV]', '2d #alpha_{P} Mtt M(t_{h}) Res. [GeV]')
                    },
                    'Reso_MTHad_vs_Gen_THadPt' : { 'Uncorrected' : ('', 'Gen p_{T}(t_{h}) [GeV]', 'Uncorrected M(t_{h}) Resolution [GeV]'),
                                                    'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{E} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{E} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{P} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{P} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{E} Mtt M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{E} Mtt M(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{P} Mtt M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{P} Mtt M(t_{h}) Res. [GeV]')
                    },
                    'Reso_MTHad_vs_Uncorrected_MTHad' : { 'Uncorrected' : ('', '173.1/Uncorrected M(t_{h}) [GeV]', 'Uncorrected M(t_{h}) Resolution [GeV]'),
                                                    'Corrected_1D_E_All' : ('_Corrected_1D_E_All', '173.1/Uncorrected M(t_{h}) [GeV]', '1d #alpha_{E} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_E_All' : ('_Corrected_2D_E_All', '173.1/Uncorrected M(t_{h}) [GeV]', '2d #alpha_{E} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_P_All' : ('_Corrected_1D_P_All', '173.1/Uncorrected M(t_{h}) [GeV]', '1d #alpha_{P} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_P_All' : ('_Corrected_2D_P_All', '173.1/Uncorrected M(t_{h}) [GeV]', '2d #alpha_{P} All M(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', '173.1/Uncorrected M(t_{h}) [GeV]', '1d #alpha_{E} Mtt M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', '173.1/Uncorrected M(t_{h}) [GeV]', '2d #alpha_{E} Mtt M(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', '173.1/Uncorrected M(t_{h}) [GeV]', '1d #alpha_{P} Mtt M(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', '173.1/Uncorrected M(t_{h}) [GeV]', '2d #alpha_{P} Mtt M(t_{h}) Res. [GeV]')
                    },
                    'Reso_THadPt_vs_Gen_THadPt' : { 'Uncorrected' : ('', 'Gen p_{T}(t_{h}) [GeV]', 'Uncorrected p_{T}(t_{h}) Resolution [GeV]'),
                                                    'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{E} All p_{T}(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{E} All p_{T}(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{P} All p_{T}(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{P} All p_{T}(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{E} Mtt p_{T}(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{E} Mtt p_{T}(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{P} Mtt p_{T}(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{P} Mtt p_{T}(t_{h}) Res. [GeV]')
                    },
                    'Reso_THadEta_vs_Gen_THadEta' : { 'Uncorrected' : ('', 'Gen #eta(t_{h}) [GeV]', 'Uncorrected #eta(t_{h}) Resolution [GeV]'),
                                                    'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 'Gen #eta(t_{h}) [GeV]', '1d #alpha_{E} All #eta(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 'Gen #eta(t_{h}) [GeV]', '2d #alpha_{E} All #eta(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 'Gen #eta(t_{h}) [GeV]', '1d #alpha_{P} All #eta(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 'Gen #eta(t_{h}) [GeV]', '2d #alpha_{P} All #eta(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 'Gen #eta(t_{h}) [GeV]', '1d #alpha_{E} Mtt #eta(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 'Gen #eta(t_{h}) [GeV]', '2d #alpha_{E} Mtt #eta(t_{h}) Res. [GeV]'),
                                                    'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 'Gen #eta(t_{h}) [GeV]', '1d #alpha_{P} Mtt #eta(t_{h}) Res. [GeV]'),
                                                    'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 'Gen #eta(t_{h}) [GeV]', '2d #alpha_{P} Mtt #eta(t_{h}) Res. [GeV]')
                    },
                    'Reso_EtaTTbar_vs_Gen_EtaTTbar' : { 'Uncorrected' : ('', 'Gen #eta(t#bar{t}) [GeV]', 'Uncorrected #eta(t#bar{t}) Resolution [GeV]'),
                                                    'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 'Gen #eta(t#bar{t}) [GeV]', '1d #alpha_{E} All #eta(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 'Gen #eta(t#bar{t}) [GeV]', '2d #alpha_{E} All #eta(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 'Gen #eta(t#bar{t}) [GeV]', '1d #alpha_{P} All #eta(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 'Gen #eta(t#bar{t}) [GeV]', '2d #alpha_{P} All #eta(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 'Gen #eta(t#bar{t}) [GeV]', '1d #alpha_{E} Mtt #eta(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 'Gen #eta(t#bar{t}) [GeV]', '2d #alpha_{E} Mtt #eta(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 'Gen #eta(t#bar{t}) [GeV]', '1d #alpha_{P} Mtt #eta(t#bar{t}) Res. [GeV]'),
                                                    'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 'Gen #eta(t#bar{t}) [GeV]', '2d #alpha_{P} Mtt #eta(t#bar{t}) Res. [GeV]')
                    },
                    'Frac_TTbar_vs_Gen_THadPt' : { 'Uncorrected' : ('', 'Gen p_{T}(t_{h}) [GeV]', 'Uncorrected M(t#bar{t}) Frac Resolution'),
                                                   'Corrected_1D_E_All' : ('_Corrected_1D_E_All', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{E} All M(t#bar{t}) Frac Res.'),
                                                   'Corrected_2D_E_All' : ('_Corrected_2D_E_All', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{E} All M(t#bar{t}) Frac Res.'),
                                                   'Corrected_1D_P_All' : ('_Corrected_1D_P_All', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{P} All M(t#bar{t}) Frac Res.'),
                                                   'Corrected_2D_P_All' : ('_Corrected_2D_P_All', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{P} All M(t#bar{t}) Frac Res.'),
                                                   'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{E} Mtt M(t#bar{t}) Frac Res.'),
                                                   'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{E} Mtt M(t#bar{t}) Frac Res.'),
                                                   'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '1d #alpha_{P} Mtt M(t#bar{t}) Frac Res.'),
                                                   'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', 'Gen p_{T}(t_{h}) [GeV]', '2d #alpha_{P} Mtt M(t#bar{t}) Frac Res.'),
                    },
         },
        'Costh' : {
                    'THad' : { 'Uncorrected' : ('', -2., 2., 'k', 'Uncorrected', 'solid'),
                               'Corrected_1D_E_All' : ('_Corrected_1D_E_All', -2., 2., 'r', '1d #alpha_{E} All', 'solid'), 'Corrected_2D_E_All' : ('_Corrected_2D_E_All', -2., 2., 'g', '2d #alpha_{E} All', 'solid'),
                               'Corrected_1D_P_All' : ('_Corrected_1D_P_All', -2., 2., 'b', '1d #alpha_{P} All', 'solid'), 'Corrected_2D_P_All' : ('_Corrected_2D_P_All', -2., 2., 'm', '2d #alpha_{P} All', 'solid'),
                               'Corrected_1D_E_Mtt' : ('_Corrected_1D_E_Mtt', -2., 2., 'r', '1d #alpha_{E} Mtt', 'dashed'), 'Corrected_2D_E_Mtt' : ('_Corrected_2D_E_Mtt', -2., 2., 'g', '2d #alpha_{E} Mtt', 'dashed'),
                               'Corrected_1D_P_Mtt' : ('_Corrected_1D_P_Mtt', -2., 2., 'b', '1d #alpha_{P} Mtt', 'dashed'), 'Corrected_2D_P_Mtt' : ('_Corrected_2D_P_Mtt', -2., 2., 'm', '2d #alpha_{P} Mtt', 'dashed'),
                  }
        }
    }

    for kvar in reso_hists.keys():
        #plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', kvar]))
   
        for obj in reso_hists[kvar].keys():

            if var_types[kvar].has_key(obj):
                xlabel = '%s Resolution' % var_types[kvar][obj]
            elif obj == 'Frac_THad':
                xlabel = 'M(t_{h}) Frac Resolution (G-R/G)'
            elif obj == 'Frac_TTbar':
                xlabel = 'M(t#bar{t}) Frac Resolution (G-R/G)'
            else:
                xlabel = ''

            to_draw_normal = []

            for corr_type in reso_hists[kvar][obj].keys():

                    ### create hists based on perm category (Correct b, wrong b, etc...)
                #print 'kvar = %s\nobj = %s\ncorr_type= %s ' % (kvar, obj, corr_type)
                if (kvar == 'Costh' and obj == 'THad') or (kvar == 'Mass' and obj == 'THad') or (kvar == 'Mass' and obj == 'TTbar'):
                    to_draw = []
                    correct_to_draw = []
                    wrong_to_draw = []
                    for cat in Correct_WJet_Categories+Wrong_WJet_Categories+['NO_MP']:
                        plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', kvar, 'Perm_Categories']))
                        #set_trace()
                        cat_hvar_extension = reso_hists[kvar][obj][corr_type][0]
                        cat_hname = '/'.join([directory, 'Post_Alpha_Correction', cat, 'Resolution', kvar, obj+cat_hvar_extension])
                        cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

                        #cat_hist = RebinView.rebin(cat_hist, new_bins)
                        #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

                        if cat_hist.Integral() == 0:
                            continue

                            ## set hist range, get hist mean and rms 
                        cat_hvar_xmin = reso_hists[kvar][obj][corr_type][1]
                        cat_hvar_xmax = reso_hists[kvar][obj][corr_type][2]
                        cat_hvar_col = reso_hists[kvar][obj][corr_type][3]
                        cat_hvar_title = reso_hists[kvar][obj][corr_type][4]

                        cat_hist.xaxis.range_user = cat_hvar_xmin, cat_hvar_xmax
                        cat_hist_mean = cat_hist.GetMean() 
                        cat_hist_rms = cat_hist.GetRMS() 
    
                        plotter.set_histo_style(cat_hist, color=styles[cat]['fillcolor'], title='%s Mean=%.2f, RMS=%.2f' % (styles[cat]['name'], cat_hist_mean, cat_hist_rms), fillstyle=styles[cat]['fillstyle'])
                        plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                        cat_hist.SetName(cat)
                        to_draw.append(cat_hist)
    

                    if to_draw:
                        stack, norm_stack, ratio = fncts.stack_plots(to_draw)
                        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='%s %s' % (cat_hvar_title, xlabel), ytitle=defyax, drawstyle='hist', x_range=(reso_hists[kvar][obj][corr_type][1], reso_hists[kvar][obj][corr_type][2]))
                        plotter.save('Reso_%s_%s_%s_Stack' % (kvar, obj, corr_type))

                        #plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                        #box1.Draw()
                        #plotter.save('Reso_%s_%s_%s_Stack_Norm' % (kvar, obj, corr_type))

                        correct_hist = sum([i for i in to_draw if i.name == 'CORRECT_WJET_CORRECT_Bs'])
                        plotter.set_histo_style(correct_hist, color='r', title='Correct Mean=%.2f, RMS=%.2f' % (correct_hist.GetMean(), correct_hist.GetRMS()), fillstyle=1001)
                        plotter.plot(correct_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')

                        wrong_hist = sum([i for i in to_draw if i.name != 'CORRECT_WJET_CORRECT_Bs'])
                        plotter.set_histo_style(wrong_hist, color='b', title='Wrong Mean=%.2f, RMS=%.2f' % (wrong_hist.GetMean(), wrong_hist.GetRMS()), fillstyle=1001)
                        plotter.plot(wrong_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')

                        
                        stack, norm_stack, ratio = fncts.stack_plots([correct_hist, wrong_hist])
                        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle='%s %s' % (cat_hvar_title, xlabel), ytitle=defyax, drawstyle='hist', x_range=(reso_hists[kvar][obj][corr_type][1], reso_hists[kvar][obj][corr_type][2]))
                        box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % ((correct_hist+wrong_hist).GetMean(), (correct_hist+wrong_hist).GetRMS()), position='NE')
                        box1.Draw()
                        plotter.save('Reso_%s_%s_%s_Correct_Wrong_Stack' % (kvar, obj, corr_type))
                    #set_trace()



                    #### normal plots
                plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', kvar]))
                hvar_extension = reso_hists[kvar][obj][corr_type][0]
                hname = '/'.join([directory, 'Post_Alpha_Correction', 'Resolution', kvar, obj+hvar_extension])

                hist = asrootpy(myfile.Get(hname)).Clone()

                if hist.Integral() == 0:
                    continue
  
                if hist.DIM == 2:
                    plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', kvar, obj]))

                    if obj ==  'Reso_MTHad_vs_Gen_MTTbar' or obj == 'Reso_MTHad_vs_Uncorrected_MTHad' or obj == 'Reso_MTHad_vs_Gen_THadPt':
                        xbins = np.linspace(hist.GetXaxis().GetBinLowEdge(1), hist.GetXaxis().GetBinUpEdge(hist.GetXaxis().GetNbins()), hist.GetXaxis().GetNbins()+1 ) # ybinning remains unchanged
                        ybins = np.linspace(-100., 98, 67)
                        hist = RebinView.newRebin2D(hist, xbins, ybins)
                        #set_trace()

                    plotter.set_histo_style(hist, xtitle=reso_hists[kvar][obj][corr_type][1], ytitle=reso_hists[kvar][obj][corr_type][2])
                    plotter.plot(hist, drawstyle='colz')

                    mean_x = hist.ProjectionX().GetMean()
                    rms_x = hist.ProjectionX().GetRMS()
                    mean_y = hist.ProjectionY().GetMean()
                    rms_y = hist.ProjectionY().GetRMS()

                    box2 = plotter.make_text_box('X Mean=%.2f, RMS=%.2f\nY Mean=%.2f, RMS=%.2f' % (mean_x, rms_x, mean_y, rms_y), position='NE')
                    box2.Draw()

                    plotter.save('%s_%s_%s' % (obj, kvar, corr_type))

                    colors = []
                    x_rebin = []
                    ybins = np.linspace(hist.GetYaxis().GetBinLowEdge(1), hist.GetYaxis().GetBinUpEdge(hist.GetYaxis().GetNbins()), hist.GetYaxis().GetNbins()+1 ) # ybinning remains unchanged
                    if 'Gen_MTTbar' in obj:
                        x_rebin = [250., 500., 700., 2000.]
                        colors = ['r', 'b', 'k']
                        if 'MTHad' in obj:
                            ybins = np.linspace(-100., 98, 67)
                            #set_trace()
                    elif 'Uncorrected_MTHad' in obj:
                        x_rebin = [0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
                        colors = ['r', 'b', 'k', 'g', 'c', 'm', '#00FFFF', '#C76114']
                        ybins = np.linspace(-100., 98, 67)
                        #set_trace()
                    elif 'Gen_THadPt' in obj:
                        x_rebin = [0., 100., 200., 1000.]
                        colors = ['r', 'b', 'k']
                    else: # eta bins
                        x_rebin = [-3., -2.4, -1.2, 0., 1.2, 2.4, 3.]
                        colors = ['r', 'b', 'k', 'g', 'c', 'm']
                    #x_rebin = mtt_bins if 'MTTbar' in obj else pt_bins
                    hist_rebinned = RebinView.newRebin2D(hist, x_rebin, ybins)

                    yprojections = []
                    yproj_norms = []
                    i = 0
                    for xbin in range(1, hist_rebinned.GetNbinsX() + 1):
                        #set_trace()
                        #if hist_rebinned.GetXaxis().GetBinLowEdge(xbin) >= hist_rebinned.xaxis.range_user[0] and hist_rebinned.GetXaxis().GetBinUpEdge(xbin) <= hist_rebinned.xaxis.range_user[1]:
                        hist_yproj = hist_rebinned.ProjectionY("", xbin, xbin).Clone()
                        hist_yproj = asrootpy(hist_yproj)
                        if hist_yproj.Integral() == 0:
                            continue

                        hmean = hist_yproj.GetMean()
                        hrms = hist_yproj.GetRMS()           

                        plotter.set_histo_style(hist_yproj, color=colors[i], title='%s-%s Mean=%.2f, RMS=%.2f' % (hist_rebinned.GetXaxis().GetBinLowEdge(xbin), hist_rebinned.GetXaxis().GetBinUpEdge(xbin), hmean, hrms))
                        plotter.plot(hist_yproj, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=reso_hists[kvar][obj][corr_type][2], ytitle=defyax, drawstyle='hist')
                        #plotter.save('test_yproj')
                        yprojections.append(hist_yproj)
                        yproj_norms.append(hist_yproj/hist_yproj.Integral())
                        i += 1
                        #if 'Reso_MTTbar' in obj: set_trace()
                        #set_trace()
        
                    #set_trace() 
                    #plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', hvar.split('/')[0], 'y_proj']))
                    plotter.overlay(yprojections, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=reso_hists[kvar][obj][corr_type][2], ytitle=defyax, drawstyle='hist')
                    plotter.save('%s_%s_yprojections' % (obj, corr_type))
                    plotter.overlay(yproj_norms, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=reso_hists[kvar][obj][corr_type][2], ytitle=defyax, drawstyle='hist')
                    plotter.save('%s_%s_yprojections_Norm' % (obj, corr_type))

                    #set_trace() 

                    continue


                    ## set hist range, get hist mean and rms 
                hvar_xmin = reso_hists[kvar][obj][corr_type][1]
                hvar_xmax = reso_hists[kvar][obj][corr_type][2]
                hvar_col = reso_hists[kvar][obj][corr_type][3]
                hvar_title = reso_hists[kvar][obj][corr_type][4]
                hvar_lstyle = reso_hists[kvar][obj][corr_type][5] # linestyle

                hist.xaxis.range_user = hvar_xmin, hvar_xmax
                hmean = hist.GetMean()
                hrms = hist.GetRMS()

                    ## set hist style, plot, and save
                plotter.set_histo_style(hist, xtitle='%s %s' % (hvar_title, xlabel), ytitle=defyax)
                plotter.plot(hist, drawstyle='hist')
                box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
                box1.Draw()
                plotter.save('Reso_%s_%s_%s' % (kvar, obj, corr_type))

                    ## set hist style to be compared with other type
                plotter.set_histo_style(hist, color=hvar_col, title='%s Mean=%.2f, RMS=%.2f' % (hvar_title, hmean, hrms), linestyle=hvar_lstyle )
                to_draw_normal.append(hist)

                #set_trace()   

                ## compare uncorreced/corrected hists
            if not to_draw_normal:
                continue
            plotter.overlay(to_draw_normal, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            plotter.save('Reso_%s_%s_Comparison' % (kvar, obj) )




def Parton_Acceptances_Plots(directory, subdir):
    #npartons = {
    #                    '2Partons' : ', 2 partons',
    #                    '3Partons' : ', 3 partons',
    #                    '4Partons' : ', 4 partons',
    #                    'All'      : ''
    #}

    #parton_acceptance_hists = {
    #                        '1D' : {
    #                                'xlabel' : 'M(t#bar{t}) Resolution [GeV]',
    #                                'title' : 'Reso_MTTbar_',
    #                                'hists' : {
    #                                            'Uncorrected' : ('', 'b'),
    #                                            'Corrected' : ('_Corrected_1D_E', 'r')
    #                                }
    #                        },
    #                        '2D' : {
    #                                'xlabel' : 'Gen M(t#bar{t}) [GeV]',
    #                                'ylabel' : 'M(t#bar{t}) Resolution [GeV]',
    #                                'hists' : {
    #                                            'Uncorrected' : 'Reso_MTTbar_vs_Gen_MTTbar_',
    #                                            'Corrected' : 'Reso_MTTbar_vs_Gen_MTTbar_'
    #                                }
    #                        }
    #}

    parton_acceptances = {
                '2Partons' :{
                            '1D' : {
                                    'xlabel' : 'M(t#bar{t}) Resolution [GeV], 2 partons',
                                    'title' : 'Reso_MTTbar_2Partons',
                                    'hists' : {
                                                'Uncorrected' : ('', 'b'),
                                                'Corrected_1D_E' : ('_Corrected_1D_E', 'r'),
                                                'Corrected_2D_E' : ('_Corrected_2D_E', 'g')
                                    }
                            },
                            '2D' : {
                                    'Gen_MTTbar': {
                                                'xlabel' : 'Gen M(t#bar{t}) [GeV], 2 partons',
                                                'ylabel' : 'M(t#bar{t}) Resolution [GeV]',
                                                'hists' : {
                                                            'Uncorrected' : 'Reso_MTTbar_vs_Gen_MTTbar_2Partons',
                                                            'Corrected_1D_E' : 'Reso_MTTbar_vs_Gen_MTTbar_2Partons_Corrected_1D_E',
                                                            'Corrected_2D_E' : 'Reso_MTTbar_vs_Gen_MTTbar_2Partons_Corrected_2D_E'
                                                }
                                    },
                                    'Gen_THadPt': {
                                                'xlabel' : 'Gen p_{T}(t_{h}) [GeV], 2 partons',
                                                'ylabel' : 'M(t#bar{t}) Frac. Resolution',
                                                'hists' : {
                                                           'Uncorrected' : 'Frac_MTTbar_vs_Gen_THadPt_2Partons',
                                                           'Corrected_1D_E' : 'Frac_MTTbar_vs_Gen_THadPt_2Partons_Corrected_1D_E',
                                                           'Corrected_2D_E' : 'Frac_MTTbar_vs_Gen_THadPt_2Partons_Corrected_2D_E'
                                                }
                                    }
                            }
                },
                '3Partons' :{
                            '1D' : {
                                    'xlabel' : 'M(t#bar{t}) Resolution [GeV], 3 partons',
                                    'title' : 'Reso_MTTbar_3Partons',
                                    'hists' : {
                                                'Uncorrected' : ('', 'b'),
                                                'Corrected_1D_E' : ('_Corrected_1D_E', 'r'),
                                                'Corrected_2D_E' : ('_Corrected_2D_E', 'g')
                                    }
                            },
                            '2D' : {
                                    'Gen_MTTbar': {
                                                'xlabel' : 'Gen M(t#bar{t}) [GeV], 3 partons',
                                                'ylabel' : 'M(t#bar{t}) Resolution [GeV]',
                                                'hists' : {
                                                            'Uncorrected' : 'Reso_MTTbar_vs_Gen_MTTbar_3Partons',
                                                            'Corrected_1D_E' : 'Reso_MTTbar_vs_Gen_MTTbar_3Partons_Corrected_1D_E',
                                                            'Corrected_2D_E' : 'Reso_MTTbar_vs_Gen_MTTbar_3Partons_Corrected_2D_E'
                                                }
                                    },
                                    'Gen_THadPt': {
                                                'xlabel' : 'Gen p_{T}(t_{h}) [GeV], 3 partons',
                                                'ylabel' : 'M(t#bar{t}) Frac. Resolution',
                                                'hists' : {
                                                           'Uncorrected' : 'Frac_MTTbar_vs_Gen_THadPt_3Partons',
                                                           'Corrected_1D_E' : 'Frac_MTTbar_vs_Gen_THadPt_3Partons_Corrected_1D_E',
                                                           'Corrected_2D_E' : 'Frac_MTTbar_vs_Gen_THadPt_3Partons_Corrected_2D_E'
                                                }
                                    }
                            }
                },
                '4Partons' :{
                            '1D' : {
                                    'xlabel' : 'M(t#bar{t}) Resolution [GeV], 4 partons',
                                    'title' : 'Reso_MTTbar_4Partons',
                                    'hists' : {
                                                'Uncorrected' : ('', 'b'),
                                                'Corrected_1D_E' : ('_Corrected_1D_E', 'r'),
                                                'Corrected_1D_E' : ('_Corrected_2D_E', 'g')
                                    }
                            },
                            '2D' : {
                                    'Gen_MTTbar': {
                                                'xlabel' : 'Gen M(t#bar{t}) [GeV], 4 partons',
                                                'ylabel' : 'M(t#bar{t}) Resolution [GeV]',
                                                'hists' : {
                                                            'Uncorrected' : 'Reso_MTTbar_vs_Gen_MTTbar_4Partons',
                                                            'Corrected_1D_E' : 'Reso_MTTbar_vs_Gen_MTTbar_4Partons_Corrected_1D_E',
                                                            'Corrected_2D_E' : 'Reso_MTTbar_vs_Gen_MTTbar_4Partons_Corrected_2D_E'
                                                }
                                    },
                                    'Gen_THadPt': {
                                                'xlabel' : 'Gen p_{T}(t_{h}) [GeV], 4 partons',
                                                'ylabel' : 'M(t#bar{t}) Frac. Resolution',
                                                'hists' : {
                                                           'Uncorrected' : 'Frac_MTTbar_vs_Gen_THadPt_4Partons',
                                                           'Corrected_1D_E' : 'Frac_MTTbar_vs_Gen_THadPt_4Partons_Corrected_1D_E',
                                                           'Corrected_2D_E' : 'Frac_MTTbar_vs_Gen_THadPt_4Partons_Corrected_2D_E'
                                                }
                                    }
                            }
                },
                'All' : {
                            '1D' : {
                                    'xlabel' : 'M(t#bar{t}) Resolution [GeV]',
                                    'title' : 'Reso_MTTbar_All',
                                    'hists' : {
                                                'Uncorrected' : ('', 'b'),
                                                'Corrected_1D_E' : ('_Corrected_1D_E', 'r'),
                                                'Corrected_2D_E' : ('_Corrected_2D_E', 'g')
                                    }
                            },
                            '2D' : {
                                    'Gen_MTTbar': {
                                                'xlabel' : 'Gen M(t#bar{t}) [GeV]',
                                                'ylabel' : 'M(t#bar{t}) Resolution [GeV]',
                                                'hists' : {
                                                           'Uncorrected' : 'Reso_MTTbar_vs_Gen_MTTbar_All',
                                                           'Corrected_1D_E' : 'Reso_MTTbar_vs_Gen_MTTbar_All_Corrected_1D_E',
                                                           'Corrected_2D_E' : 'Reso_MTTbar_vs_Gen_MTTbar_All_Corrected_2D_E'
                                                }
                                    },
                                    'Gen_THadPt': {
                                                'xlabel' : 'Gen p_{T}(t_{h}) [GeV]',
                                                'ylabel' : 'M(t#bar{t}) Frac. Resolution',
                                                'hists' : {
                                                           'Uncorrected' : 'Frac_MTTbar_vs_Gen_THadPt_All',
                                                           'Corrected_1D_E' : 'Frac_MTTbar_vs_Gen_THadPt_All_Corrected_1D_E',
                                                           'Corrected_2D_E' : 'Frac_MTTbar_vs_Gen_THadPt_All_Corrected_2D_E'
                                                }
                                    }
                            }
                }
    }

    #test_vars = [ # input_dir, variable, xlabel, ylabel, title, color, output dir
    #    ('Reso_MTTbar_All', 'Uncorrected M(t#bar{t}) Resolution [GeV]', defyax, '', 'b')
    #]

    #for var, xlabel, ylabel, htitle, color in test_vars:
    #
    #                    #### normal plots
    #    varname = title+parton_acceptances[nparts][ndim]['hists'][corr_type][0]
    #    hname = '/'.join([directory, 'Post_Alpha_Correction', 'Resolution', 'Parton_Acceptance', varname])
    #    
    #    hist = asrootpy(myfile.Get(hname)).Clone()
    #    
    #    if hist.Integral() == 0:
    #        continue
    #    
    #    hvar_col = parton_acceptances[nparts][ndim]['hists'][corr_type][1]
    #    
    #    #hist.xaxis.range_user = hvar_xmin, hvar_xmax
    #    hmean = hist.GetMean()
    #    hrms = hist.GetRMS()
    #    
    #        ## set hist style, plot, and save
    #    plotter.set_histo_style(hist, xtitle='%s %s' % (corr_type, xlabel), ytitle=defyax)
    #    plotter.plot(hist, drawstyle='hist')
    #    box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
    #    box1.Draw()
    #    plotter.save(varname)

    #set_trace()
    for nparts in parton_acceptances.keys():
 
        for ndim in parton_acceptances[nparts].keys():

            if ndim == '1D':

                if nparts == 'All' :
                    plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances']))
                else:
                    plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances', nparts]))

                xlabel = parton_acceptances[nparts][ndim]['xlabel']
                title = parton_acceptances[nparts][ndim]['title']
                to_draw_normal = []

                for corr_type in parton_acceptances[nparts][ndim]['hists'].keys():


                        ### create hists based on perm category (Correct b, wrong b, etc...)
                    to_draw = []
                    for cat in Correct_WJet_Categories+Wrong_WJet_Categories+['NO_MP']:
                        plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances', nparts, 'Perm_Categories', corr_type]))
                        evt_col=styles[cat]['fillcolor']
                        evt_type = styles[cat]['name']
                        varname = title+parton_acceptances[nparts][ndim]['hists'][corr_type][0]
                        cat_hname = '/'.join([directory, 'Post_Alpha_Correction', cat, 'Resolution', 'Parton_Acceptance', varname])
                        cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

                        #set_trace()   
                        #cat_hist = RebinView.rebin(cat_hist, new_bins)
                        #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                        if cat_hist.Integral() == 0:
                            continue
    
                        cat_hist_mean = cat_hist.GetMean() 
                        cat_hist_rms = cat_hist.GetRMS() 

                        plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                        cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                        plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                        to_draw.append(cat_hist)

                        box = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (cat_hist_mean, cat_hist_rms), position='NE')
                        box.Draw()
                        plotter.save('Post_Alpha_%s_%s_%s' % (corr_type, title, cat))
 
                    if not to_draw:
                        continue
                    stack, norm_stack, ratio = fncts.stack_plots(to_draw)
                    plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                    #set_trace()
                    
                    box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (sum(to_draw).GetMean(), sum(to_draw).GetRMS()), position='NE')
                    box1.Draw()

                    plotter.save('Post_Alpha_%s_%s_Stack' % (corr_type, title))
    
                    plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                    box1.Draw()
                    plotter.save('Post_Alpha_%s_%s_Stack_Norm' % (corr_type, title))


                    #### create hists based on gen category (Correct bhad, correct blep, wrong b, etc...)
                    #to_draw = []
                    #for cat in Gen_Categories:
                    #    plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances', nparts, 'Gen_Categories', corr_type]))
                    #    evt_col=styles[cat]['fillcolor']
                    #    evt_type = styles[cat]['name']
                    #    varname = title+parton_acceptances[nparts][ndim]['hists'][corr_type][0]
                    #    cat_hname = '/'.join([directory, 'Post_Alpha_Correction', cat, 'Resolution', 'Parton_Acceptance', varname])
                    #    cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                    #    #cat_hist = RebinView.rebin(cat_hist, new_bins)
                    #    #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                    #    if cat_hist.Integral() == 0:
                    #        continue
    
                    #    cat_hist_mean = cat_hist.GetMean() 
                    #    cat_hist_rms = cat_hist.GetRMS() 
    
                    #    #set_trace()

                    #    plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                    #    cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                    #    plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                    #    to_draw.append(cat_hist)

                    #    box = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (cat_hist_mean, cat_hist_rms), position='NE')
                    #    box.Draw()
                    #    plotter.save('Post_Alpha_%s_%s_%s' % (corr_type, title, cat))

    
                    #if not to_draw:
                    #    continue
                    #stack, norm_stack, ratio = fncts.stack_plots(to_draw)
                    #plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                    #box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (sum(to_draw).GetMean(), sum(to_draw).GetRMS()), position='NE')
                    #box1.Draw()

                    #plotter.save('Post_Alpha_%s_%s_Stack' % (corr_type, title))
    
                    #plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                    #box1.Draw()
                    #plotter.save('Post_Alpha_%s_%s_Stack_Norm' % (corr_type, title))

                    #set_trace()

                        #### normal plots
                    plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances', nparts]))
                    varname = title+parton_acceptances[nparts][ndim]['hists'][corr_type][0]
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
                    to_draw_normal.append(hist)


                    ## compare uncorreced/corrected hists
                if not to_draw_normal:
                    continue
                plotter.overlay(to_draw_normal, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
                plotter.save('%s_Comparison' % title )


            elif ndim == '2D':
                
                #set_trace()
                for gen_var_type in parton_acceptances[nparts][ndim].keys():

                    if gen_var_type == 'Gen_MTTbar':
                        yproj_label_range = 'M(t#bar{t})'
                    else:
                        yproj_label_range = 'p_{T}(t_{h})'

                    xlabel = parton_acceptances[nparts][ndim][gen_var_type]['xlabel']
                    ylabel = parton_acceptances[nparts][ndim][gen_var_type]['ylabel']

                    for corr_type in parton_acceptances[nparts][ndim][gen_var_type]['hists'].keys():
                        if nparts == 'All' :
                            plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances']))
                        else:
                            plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances', nparts]))

                        varname = parton_acceptances[nparts][ndim][gen_var_type]['hists'][corr_type]

                        if nparts == 'All':
                            hname = '/'.join([directory, 'Post_Alpha_Correction', 'Resolution', 'Parton_Acceptance', varname])
                        else:
                            hname = '/'.join([directory, 'Post_Alpha_Correction', 'Resolution', 'Parton_Acceptance', varname])
 
                        hist = asrootpy(myfile.Get(hname)).Clone()

                        if hist.Integral() == 0:
                            continue

                        plotter.set_histo_style(hist, xtitle=xlabel, ytitle='%s %s' % (corr_type, ylabel))
                        plotter.plot(hist, drawstyle='colz')

                        mean_x = hist.ProjectionX().GetMean()
                        rms_x = hist.ProjectionX().GetRMS()
                        mean_y = hist.ProjectionY().GetMean()
                        rms_y = hist.ProjectionY().GetRMS()

                        box2 = plotter.make_text_box('X Mean=%.2f, RMS=%.2f\nY Mean=%.2f, RMS=%.2f' % (mean_x, rms_x, mean_y, rms_y), position='NE')
                        box2.Draw()

                        plotter.save( varname )

                        #set_trace()

                            #### split 2D plots into bins of gen M(ttbar) or pT(thad)
                        plotter.set_subdir('/'.join([subdir, 'Post_Alpha', 'Resolution', 'Partons_Acceptances', nparts, 'yproj', corr_type]))

                            ## make list of cumulative integrals for each bin
                        max_nbins = 5
                        round_val = 0
                        if hist.Integral() < 10000:
                            round_val = -2
                        elif hist.Integral() < 100000:
                            round_val = -3
                        else:
                            round_val = -4

                        starting_min_binval = round( hist.Integral(), round_val )/max_nbins # minimum number of allowed events for bins
                        min_binval =  starting_min_binval # minimum number of allowed events for bins
                        binval = 0
                        bins = []

                        bin_cum_ints = [ hist.Integral(1, xbin, 1, hist.GetNbinsY()) for xbin in range(1, hist.GetNbinsX()+1) ]
                        first_nonzero_indx = [ n for n,i in enumerate(bin_cum_ints) if i > 0 ][0]+1
                        last_nonzero_indx = [ n for n,i in enumerate(bin_cum_ints) if i == hist.Integral() ][0]+2

                        bins.append(first_nonzero_indx)

                            ### start loop
                        while round(hist.Integral()-min_binval, round_val) > starting_min_binval:
                        #while round(hist.Integral()-binval, round_val) > starting_min_binval:
                            Bin = [ n for n,i in enumerate(bin_cum_ints) if i > min_binval ][0]+2
                            bins.append(Bin)
                            binval = bin_cum_ints[Bin-2] #find integral up to that bin
                            min_binval = starting_min_binval + binval # update minimum integral needed

                            #set_trace()

                        bins.append(last_nonzero_indx)


                        #set_trace()
                        #bin_ints = [ bin_cum_ints[i-2] for i in bins ]
                        xbins = np.array( [hist.GetXaxis().GetBinLowEdge(i) for i in bins ] )
                        ybins = np.linspace(hist.GetYaxis().GetBinLowEdge(1), hist.GetYaxis().GetBinUpEdge(hist.GetYaxis().GetNbins()), hist.GetYaxis().GetNbins()+1 ) # ybinning remains unchanged
                        hist_rebinned = RebinView.newRebin2D(hist, xbins, ybins)

                            ### split 2D gen pT(thad) plots at gen pT(thad) > 450
                        if gen_var_type == 'Gen_THadPt':
                            hist_low_xbin_edges = [hist.GetXaxis().GetBinLowEdge(i+1) for i in range(hist.GetNbinsX()) ]
                            pt_cut = 350
                            pt_bin_350 = [ n for n,i in enumerate(hist_low_xbin_edges) if i >= pt_cut ][0]+1
                            xbin_indices = [ first_nonzero_indx, pt_bin_350, last_nonzero_indx ]
                            xbin = np.array( [hist.GetXaxis().GetBinLowEdge(i) for i in xbin_indices ] )

                            pt_350_hist = RebinView.newRebin2D(hist, xbin, ybins)

                            #set_trace()

                                ## yprojections 
                            for binx in range(1, pt_350_hist.GetNbinsX() + 1):
                                hist_yproj = pt_350_hist.ProjectionY("", binx, binx)

                                if hist_yproj.Integral() == 0:
                                    continue

                                #set_trace()

                                hist_yproj.Draw('hist')
                                hist_yproj.SetXTitle('%s %s' % (corr_type, ylabel))
                                hist_yproj.SetYTitle(defyax)

                                mean_proj = hist_yproj.GetMean()
                                rms_proj = hist_yproj.GetRMS()
                                
                                genmttbar_range = '%.0f #leq Gen %s #leq %.0f' % (pt_350_hist.GetXaxis().GetBinLowEdge(binx), yproj_label_range, pt_350_hist.GetXaxis().GetBinUpEdge(binx))

                                box1 = plotter.make_text_box('%s\nMean=%.2f\nRMS=%.2f' % (genmttbar_range, mean_proj, rms_proj), position='NE')
                                box1.Draw()

                                plotter.save('%s_%sSplit_Part%s' % (varname, pt_cut, binx) )

                            #set_trace()


                            ## yprojections for all 2D hists
                        for xbin in range(hist_rebinned.GetNbinsX() + 1):
                            hist_yproj = hist_rebinned.ProjectionY("", xbin, xbin)

                            if hist_yproj.Integral() == 0:
                                continue

                            #set_trace()

                            hist_yproj.Draw('hist')
                            hist_yproj.SetXTitle('%s %s' % (corr_type, ylabel))
                            hist_yproj.SetYTitle(defyax)

                            mean_proj = hist_yproj.GetMean()
                            rms_proj = hist_yproj.GetRMS()
                            
                            genmttbar_range = '%.0f #leq Gen %s #leq %.0f' % (hist_rebinned.GetXaxis().GetBinLowEdge(xbin), yproj_label_range, hist_rebinned.GetXaxis().GetBinUpEdge(xbin))

                            box1 = plotter.make_text_box('%s\nMean=%.2f\nRMS=%.2f' % (genmttbar_range, mean_proj, rms_proj), position='NE')
                            box1.Draw()

                            plotter.save('%s_%s' % (varname, hist_rebinned.GetXaxis().GetBinLowEdge(xbin)) )
                        #set_trace()



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


def Post_Alpha_Plots( plot ):

    reco_dir, reco_subdir = '', ''


    if args.perm == 'Event':
        reco_dir = '3J/nosys/Lost_BP'
        reco_subdir = '3J_Event_Plots/Lost_BP'

    elif args.perm == 'Matched':
        reco_dir = '3J/nosys/Matched_Perm'
        reco_subdir = '3J_Event_Plots/Matched_Perm'


    if plot == 'Reconstruction':
        Reconstruction_Plots( reco_dir , reco_subdir )
    elif plot == 'Resolution':
        Reso_Plots( reco_dir , reco_subdir )

#    if plot == 'Post_Alpha':
#        ### create post alpha correction plots for lost-jet events
#        post_alpha_corrections( reco_dir , reco_subdir, classes )
#
        #set_trace()
    


#####################################################################################################

if args.plot == 'All':
    Post_Alpha_Plots('Reconstruction')
    Post_Alpha_Plots('Resolution')
    Post_Alpha_Plots('Partons_Acceptance')
else:
    Post_Alpha_Plots(args.plot)

