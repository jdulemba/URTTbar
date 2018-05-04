import rootpy.io as io
import rootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
rootpy.log["/data_views"].setLevel(rootpy.log.WARNING)
log = rootpy.log["/make_postfit_plots"]
log.setLevel(rootpy.log.INFO)
from rootpy.plotting import Hist
from argparse import ArgumentParser
import os, glob, sys
import re
import ROOT
from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.Utilities.prettyjson as prettyjson
from URAnalysis.Utilities.quad import quad
from styles import styles
from pdb import set_trace
from labels import set_pretty_label
from URAnalysis.Utilities.tables import Table
from argparse import ArgumentParser
import numpy as np
import logging
import matplotlib.pyplot as plt

parser = ArgumentParser()
parser.add_argument('--wp')
parser.add_argument('--eras', default='')

args = parser.parse_args()

dirnames = {'All' : 'All_Runs', 'B' : 'Run_B', 'CtoE' : 'Run_CtoE', 'EtoF' : 'Run_EtoF'}

dir_name = dirnames.get(args.eras, 'Exit')

if dir_name == 'Exit':
    logging.error('Not a valid era to choose from.')
    sys.exit()

#set_trace()

input_dir = '%s/plots/%s/ctageff/%s/mass_discriminant/%s' % (os.environ['URA_PROJECT'],os.environ['jobid'], dir_name, args.wp)

if args.eras == 'All':
    plotter = BasePlotter(
        '%s/datacard' % input_dir,
        defaults = {'save' : {'png' : True, 'pdf' : False}},
        styles = {
            'right_whad*' : styles['right_whad *'],
            'qcd*' : styles['QCD*'],
            'single_top*' : styles['single*'],
            'wrong_whad*' : styles['wrong_whad *'],
            'nonsemi_tt*' : styles['nonsemi_tt *'],
            'vjets*' : styles['[WZ]Jets*'],
            'data_obs' : styles['data*'],
            'Total signal+background *' : {
                'legendstyle' : 'f',
                'drawstyle' : 'PE2',
                'linewidth' : 0,
                'title' : "uncertainty",
                'markersize' : 0,
                'fillcolor' : 1,
                'fillstyle' : 3013
                }
            }
    )
#elif args.eras == 'B':
#    plotter = BasePlotter(
#        '%s/datacard' % input_dir,
#        defaults = {'save' : {'png' : True, 'pdf' : True}},
#        styles = {
#            'right_whad' : styles['right_whad *'],
#            'qcd' : styles['QCD*'],
#            'single_top' : styles['single*'],
#            'wrong_whad' : styles['wrong_whad *'],
#            'nonsemi_tt' : styles['nonsemi_tt *'],
#            'vjets' : styles['[WZ]Jets*'],
#            'data_obs' : styles['data*B*'],
#            'Total signal+background *' : {
#                'legendstyle' : 'f',
#                'drawstyle' : 'PE2',
#                'linewidth' : 0,
#                'title' : "uncertainty",
#                'markersize' : 0,
#                'fillcolor' : 1,
#                'fillstyle' : 3013
#                }
#            }
#    )
#elif args.eras == 'CtoE':
#    plotter = BasePlotter(
#        '%s/datacard' % input_dir,
#        defaults = {'save' : {'png' : True, 'pdf' : True}},
#        styles = {
#            'right_whad' : styles['right_whad *'],
#            'qcd' : styles['QCD*'],
#            'single_top' : styles['single*'],
#            'wrong_whad' : styles['wrong_whad *'],
#            'nonsemi_tt' : styles['nonsemi_tt *'],
#            'vjets' : styles['[WZ]Jets*'],
#            'data_obs' : styles['data*CtoE*'],
#            'Total signal+background *' : {
#                'legendstyle' : 'f',
#                'drawstyle' : 'PE2',
#                'linewidth' : 0,
#                'title' : "uncertainty",
#                'markersize' : 0,
#                'fillcolor' : 1,
#                'fillstyle' : 3013
#                }
#            }
#    )
#elif args.eras == 'EtoF':
#    plotter = BasePlotter(
#        '%s/datacard' % input_dir,
#        defaults = {'save' : {'png' : True, 'pdf' : True}},
#        styles = {
#            'right_whad' : styles['right_whad *'],
#            'qcd' : styles['QCD*'],
#            'single_top' : styles['single*'],
#            'wrong_whad' : styles['wrong_whad *'],
#            'nonsemi_tt' : styles['nonsemi_tt *'],
#            'vjets' : styles['[WZ]Jets*'],
#            'data_obs' : styles['data*EtoF*'],
#            'Total signal+background *' : {
#                'legendstyle' : 'f',
#                'drawstyle' : 'PE2',
#                'linewidth' : 0,
#                'title' : "uncertainty",
#                'markersize' : 0,
#                'fillcolor' : 1,
#                'fillstyle' : 3013
#                }
#            }
#    )
else:
    logging.error('Not a valid era to choose from.')
    sys.exit()

datacard_file = io.root_open(
    '%s/datacard.root' % input_dir
)

right_whad_yields_file = prettyjson.loads(open('%s/datacard.json' % input_dir).read())
if not 'right_whad_yields' in right_whad_yields_file.keys():
    logging.error('Raw right_whad yields not in datacard.json.')
    sys.exit()

## things for plotting data/mc ratios
sys_styles = {
            'MTOPUp' : {
                'linestyle' : '-',
                'color' : 'red',
                },
            'MTOPDown' : {
                'linestyle' : '-.',
                'color' : 'red',
                },
            'JERUp' : {
                'linestyle' : '-',
                'color' : 'blue',
                },
            'JERDown' : {
                'linestyle' : '-.',
                'color' : 'blue',
                },
            'JESUp' : {
                'linestyle' : '-',
                'color' : 'green',
                },
            'JESDown' : {
                'linestyle' : '-.',
                'color' : 'green',
                },
            'PDFUp' : {
                'linestyle' : '-',
                'color' : 'magenta',
                },
            'PDFDown' : {
                'linestyle' : '-.',
                'color' : 'magenta',
                },
}
mc_xbins = np.linspace(6,20, 100)
##

data_mc_ratios = { #dict of data/MC ratios
    'ditag' : {},
    'notag' : {},
    'subtag' : {},
    'leadtag' : {},
    'Inclusive' : {},
}
tdirs = ['ditag', 'notag', 'subtag', 'leadtag']


    ## create stack plot for nominal data/MC plots
def stack_nom_MC(MC_hists, data_hist, tdir): # (array of MC hists for each topology, data hist, tag category)
        ## rescale right_whad hists
    rwhad_indx = [i for i,j in enumerate(MC_hists) if j.name == 'right_whad']
    MC_hists[rwhad_indx[0]] = MC_hists[rwhad_indx[0]]*right_whad_yields
    MC_hists[rwhad_indx[0]].name = 'right_whad'
    
    #non_sys_samples.sort(key=ordering)
    stack = plotter.create_stack(*MC_hists, sort=False)
    legend = LegendDefinition(position='NE')
    plotter.overlay_and_compare(
        [stack], data_hist,
        legend_def = legend,
        lower_y_range=0.5,
        method='datamc',
        xtitle='%s #lambda_{M}' % tdir,
        ytitle='Events',
        writeTo=tdir,
        )

    #set_trace()
    data_mc_ratios[tdir]['Nosys'] = data_hist.Integral()/stack.Integral()
    #plt.plot(12., data_mc_ratios[tdir]['Nosys'], label='Nominal', marker='.', color='black')


    ## create stack plot for systematic data/MC plots
def stack_systematic_MC(mc_hist, data_hist, tdir, title): # (array of MC hists for each topology, data hist, tag category, shift title)
        ## rescale right_whad hists
    rwhad_indx = [i for i,j in enumerate(mc_hist) if j.name == 'right_whad_%s' % title]
    mc_hist[rwhad_indx[0]] = mc_hist[rwhad_indx[0]]*right_whad_yields
    mc_hist[rwhad_indx[0]].name = 'right_whad_%s' % title
    
    mc_stack = plotter.create_stack(*mc_hist, sort=False)
    legend = LegendDefinition(position='NE')
    plotter.overlay_and_compare(
        [mc_stack], data_hist,
        legend_def = legend,
        lower_y_range=0.0,
        method='datamc',
        xtitle='%s %s #lambda_{M}' % (tdir, title),
        ytitle='Events',
        writeTo='Sys_Shifts/%s/%s_%s' % (tdir, tdir, title)
        )
    
    #plt.plot(mc_xbins, np.ones(100)*data_hist.Integral()/mc_stack.Integral(), label=title, linestyle=sys_styles[title]['linestyle'], color=sys_styles[title]['color'])
    data_mc_ratios[tdir][title] = data_hist.Integral()/mc_stack.Integral()

def get_sys_hists(categories, sys):
    sys_samples = [cat for cat in categories if sys in cat]

    ## lists of non-ttjets MC sammples
    qcd_vjets_samples = [cat for cat in categories if 'qcd' in cat or 'vjets' in cat] #qcd and v+jets
    singlet_samples = [cat for cat in categories if cat=='single_top'] #single top

    ## get hists for systematics samples (systematics+MC samples they don't include)
    sys_hists = []
    if 'MTOP' in sys or 'PDF' in sys:
        sys_hists = [shapes.Get(i) for i in sys_samples+singlet_samples+qcd_vjets_samples]
    if 'JES' in sys or 'JER' in sys:
        sys_hists = [shapes.Get(i) for i in sys_samples+qcd_vjets_samples]

    return sys_hists


    ### makes plots for combined (inclusive) plots
fig = plt.figure()
#sys_list = ['MTOPUp', 'Nosys']
#for sys in sys_list:
for sys in sys_styles.keys()+['Nosys']:
    all_indiv_MC_hists = []
    all_indiv_data_hists = []
    inclusive_MC_hists = []

    for tdir in tdirs:
        shapes = datacard_file.Get(tdir)
        categories = [i.name for i in shapes.keys()]
    
        right_whad_yields = right_whad_yields_file['right_whad_yields'][tdir]
    
        ## check if directory for sys shifts exists, make if it doesn't
        if not os.path.exists('%s/datacard/Sys_Shifts/%s' % (input_dir, tdir) ):
            os.makedirs('%s/datacard/Sys_Shifts/%s' % (input_dir, tdir))

        ## list of MC hists with no systematics
        if sys == 'Nosys':
            MC_samples = [cat for cat in categories if 'Up' not in cat and 'Down' not in cat and 'data' not in cat]
            MC_hists = [shapes.Get(i) for i in MC_samples]
                ## rescale right_whad hists
            rwhad_indx = [i for i,j in enumerate(MC_hists) if j.name == 'right_whad']
            MC_hists[rwhad_indx[0]] = MC_hists[rwhad_indx[0]]*right_whad_yields
            MC_hists[rwhad_indx[0]].name = 'right_whad'
            all_indiv_MC_hists.append(MC_hists)

        else:
            MC_hists = get_sys_hists(categories, sys)
                ## rescale right_whad hists
            rwhad_indx = [i for i,j in enumerate(MC_hists) if j.name == 'right_whad_%s' % sys]
            MC_hists[rwhad_indx[0]] = MC_hists[rwhad_indx[0]]*right_whad_yields
            MC_hists[rwhad_indx[0]].name = 'right_whad_%s' % sys
            all_indiv_MC_hists.append(MC_hists)


        data_samples = [cat for cat in categories if 'data' in cat]

        ## get hists for data and MC samples
        data = shapes.Get('data_obs')
        
        all_indiv_data_hists.append(data)
        #set_trace()

    for i in range(len(all_indiv_MC_hists[0])):
        inclusive_MC_hists.append(sum([all_indiv_MC_hists[j][i]  for j in range(len(all_indiv_MC_hists))]))
        inclusive_MC_hists[i].name = inclusive_MC_hists[i].title
        #set_trace()

    inclusive_data = sum(all_indiv_data_hists)
    inclusive_data.name = inclusive_data.title

    stack = plotter.create_stack(*inclusive_MC_hists, sort=False)
    legend = LegendDefinition(position='NE')
    plotter.overlay_and_compare(
        [stack], inclusive_data,
        legend_def = legend,
        lower_y_range=0.5,
        method='datamc',
        xtitle='%s Inclusive #lambda_{M}' % sys,
        ytitle='Events',
        writeTo='Sys_Shifts/%s_Inclusive' % sys,
        )

    ## get data/mc ratio for each bin
    data_mc_ratios['Inclusive'][sys] = [inclusive_data.Integral(bins+1,bins+1)/sum( stack[i].Integral(bins+1,bins+1) for i in range(len(stack)) ) for bins in range(inclusive_data.GetXaxis().GetNbins())]

    if sys == 'Nosys':
        xvals = [inclusive_data.GetXaxis().GetBinCenter(bins+1) for bins in range(inclusive_data.GetXaxis().GetNbins())]
        yvals = data_mc_ratios['Inclusive'][sys]
        plt.plot(xvals, yvals, label=sys, linestyle='None', marker='.', color='black')
    else:
        xvals = np.linspace(stack[0].GetXaxis().GetXmin(), stack[0].GetXaxis().GetXmax(), 100)
            ## get number of points within each bin
        ybin_ranges = [len([i for (i,j) in enumerate(xvals) if xvals[i] <= stack[0].GetXaxis().GetBinUpEdge(bins+1) and xvals[i] >= stack[0].GetXaxis().GetBinLowEdge(bins+1)]) for bins in range(stack[0].GetXaxis().GetNbins())]
        ybin_values = [np.ones(ybin_ranges[bins])*data_mc_ratios['Inclusive'][sys][bins] for bins in range(len(ybin_ranges))]
        yvals = np.concatenate(ybin_values).ravel()
        plt.plot(xvals, yvals, label=sys, linestyle=sys_styles[sys]['linestyle'], color=sys_styles[sys]['color'])
        #set_trace()
    #set_trace()

plt.xlabel('Inclusive $\lambda_{M}$')
plt.ylabel('data/MC')
plt.xlim(6.0, 20.0)
plt.grid()
plt.legend(loc='upper right',fontsize=8, numpoints=1)
fig.savefig('%s/datacard/Sys_Shifts/Inclusive_sys_ratios' % input_dir)
plt.close(fig)
#set_trace()


    ### makes plots for individual tags
for tdir in tdirs:
    shapes = datacard_file.Get(tdir)
    categories = [i.name for i in shapes.keys()]

    right_whad_yields = right_whad_yields_file['right_whad_yields'][tdir]

    ## list of MC hists with no systematics and then get hists
    non_sys_samples = [cat for cat in categories if 'Up' not in cat and 'Down' not in cat and 'data' not in cat]
    nom_MC_hists = [shapes.Get(i) for i in non_sys_samples]

    ## list of hists for data and then get hists
    data_samples = [cat for cat in categories if 'data' in cat]
    data = shapes.Get('data_obs')

    ### create stack plot for nominal data/MC plots
    stack_nom_MC(nom_MC_hists, data, tdir)
    #set_trace()

    fig = plt.figure()

    ## create stack plots for all systematics
    for sys in sys_styles.keys():
        sys_hists = get_sys_hists(categories, sys)
        stack_systematic_MC(sys_hists, data, tdir, sys) 
        #set_trace()

    plt.xlabel('%s $\lambda_{M}$' % tdir)
    plt.ylabel('data/MC')
    plt.xlim(6.0, 20.0)
    plt.grid()
    plt.legend(loc='upper right',fontsize=8, numpoints=1)
    fig.savefig('%s/datacard/%s_sys_ratios' % (input_dir, tdir))
    plt.close(fig)

    #set_trace()
    with open('%s/datacard/sys_ratios.json' % input_dir, 'w') as f:
        f.write(prettyjson.dumps(data_mc_ratios))
