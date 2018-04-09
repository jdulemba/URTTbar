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
        defaults = {'save' : {'png' : True, 'pdf' : True}},
        styles = {
            'right_whad' : styles['right_whad *'],
            'qcd' : styles['QCD*'],
            'single_top' : styles['single*'],
            'wrong_whad' : styles['wrong_whad *'],
            'nonsemi_tt' : styles['nonsemi_tt *'],
            'vjets' : styles['[WZ]Jets*'],
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


tdirs = ['ditag', 'notag', 'subtag', 'leadtag']

for tdir in tdirs:
    shapes = datacard_file.Get(tdir)
    categories = [i.name for i in shapes.keys()]

    right_whad_yields = right_whad_yields_file['right_whad_yields'][tdir]

    non_sys_samples = []
    data_samples = []
    for cat in categories:

        ### only get hists that aren't systematics or data
        if 'Up' not in cat and 'Down' not in cat and 'data' not in cat:
            non_sys_samples.append(cat)

        if 'data' in cat:
            data_samples.append(cat)
    data = shapes.Get('data_obs')
    samples = [shapes.Get(i) for i in non_sys_samples]

    rwhad_indx = [i for i,j in enumerate(samples) if j.name == 'right_whad']
    samples[rwhad_indx[0]] = samples[rwhad_indx[0]]*right_whad_yields
    samples[rwhad_indx[0]].name = 'right_whad'
        
    #non_sys_samples.sort(key=ordering)
    stack = plotter.create_stack(*samples, sort=False)
    legend = LegendDefinition(position='NE')
    plotter.overlay_and_compare(
        [stack], data,
        legend_def = legend,
        lower_y_range=0.5,
        method='datamc',
        xtitle='%s #lambda_{M}' % tdir,
        ytitle='Events',
        writeTo=tdir,
        )

