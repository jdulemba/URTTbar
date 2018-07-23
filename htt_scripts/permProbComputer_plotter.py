'''
permProbComputer Analyzer Plotter macro
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os, glob, sys
import logging
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
import rootpy.plotting as plotting
#from rootpy.plotting import views, Graph, Hist
from argparse import ArgumentParser
import URAnalysis.Utilities.roottools as urtools
import ROOT
from URAnalysis.Utilities.tables import latex_table
from URAnalysis.PlotTools.views.RebinView import RebinView
import argparse
import matplotlib.pyplot as plt
import functions as fncts
views = plotting.views


analyzer = 'permProbComputer'

parser = argparse.ArgumentParser(description='Create plots using files from permProbComputer.')
parser.add_argument('--fname', help='Choose file to use')
args = parser.parse_args()


jobid = jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']

prob_file = args.fname
#prob_file = 'prob_prob_test.root' #inputs/2017Aug24/INPUT/prob_ttJets_3J.root
#prob_file = 'prob_htt_perm_baseline_j20_ttJetsAll.root'
#prob_file = 'prob_ttJets_3J.root' #inputs/2017Aug24/INPUT/prob_ttJets_3J.root

##### check if correct prob_file is present
dir_f_list = []
[dir_f_list.append(f) for f in os.listdir('%s/inputs/%s/INPUT' % (project, jobid))]
if not prob_file in dir_f_list:
    print "File %s isn't in the directory inputs/%s/INPUT. It has to be to run this!" % (prob_file, jobid)
    sys.exit()

myfile = root_open('%s/inputs/%s/INPUT/%s' % (project, jobid, prob_file), 'read')
normfile = views.NormalizeView(root_open('%s/inputs/%s/INPUT/%s' % (project, jobid, prob_file), 'read'))

lumis = glob.glob('inputs/%s/*.lumi' % jobid)
intlumi = (sum([float(open(lumi, 'read').readline()) for lumi in lumis if 'data_Single' in lumi])/2)/1000 #integrated lumi for 2016 in fb-1

plotter = BasePlotter(
    '%s/plots/%s/permutations/Comparisons' % (project, jobid),
	#defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}}
	defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['%.1f fb^{-1} (13 TeV, 25ns)' % intlumi, True]}
)

#def plot_dir():
#    print '\ncp -r /uscms/home/jdulemba/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/%s/%s/ .\n' % (analyzer, jobid)
print "File %s chosen for likelihood distributions." % args.fname

##############################################################################################

plot_types = [
    ('3J_mbpjet_vs_maxmjet', 'combination', 'm_{j}^{max} [GeV] 3 jets'),
    ('3J_mbpjet', ' b_{h}', 'm_{j_{i}+j_{j}} [Gev] 3 jets'),
    ('3J_nusolver_chi2', ' b_{l}', '#chi^{2} 3 jets'),
    ('3J_nusolver_dist', ' b_{l}', 'D_{#nu, min} 3 jets')
]

evt_types = [
    ('merged', 'Merged Jet Events'),
    ('lost', 'Lost Jet Events')
]

categories = [
    ('right', 'Correct', 'red'),
    ('wrong', 'Wrong', 'blue')
]

#set_trace()
for var, obj, xlabel in plot_types:
    for evt, leg in evt_types:
        to_draw = []
        for cat, txt_label, col in categories:
            hist = asrootpy(normfile.Get('nosys/'+var+'_%s_%s' % (evt, cat))).Clone()

            if hist.Integral() == 0:
                continue

            #set_trace()
            if hist.DIM > 1:
                plotter.plot(hist)
                hist.Draw('colz')
                hist.set_x_title(xlabel+', '+evt)
                box=plotter.make_text_box(txt_label, position='NE')
                box.Draw()
                plotter.save(var+'_%s_%s' % (evt, cat))
            else:
                plotter.set_histo_style(hist, color=col, title=txt_label+obj)
                plotter.plot(hist, drawstyle='hist')
                to_draw.append(hist)

        if len(to_draw) > 0:
            plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), drawstyle='hist', xtitle=xlabel+', '+evt)
            #box = plotter.make_text_box(leg, position='N')
            #box.Draw()
        
            plotter.save(var+'_%s' % evt)

