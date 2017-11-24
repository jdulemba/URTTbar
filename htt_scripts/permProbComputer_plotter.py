'''
permProbComputer Analyzer Plotter macro
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


analyzer = 'permProbComputer'

parser = argparse.ArgumentParser(description='Create plots using files from permProbComputer.')

jobid = jobid = os.environ['jobid']

prob_file = 'prob_ttJets_3J.root' #inputs/2017Aug24/INPUT/prob_ttJets_3J.root

##### check if correct prob_file is present
dir_f_list = []
[dir_f_list.append(f) for f in os.listdir('../inputs/%s/INPUT' % jobid)]
if not prob_file in dir_f_list:
    print "File %s isn't in the directory inputs/%s/INPUT. It has to be to run this!" % (prob_file, jobid)
    sys.exit()

myfile = root_open('../inputs/%s/INPUT/%s' % (jobid, prob_file), 'read')
normfile = views.NormalizeView(root_open('../inputs/%s/INPUT/%s' % (jobid, prob_file), 'read'))

plotter = BasePlotter(
	'plots/%s/%s' % (analyzer, jobid),
	defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}}
	#defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
)


def plot_dir():
    print '\ncp -r /uscms/home/jdulemba/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/%s/%s/ .\n' % (analyzer, jobid)


##############################################################################################

plot_types = [
    ('3J_mbpjet_vs_maxmjet', 'combination'),
    ('3J_mbpjet', ' b_{h}'),
    ('3J_nusolver_chi2', ' b_{l}'),
    ('3J_nusolver_dist', ' b_{l}')
]

evt_types = [
    ('merged', 'Merged Jet Events'),
    ('lost', 'Lost Jet Events')
]

categories = [
    ('right', 'Correct', 'red'),
    ('wrong', 'Wrong', 'blue')
]

for var, obj in plot_types:
    for evt, leg in evt_types:
        to_draw = []
        for cat, txt_label, col in categories:
#            print var+'_%s_%s' % (evt, cat)
#            set_trace()
            hist = asrootpy(normfile.Get('nosys/'+var+'_%s_%s' % (evt, cat))).Clone()
            if hist.Integral() == 0:
                continue
            if hist.DIM > 1:
                plotter.plot(hist)
                hist.Draw('colz')
                box=plotter.make_text_box(txt_label, position='NE')
                plotter.save(var+'_%s_%s' % (evt, cat))
            else:
                plotter.set_histo_style(hist, color=col, title=txt_label+obj)
                plotter.plot(hist, drawstyle='hist')
                to_draw.append(hist)

        if len(to_draw) > 0:
            plotter.overlay(to_draw, legend_def=LegendDefinition(position='NW'), drawstyle='hist')
            box = plotter.make_text_box(leg, position='NE')
            box.Draw()
        
            plotter.save(var+'_%s' % evt)

plot_dir()
