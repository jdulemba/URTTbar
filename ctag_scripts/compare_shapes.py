#! /bin/env python

import rootpy as rpy
import rootpy.plotting as plotting
import rootpy.io as io
import sys
import glob
import os
from pdb import set_trace
rpy.log["/"].setLevel(rpy.log.INFO)
import ROOT
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
from argparse import ArgumentParser
from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
from styles import styles
from fnmatch import fnmatch

parser = ArgumentParser()
parser.add_argument('shapes')
parser.add_argument('--filter', default='')
args = parser.parse_args()

common_style = {
   'legendstyle' : 'l',
   'drawstyle' : 'hist',
   'fillstyle': 'hollow',
   'linewidth' : 3,
}
getstyle = lambda x: {'linecolor' : styles[x]['linecolor'], 'title' : styles[x]['name']}
plotter = BasePlotter(
   os.path.join(
      os.path.dirname(args.shapes),
      'shapes'
      ),
   defaults={'clone' : False},
   styles = {
      'right_whad' : {
         'linecolor' : '#6666b3',
         'title' : "tt, right W_{h}",
         },
      'qcd' : getstyle('QCD*'),
      'single_top' : getstyle('single*'),
      'wrong_whad' : {
         'linecolor' : '#ab5555',
         'title' : "tt, wrong W_{h}",
         },
      'nonsemi_tt' : {
         'linecolor' : '#668db3',
         'title' : "tt, non semi-lep.",
         },      
      'vjets' : getstyle('[WZ]Jets*'),
      }
   )

shapes = io.root_open(args.shapes)
categories = [i.name for i in shapes.keys()]
to_plot = {'wrong_whad', 'nonsemi_tt', 'single_top', 'qcd', 'right_whad', 'vjets'}

sys_colors = {
   'JES' : ('JES', '#0000ff'),
   'JER' : ('JER', '#ff0000'),
   'HScale' : ('H. scale', '#008000'),
   'PDF' : ('PDF', '#660066'),
   'MTOP' : ('M top', '#ff9b02'),
}

sys_lstyles = ['solid', 'dashed']

for name in categories:  
   category = shapes.Get(name)
   all_keys = [i.name for i in category.keys()]
   samples = [i for i in all_keys if i in to_plot]
   histos = [category.Get(i).Clone() for i in samples]
   linecolors = [i.linecolor for i in histos]
   for i in histos:
      i.xaxis.title = '#lambda_{M}'
      i.yaxis.title = 'Events'
      integral = i.Integral()
      i.Scale(1./integral if integral else 1.)
   
   plotter.overlay(histos, LegendDefinition(position='NE'), **common_style)
   plotter.save(name)

   systematics = [i for i in all_keys if i.endswith('Down') and '_bin_' not in i and not i.endswith('_MCStatDown')]
   if args.filter:
      systematics = [i for i in systematics if fnmatch(i, args.filter)]

   for sample in samples:
      base = category.Get(sample).Clone()
      sys_dev = [i for i in systematics if i.startswith(sample)]
      if len(sys_dev) > len(sys_colors):
         raise RuntimeError('you have more systematics than colors to assign!')
      to_draw = [base]
      to_legend = [base]
      lstyles = ['solid']+sys_lstyles*len(sys_dev)      
      lcolors = ['black']
      for dev in sys_dev:
         sys_name = dev[len(sample):-len('Down')].strip('_')
         title, col = sys_colors[sys_name]
         lcolors.extend([col, col])
         hd = category.Get(dev)
         hd.title = title
         hu = category.Get(dev.replace('Down', 'Up'))
         hu.title = title
         to_draw.extend([hu, hd])
         to_legend.append(hu)

      legend = LegendDefinition(position='NE', entries=to_legend)
      plotter.overlay(to_draw, legend, linestyle=lstyles, linecolor=lcolors, **common_style)
      #plotter.add_legend(to_legend, False)
      plotter.save('%s_%s' % (name, sample))
      
   
