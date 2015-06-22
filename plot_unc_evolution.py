#! /bin/env python                                                                                                                           

import URAnalysis.Utilities.prettyjson as prettyjson
from argparse import ArgumentParser
from rootpy.plotting import Canvas, Hist, Legend, Pad
import os
from pdb import set_trace
import ROOT
from math import ceil
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

parser = ArgumentParser()
parser.add_argument('bin_directory', help='where all the trial bins are stored')
parser.add_argument('binning', nargs='+', help='input json files')
parser.add_argument('-o', dest='out', help='output filename', default='out.png')

args = parser.parse_args()
npoints = len(args.binning) -1
input_dirs = ['_'.join(args.binning[:idx]) for idx in range(2, len(args.binning)+1)]
inputs = [os.path.join(args.bin_directory, i, 'toys/shapes/tt_right.json') for i in input_dirs]

jsons = []
for i in inputs:
   jsons.append(
      prettyjson.loads(open(i).read()) 
      )
jsons.sort(key=lambda x: len(x))
max_bins = len(jsons[-1])

colors = ['#FFCC66', '#2aa198', '#9999CC', '#5555ab', 
          '#aaaad5', '#ab5555', '#d5aaaa', '#668db3',
          '#6699ff', '#ff8066', '#a12a33', '#2aa15d',
          '#a15d2a', '#000000', '#a206f4', '#06f4a2',
          '#f406cf', '#cff406']
hists = [Hist(max_bins, 0.5, max_bins+0.5, title='Bin %i' % (i+1)) for i in range(max_bins)]
for hist, color in zip(hists, colors):
   hist.markercolor = color
   hist.fillcolor = color
   hist.drawstyle = 'P'
   hist.inlegend = True
   hist.legendstyle = 'p'

for idx, json in enumerate(jsons):
   for hist, fit_bin in zip(hists, json):
      hist[idx+1].value = fit_bin['one_sigma']['relative']

canvas = Canvas(800, 800)
pad = Pad(0, 0., 1., 0.33)
pad.SetPad(0, 0.33, 1., 1.)
pad.Draw()
lower_pad = Pad(0, 0., 1., 0.33)
lower_pad.Draw()

nhists = len(hists)
nlegends = int(ceil(float(nhists) / 5))
nhist_for_leg = int(ceil(float(nhists)/nlegends))
legends = []
left_margin = 0.1
right_margin = 0.8
for _ in range(nlegends-1):
   legends.append(
      Legend(
         nhist_for_leg,
         leftmargin=left_margin,
         rightmargin=right_margin
         )
      )
   left_margin += 0.1
   right_margin -= 0.1

legends.append(
   Legend(
      nhists - nhist_for_leg*(nlegends-1),
      leftmargin=left_margin,
      rightmargin=right_margin
      )
   )

for idx, i in enumerate(hists):
   legend_id = idx / nhist_for_leg
   try:
      legends[legend_id].AddEntry(i)
   except:
      set_trace()

maxval = max(max(i.value for i in hist) for hist in hists)
hists[0].GetYaxis().SetRangeUser(0, maxval*1.2)
pad.cd()
hists[0].Draw()
for i in hists[1:]:
   i.Draw('same')
pad.Draw()

lower_pad.cd()
for legend in legends:
   legend.Draw()
lower_pad.Draw()

canvas.Update()
canvas.SaveAs(args.out)
