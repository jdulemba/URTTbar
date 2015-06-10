#! /bin/env python

import ROOT
import rootpy.plotting as plotting
import rootpy
import rootpy.io as io
from URAnalysis.Utilities.roottools import ArgSet, ArgList
from pdb import set_trace
import logging
import os
from URAnalysis.PlotTools.Plotter import Plotter
import math

asrootpy = rootpy.asrootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptFit(11111)

import URAnalysis.Utilities.prettyjson as prettyjson
from argparse import ArgumentParser
import uuid

parser = ArgumentParser()
parser.add_argument('var')
parser.add_argument('outdir')
parser.add_argument('results', type=str, nargs='+')
                    
args = parser.parse_args()

results = [prettyjson.loads(open(i).read())[0] for i in args.results]
#set_trace()
results.sort(key=lambda x: x['median'])

nevts_graph = plotting.Graph(len(results))
upbound_graph = plotting.Graph(len(results))
max_unc = 0.
bound_range = results[-1]['up_edge'] - results[0]['up_edge']

for idx, info in enumerate(results):
   nevts_graph.SetPoint(idx, info['median'], info["one_sigma"]["relative"])
   upbound_graph.SetPoint(idx, info['up_edge'], info["one_sigma"]["relative"])
   if info["one_sigma"]["relative"] > max_unc:
      max_unc = info["one_sigma"]["relative"]

canvas = plotting.Canvas()
nevts_graph.Draw('APL')
nevts_graph.GetXaxis().SetTitle('average number of events')
nevts_graph.GetYaxis().SetTitle('relative fit uncertainty')
canvas.SaveAs(
   os.path.join(args.outdir, 'nevts_%s.png' % args.var)
   )

upbound_graph.Draw('APL')
upbound_graph.GetXaxis().SetTitle('upper bin edge')
upbound_graph.GetYaxis().SetTitle('relative fit uncertainty')

tf1 = Plotter.parse_formula('scale / (x - shift)', 'scale[1,0,10000],shift[0,0,1000]')
# ROOT.TF1('ret', '[0]/(x - [1])', -3, 1000)
tf1.SetRange(0, 1000)
tf1.SetLineColor(ROOT.EColor.kAzure)
tf1.SetLineWidth(3)
result = upbound_graph.Fit(tf1, 'MES') #WL

scale = tf1.GetParameter('scale') 
shift = tf1.GetParameter('shift') 

distance = float('inf')
min_idx = -1
distances = []
for idx, info in enumerate(results):
   new_dist = math.sqrt(
      ((info['up_edge'] - shift)/bound_range)**2 + 
      (info["one_sigma"]["relative"]/max_unc)**2
      )
   distances.append(new_dist)
   #print idx, new_dist, distance
   if new_dist < distance:
      distance = new_dist
      min_idx  = idx

#max_unc = float('inf')
#bound_range = results[-1]['up_edge'] - results[0]['up_edge']
#print results[min_idx]['up_edge'], results[min_idx]["one_sigma"]["relative"]

upbound_best = plotting.Graph(1)
upbound_best.SetPoint(0, results[min_idx]['up_edge'], results[min_idx]["one_sigma"]["relative"])
upbound_best.markerstyle = 29
upbound_best.markersize  = 3
upbound_best.markercolor = 2
upbound_best.Draw('P same')

canvas.SaveAs(
   os.path.join(args.outdir, 'upbound_%s.png' % args.var)
   )
