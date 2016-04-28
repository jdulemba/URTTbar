#! /bin/env python

import ROOT
ROOT.gROOT.SetBatch(True)
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
ROOT.gStyle.SetOptFit(11111)

import URAnalysis.Utilities.prettyjson as prettyjson
from argparse import ArgumentParser
import uuid
from URAnalysis.Utilities.struct import Struct

def run_module(**kwargs):
   args = Struct(**kwargs)
   if not args.name:
      args.name = args.var
   results = [prettyjson.loads(open(i).read())[-1] for i in args.results]
   #set_trace()
   results.sort(key=lambda x: x['median'])

   nevts_graph = plotting.Graph(len(results))
   upbound_graph = plotting.Graph(len(results))
   max_unc = 0.
   bound_range = args.vrange #results[-1]['up_edge'] - results[0]['up_edge']
   step = results[1]['up_edge'] - results[0]['up_edge']
   #bound_range += step
   bound_min = results[0]['up_edge'] - step

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
      os.path.join(args.outdir, 'nevts_%s.png' % args.name)
      )

   tf1 = Plotter.parse_formula(
      'scale / (x - shift) + offset', 
      'scale[1,0,10000],shift[%.2f,%.2f,%.2f],offset[0, 0, 1]' % (bound_min, bound_min-2*step, bound_min+step)
      )
   # ROOT.TF1('ret', '[0]/(x - [1])', -3, 1000)
   tf1.SetRange(0, 1000)
   tf1.SetLineColor(ROOT.EColor.kAzure)
   tf1.SetLineWidth(3)
   result = upbound_graph.Fit(tf1, 'MES') #WL

   scale = tf1.GetParameter('scale') 
   shift = tf1.GetParameter('shift') 
   offset = tf1.GetParameter('offset')

   upbound_graph.Draw('APL')
   upbound_graph.GetXaxis().SetTitle('upper bin edge')
   upbound_graph.GetYaxis().SetTitle('relative fit uncertainty')
   if args.fullrange:
      upbound_graph.GetYaxis().SetRangeUser(offset, max_unc*1.2)
   upbound_graph.GetXaxis().SetRangeUser(shift, (shift+bound_range)*1.2)


   delta = lambda x, y: ((x-shift)/bound_range)**2+((y-offset)/(max_unc-offset))
   points = ((bound_min+step)+i*(bound_range/100.) for i in xrange(100))
   math_best_x = min(points, key=lambda x: delta(x, tf1.Eval(x)))
   math_best_y = tf1.Eval(math_best_x)

   ## math_best_x = math.sqrt(scale*bound_range*(max_unc-offset))+shift
   ## math_best_y = tf1.Eval(math_best_x) #max_unc*math.sqrt(scale)+offset
   best = min(results, key=lambda x: abs(x['up_edge']-math_best_x))

   upbound_best = plotting.Graph(1)
   upbound_best.SetPoint(0, best['up_edge'], best["one_sigma"]["relative"])
   upbound_best.markerstyle = 29
   upbound_best.markersize  = 3
   upbound_best.markercolor = 2
   upbound_best.Draw('P same')

   print math_best_x, math_best_y
   print best['up_edge'], best["one_sigma"]["relative"]

   math_best = plotting.Graph(1)
   math_best.SetPoint(0, math_best_x, math_best_y)
   math_best.markerstyle = 29
   math_best.markersize  = 3
   math_best.markercolor = ROOT.EColor.kAzure
   math_best.Draw('P same')

   canvas.SaveAs(
      os.path.join(args.outdir, 'upbound_%s.png' % args.name)
      )
   json = {
      'best' : best['up_edge'],
      'unc'  : best["one_sigma"]["relative"]
      }
   with open(os.path.join(args.outdir, '%s.json' % args.name), 'w') as jfile:
      jfile.write(prettyjson.dumps(json))


if __name__ == '__main__':
   parser = ArgumentParser()
   parser.add_argument('var')
   parser.add_argument('outdir')
   parser.add_argument('vrange', type=float)
   parser.add_argument('results', type=str, nargs='+')
   parser.add_argument('--name', help='name to give to the output files')
   parser.add_argument('--fullrange', action='store_false', help='infer the Y-axis range from the points, rather than from the fit')
   
   args = parser.parse_args()
   run_module(**dict(args._get_kwargs()))
