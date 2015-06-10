#! /bin/env python

import ROOT
import rootpy.plotting as plotting
import rootpy
import rootpy.io as io
from URAnalysis.Utilities.roottools import ArgSet, ArgList
from pdb import set_trace
import logging
import os

asrootpy = rootpy.asrootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

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

for idx, info in enumerate(results):
   nevts_graph.SetPoint(idx, info['median'], info["one_sigma"]["relative"])
   upbound_graph.SetPoint(idx, info['up_edge'], info["one_sigma"]["relative"])

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
canvas.SaveAs(
   os.path.join(args.outdir, 'upbound_%s.png' % args.var)
   )
