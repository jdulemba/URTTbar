#! /bin/env python

import ROOT
import rootpy.plotting as plotting
import rootpy
import rootpy.io as io
from URAnalysis.Utilities.roottools import ArgSet
from pdb import set_trace
import logging

asrootpy = rootpy.asrootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)

import URAnalysis.Utilities.prettyjson as prettyjson
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('fitresult')
parser.add_argument('binning')
parser.add_argument(
   '-o', dest='out', type=str, default='out.root',
   help='output')
args = parser.parse_args()

binning = {}
with open(args.binning) as bins:
   binning = prettyjson.loads(bins.read())

for info in binning.itervalues():
   edges = set(
      i['low_edge'] for i in info.itervalues()
      )
   edges.update(
      set(
         i['up_edge'] for i in info.itervalues()
         )
      )
   edges = sorted(list(edges))
   info['edges'] = edges

with io.root_open(args.fitresult) as results:
   norms = ArgSet(results.Get('norm_fit_s'))
   norms = [i for i in norms]
   fit_result = results.Get('fit_s')
   
   with io.root_open(args.out, 'recreate') as output:
      hcorr = fit_result.correlationHist()
      hcorr.Write()
      for obs, info in binning.iteritems():
         var_dir = output.mkdir(obs)
         var_dir.cd()
         hists = {}
         for norm in norms:
            category, sample = tuple(norm.GetName().split('/'))
            if category not in info: continue
            if sample not in hists:
               hists[sample] = plotting.Hist(
                  info['edges'],
                  name = sample
                  )
            hists[sample].SetBinContent(
               info[category]['idx']+1, norm.getVal()
               )
            hists[sample].SetBinError(
               info[category]['idx']+1, norm.getError()
               )

         for h in hists.itervalues():
            logging.debug( h.Write() )
         #var_dir.Write()
