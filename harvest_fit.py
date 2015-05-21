#! /bin/env python

import ROOT
import rootpy.plotting as plotting
import rootpy
import rootpy.io as io
from URAnalysis.Utilities.roottools import ArgSet, ArgList
from pdb import set_trace
import logging

asrootpy = rootpy.asrootpy
rootpy.log["/"].setLevel(rootpy.log.DEBUG)
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
   pars = ArgList(fit_result.floatParsFinal())
   pars = dict((i.GetName(), i) for i in pars)
   
   with io.root_open(args.out, 'recreate') as output:
      hcorr = asrootpy(fit_result.correlationHist())
      fit_pars = set(
         hcorr.xaxis.GetBinLabel(i) for i in range(1,hcorr.xaxis.GetNbins()+1)
         )
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
            idx = info[category]['idx']+1
            hists[sample].SetBinContent(
               idx, norm.getVal()
               )
            hists[sample].SetBinError(
               idx, norm.getError()
               )
            fit_par_name = '%sYieldSF_%s' % (category, sample)
            if fit_par_name in fit_pars:
               fit_par = pars[fit_par_name]
               logging.debug('%s' % fit_par_name)
               logging.debug('norm %.2f +/- %.2f' % (norm.getVal(), norm.getError()))
               logging.debug('SF   %.4f +/- %.4f' % (fit_par.getVal(), fit_par.getError()))
               logging.debug('norm unc. %.6f, SF unc. %.6f' % (
                     norm.getError()/norm.getVal(),
                     fit_par.getError()/fit_par.getVal()
                     )
                             )
               # assert(norm.getError()/norm.getVal() == fit_par.getError()/fit_par.getVal())
               logging.debug(
                  'Assigning label %s to bin %i for %s/%s' % (fit_par_name, idx, category, sample)
                  )
               hists[sample].xaxis.SetBinLabel(idx, fit_par_name)

         for h in hists.itervalues():
            logging.debug( h.Write() )
         #var_dir.Write()
