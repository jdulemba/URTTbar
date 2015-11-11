#! /bin/env python

import ROOT
import rootpy.plotting as plotting
import rootpy
import rootpy.io as io
from URAnalysis.Utilities.roottools import ArgSet, ArgList
from pdb import set_trace
import logging
from os.path import join
import URAnalysis.Utilities.prettyjson as prettyjson
from argparse import ArgumentParser
import uuid
from URAnalysis.Utilities.struct import Struct
import re
import copy

asrootpy = rootpy.asrootpy
rootpy.log["/"].setLevel(rootpy.log.ERROR)
rootpy.log["/rootpy"].setLevel(rootpy.log.ERROR)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

def run_module(**kwargs):
   args = Struct(**kwargs)
   redundant_binning = {}
   with open(args.binning) as bins:
      redundant_binning = prettyjson.loads(bins.read())

   #group binning to merge jet categories
   grouping = re.compile('^(?P<base_category>[A-Za-z0-9]+)_\d+Jets$')
   binning = {}
   for var, categories in redundant_binning.iteritems():
      if var not in binning: binning[var] = {}
      for category, bin_info in categories.iteritems():
         m = grouping.match(category)
         if not m:
            raise ValueError('Category name %s did not match the regex!' % category)
         base = m.group('base_category')
         if base not in binning[var]:
            binning[var][base] = copy.deepcopy(bin_info)
         else:
            #make sure that all jet categories have the same bin edges
            assert(binning[var][base] == bin_info)

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

   prefit_norms = {}
   with io.root_open(args.input_shape) as shapes:
      for key in shapes.keys():
         obj = key.ReadObj()
         if not obj.InheritsFrom('TDirectory'): continue
         
         for hkey in obj.GetListOfKeys():
            hist = hkey.ReadObj()
            if not hist.InheritsFrom('TH1'): continue
            
            err = ROOT.Double()
            integral = hist.IntegralAndError(
               1,
               hist.GetNbinsX(),
               err
               )
            val_id = uuid.uuid4().hex
            val =  ROOT.RooRealVar(val_id, val_id, integral)
            val.setError(err)
            prefit_norms['%s/%s' % (obj.GetName(), hist.GetName())] = val
               

   with io.root_open(args.fitresult) as results:
      dirs = ['']
      if args.toys:
         dirs = [i.GetName() for i in results.GetListOfKeys() if i.GetName().startswith('toy_')]
      
      with io.root_open(args.out, 'recreate') as output:
         is_prefit_done = False
         for dirname in dirs:
            input_dir = results.Get(dirname) if dirname else results
            if not hasattr(input_dir, 'fit_s'):
               continue

            fit_result = input_dir.fit_s
            pars = asrootpy(fit_result.floatParsFinal())
            prefit_pars = asrootpy(fit_result.floatParsInit())
            tdir = output
            if dirname: 
               tdir = output.mkdir(dirname)
               tdir.cd()
            hcorr = asrootpy(fit_result.correlationHist())
            par_names = set([i.name for i in pars])
            yield_par_names = filter(lambda x: '_FullYield_' in x, par_names)
            hcorr.Write()
            for observable, info in binning.iteritems():
               var_dir = tdir.mkdir(observable)
               var_dir.cd()
               hists = {}
               hists_prefit = {}
               for rvar_name in yield_par_names:
                  category, sample = tuple(rvar_name.split('_FullYield_'))
                  if category not in info: continue
                  if sample not in hists:
                     hists[sample] = plotting.Hist(
                        info['edges'],
                        name = sample
                        )
                     if not is_prefit_done:
                        hists_prefit[sample] = plotting.Hist(
                           info['edges'],
                           name = sample
                           )
    
                  idx = info[category]['idx']+1
                  hists[sample][idx].value = pars[rvar_name].value
                  error = pars[rvar_name].error
                  hists[sample][idx].error = max(abs(i) for i in error) if isinstance(error, tuple) else error #get max of asym error

                  if not is_prefit_done:
                     hists_prefit[sample][idx].value = prefit_pars[rvar_name].value
                  ## Pre-fit floating parameters have no uncertainties
                  ## hists_prefit[sample][idx].error = max(prefit_pars[rvar_name].error)
                  logging.debug(
                     'Assigning label %s to bin %i for %s/%s' % (rvar_name, idx, category, sample)
                     )
                  hists[sample].xaxis.SetBinLabel(idx, rvar_name)

               for h in hists.itervalues():
                  logging.debug( h.Write() )

               if not is_prefit_done:
                  is_prefit_done = True
                  output.mkdir('prefit').cd()
                  for h in hists_prefit.itervalues():
                     logging.debug( h.Write() )
               #var_dir.Write()

if __name__ == '__main__':
   parser = ArgumentParser()
   parser.add_argument('fitresult') #, nargs='+')
   parser.add_argument('binning')
   parser.add_argument('input_shape')
   parser.add_argument(
      '-o', dest='out', type=str, default='out.root',
      help='output')
   parser.add_argument(
      '-t','--toys', dest='toys', action='store_true',
      help='configure script to access toys'
      )
   
   args = parser.parse_args()
   run_module(**dict(args._get_kwargs()))
