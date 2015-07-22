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

asrootpy = rootpy.asrootpy
rootpy.log["/"].setLevel(rootpy.log.ERROR)
rootpy.log["/rootpy"].setLevel(rootpy.log.ERROR)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

def run_module(**kwargs):
   args = Struct(**kwargs)
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
      
      is_prefit_done = False
      with io.root_open(args.out, 'recreate') as output:
         for dirname in dirs:
            input_dir = results.Get(dirname) if dirname else results
            keys = set([i.GetName() for i in input_dir.GetListOfKeys()])
            if 'norm_fit_s' not in keys or 'fit_s' not in keys:
               continue
            norms = ArgSet(
               input_dir.Get(
                  'norm_fit_s'
                  )
               )
            norms = [i for i in norms]

            fit_result = input_dir.Get(
               'fit_s'
               )
            pars = ArgList(fit_result.floatParsFinal())
            pars = dict((i.GetName(), i) for i in pars)
      
            tdir = output
            if dirname: 
               tdir = output.mkdir(dirname)
               tdir.cd()
            hcorr = asrootpy(fit_result.correlationHist())
            fit_pars = set(
               hcorr.xaxis.GetBinLabel(i) for i in range(1,hcorr.xaxis.GetNbins()+1)
               )
            hcorr.Write()
            for obs, info in binning.iteritems():
               var_dir = tdir.mkdir(obs)
               var_dir.cd()
               hists = {}
               hists_prefit = {}
               for norm in norms:
                  category, sample = tuple(norm.GetName().split('/'))
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
                  hists[sample].SetBinContent(
                     idx, norm.getVal()
                     )
                  hists[sample].SetBinError(
                     idx, norm.getError()
                     )
                  if not is_prefit_done:
                     hists_prefit[sample].SetBinContent(
                        idx, prefit_norms[norm.GetName()].getVal()
                        )
                     hists_prefit[sample].SetBinError(
                        idx, prefit_norms[norm.GetName()].getError()
                        )                  

                  fit_par_name = '%sYieldSF_%s' % (category, sample)
                  if fit_par_name in fit_pars:
                     fit_par = pars[fit_par_name]
                     logging.debug('%s' % fit_par_name)
                     logging.debug('norm %.2f +/- %.2f' % (norm.getVal(), norm.getError()))
                     logging.debug('SF   %.4f +/- %.4f' % (fit_par.getVal(), fit_par.getError()))
                     rel_y_err = norm.getError()/norm.getVal() if norm.getVal() else -1
                     rel_f_err = fit_par.getError()/fit_par.getVal() if fit_par.getVal() else -1
                     ## logging.debug(
                     ##    'norm unc. %.6f, SF unc. %.6f, delta %.3f' % (
                     ##       rel_y_err,
                     ##       rel_f_err,
                     ##       abs(rel_y_err-rel_f_err)*2/(rel_y_err+rel_f_err)
                     ##       )
                     ##    )
                     # assert(norm.getError()/norm.getVal() == fit_par.getError()/fit_par.getVal())
                     logging.debug(
                        'Assigning label %s to bin %i for %s/%s' % (fit_par_name, idx, category, sample)
                        )
                     hists[sample].xaxis.SetBinLabel(idx, fit_par_name)

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
