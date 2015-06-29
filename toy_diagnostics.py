#! /bin/env python

import re
import ROOT
import rootpy.plotting as plotting
import rootpy
import rootpy.io as io
from URAnalysis.PlotTools.views import EnvelopeView
from URAnalysis.Utilities.roottools import ArgSet, ArgList
from pdb import set_trace
import logging
from os.path import join
import os
import re
from URAnalysis.Utilities.struct import Struct

asrootpy = rootpy.asrootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
rootpy.log["/rootpy"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()
log = rootpy.log["/toy_diagnostics"]
rootpy.log["/toy_diagnostics"].setLevel(rootpy.log.INFO)

import URAnalysis.Utilities.prettyjson as prettyjson
from argparse import ArgumentParser
canvas = plotting.Canvas()

def mkdir(path):
   if not os.path.isdir(path):
      os.makedirs(path)


def fill_hist(vals):
   low = min(vals)
   low = low*0.8 if low > 0 else low*1.2
   hi = max(vals)
   hi *= 1.2 if hi > 0 else 0.8
   if low == hi: #0 == 0
      low, hi = -0.1, 0.1
   hist = plotting.Hist(50, low, hi, title='')
   for v in vals:
      hist.Fill(v)
   return hist

def make_hist(key, rfit_vals, out, prefix=''):
   canvas = plotting.Canvas()
   vals = [i.getVal() for i in rfit_vals]
   hist = fill_hist(vals)
   hist.Draw()
   log.debug('%s has %i entries, %.0f integral' % (key, hist.GetEntries(), hist.Integral()))
   canvas.SaveAs('%s/%s%s.png' % (out, prefix, key.replace('/','_')))
   #canvas.SaveAs('%s/%s.pdf' % (out, key.replace('/','_')))

def make_post_distributions(key, vals, out, mean_summary, sigma_summary, 
                            dist='pull', fitfunc='gaus', prefix='', prefit=None, 
                            tdir=None, skipFit=False):
   canvas = plotting.Canvas()
   values = []
   nvals  = len(vals)
   prefit_val = 1. if 'YieldSF' in key else 0.
   if prefit:
      parname = vals[0].GetName()
      prefit_val = prefit[parname].getVal() if parname in prefit else 1. #FIXME 1. assumption quite arbitrary
   if dist == 'pull':
      for val in vals:
         isValid = (val.getError() and val.hasError()) or \
            (val.getErrorLo() and val.getErrorHi() and val.hasAsymError())
         if not isValid: continue
         pull = val.getVal()-prefit_val
         err  = val.getError()
         if val.hasAsymError():
            err = val.getErrorHi() if pull <= 0 else -val.getErrorLo()
         values.append(pull/err)

      ninf = len(vals) - len(values)
      if ninf:
         log.error('%s: found %i with no paramenter error' % (key, ninf))
   elif dist == 'delta':
      prefix += dist
      values = [(i.getVal()-1) for i in vals]
   else:
      log.error('Distribution named "%s" is not supported!' % dist)
      raise ValueError('Distribution named "%s" is not supported!' % dist)

   #print key
   if filter(str.isdigit, key) != '':
      index = int(filter(str.isdigit, key))+1
   else:
      index = 1
   #print index
   
   singlekey=key.replace("YieldSF","")
   singlekey=singlekey.replace("Bin","")
   singlekey=singlekey.replace("MCStat","")
   singlekey=singlekey.strip("_1234567890")
   hist = fill_hist(values)
   function = None
   if skipFit:
      return
   if nvals > 1:
      hist.Fit(fitfunc, 'QIMES')
      function = hist.GetFunction("gaus")
      
   hist.Draw()
   canvas.SaveAs('%s/%s%s.png' % (out, prefix, key.replace('/','_')))
   tdir.WriteObject(hist, '%s%s' % (prefix, key.replace('/','_')))
   if (not hist.GetFunction("gaus") and nvals > 1):
      log.warning("Function not found for histogram %s" % name)
   else:
      mean = function.GetParameter(1) if nvals > 1 else vals[0].getVal()
      meanerror  = function.GetParError(1) if nvals > 1 else vals[0].getError()
      sigma      = function.GetParameter(2) if nvals > 1 else 0
      sigmaerror = function.GetParError(2) if nvals > 1 else 0

      mean_summary[singlekey].SetBinContent(index,mean)
      mean_summary[singlekey].SetBinError(index,meanerror)
      sigma_summary[singlekey].SetBinContent(index,sigma)
      sigma_summary[singlekey].SetBinError(index,sigmaerror)

      index = min(
         i for i in range(1, mean_summary['all'].GetNbinsX()+1) 
         if mean_summary['all'].GetBinContent(i) == 0
         )
      mean_summary['all'].SetBinContent(index,mean)
      mean_summary['all'].SetBinError(index,meanerror)
      mean_summary['all'].GetXaxis().SetBinLabel(index,key)
      sigma_summary['all'].SetBinContent(index,sigma)
      sigma_summary['all'].SetBinError(index,sigmaerror)
      sigma_summary['all'].GetXaxis().SetBinLabel(index,key)


def run_module(**kwargs):
   args = Struct(**kwargs)
   mkdir(args.out)
   canvas = plotting.Canvas()

   pars_regex = None
   if args.pars_regex:
      pars_regex = re.compile(args.pars_regex)
      
   sample_regex = None
   if args.sample_regex:
      sample_regex = re.compile(args.sample_regex)

   pars_out_regex = None
   if args.pars_out_regex:
      pars_out_regex = re.compile(args.pars_out_regex)
      
   sample_out_regex = None
   if args.sample_out_regex:
      sample_out_regex = re.compile(args.sample_out_regex)

   output_file = io.root_open('%s/output.root' % args.out, 'recreate')
   fpars_tdir = output_file.mkdir('floating_pars')
   pulls_tdir = output_file.mkdir('postfit_pulls')

   failed_fits = set()
   fit_statuses = plotting.Hist(10, -1.5, 8.5)
   with io.root_open(args.mlfit) as mlfit:
      failed_results = []
      passes_results = []
      pars = {}
      yields = {}
      first = True
      toys = [i.GetName() for i in mlfit.keys() if i.GetName().startswith('toy_')] if not args.oneshot else [None]
      log.info('examining %i toys' % len(toys))
      prefit_nuis = None
      if args.useprefit:
         prefit_nuis = ArgSet(mlfit.nuisances_prefit)

      nfailed = 0
      for toy in toys:
         toy_dir = mlfit.Get(toy) if not args.oneshot else mlfit
         keys = set([i.GetName() for i in toy_dir.GetListOfKeys()])
         if 'norm_fit_s' not in keys or 'fit_s' not in keys:
            log.error('Fit %s failed to produce output!' % toy)
            failed_fits.add(toy)
            continue
         norms = ArgSet(
            toy_dir.Get(
               'norm_fit_s'
               )
            )
         norms = [i for i in norms]

         fit_result = toy_dir.Get(
            'fit_s'
            )
         fit_pars = ArgList(fit_result.floatParsFinal())
         
         if first:
            first = False
            for i in fit_pars:
               if pars_regex and not pars_regex.match(i.GetName()): continue
               if pars_out_regex and pars_out_regex.match(i.GetName()): continue
               pars[i.GetName()] = []

            for i in norms:
               if sample_regex and not sample_regex.match(i.GetName()): continue
               if sample_out_regex and sample_out_regex.match(i.GetName()): continue
               yields[i.GetName()] = []

         fit_statuses.Fill(fit_result.status())
         fit_failed = any(i.getError() == 0 for i in fit_pars) or fit_result.status() != 0
         if fit_failed:
            log.error('Fit %s failed to converge properly. It has status %i!' % (toy, fit_result.status()))
            nfailed+=1
            failed_fits.add(toy)
            failed_results.append(fit_result)
            continue

         passes_results.append(fit_result)

         for i in norms:
            if i.GetName() in yields:
               yields[i.GetName()].append(i)
            
         for i in fit_pars:
            if i.GetName() in pars:
               pars[i.GetName()].append(i)

      if nfailed:
         log.error('There were %i fit failed!' % nfailed)
      with open('%s/info.txt' % args.out, 'w') as info:
         info.write('There were %i fit failed!\n' % nfailed)
      fit_statuses.Draw()
      canvas.SaveAs('%s/fit_status.png' % args.out)

      if not args.nopars:
         #Plots the post-fit distribution of the POI and nuisances
         out = os.path.join(args.out, 'floating_parameters')
         mkdir(out)
         for i, j in yields.iteritems():
            make_hist(i, j, out, prefix='yield_')
            
         for i, j in pars.iteritems():
            make_hist(i, j, out, prefix='par_')

      if not args.postpulls:
         #Plots the post-fit pulls (nuisance(post) - nuisance(pre))/unc(post)
         pulls_dir = os.path.join(args.out, 'postfit_pulls')
         mkdir(pulls_dir)

         ROOT.gStyle.SetOptFit(11111)
         singlenames=set()
         for name,value in pars.iteritems():
            if pars_regex and not pars_regex.match(name): continue
            if pars_out_regex and pars_out_regex.match(i): continue
            name=name.replace("YieldSF","")
            name=name.replace("Bin","")
            name=name.replace("MCStat","")
            name=name.strip("_1234567890")
            singlenames.add(name)
         
         pulls_mean_summary={}
         pulls_sigma_summary={}
         deltas_mean_summary={}
         deltas_sigma_summary={}
         for name in singlenames:
            nbins = 0
            for fullname in pars:
               if name in fullname:
                  nbins = nbins + 1
            #print name, nbins
            hist = plotting.Hist(nbins, 0.5,nbins+0.5, name = "%s_pull_mean_summary" %name)
            pulls_mean_summary[name] = hist
            hist = plotting.Hist(nbins, 0.5,nbins+0.5, name = "%s_pull_sigma_summary" %name)
            pulls_sigma_summary[name] = hist
            hist = plotting.Hist(nbins, 0.5,nbins+0.5, name = "%s_delta_mean_summary" %name)
            deltas_mean_summary[name] = hist
            hist = plotting.Hist(nbins, 0.5,nbins+0.5, name = "%s_delta_sigma_summary" %name)
            deltas_sigma_summary[name] = hist

         pulls_mean_summary[  'all'] = plotting.Hist(len(pars), 0.5, len(pars)+0.5, name = "all_pull_mean_summary"  )
         pulls_sigma_summary[ 'all'] = plotting.Hist(len(pars), 0.5, len(pars)+0.5, name = "all_pull_sigma_summary" )
         deltas_mean_summary[ 'all'] = plotting.Hist(len(pars), 0.5, len(pars)+0.5, name = "all_delta_mean_summary" )
         deltas_sigma_summary['all'] = plotting.Hist(len(pars), 0.5, len(pars)+0.5, name = "all_delta_sigma_summary")

         
         for i, j in pars.iteritems():
            make_post_distributions(i, j, pulls_dir, pulls_mean_summary, pulls_sigma_summary, prefix='pull_',
                                    dist='pull', prefit=prefit_nuis, tdir=pulls_tdir, skipFit=args.skipFit)
            make_post_distributions(i, j, pulls_dir, deltas_mean_summary, deltas_sigma_summary, prefix='delta_',
                                    dist='delta', prefit=prefit_nuis, tdir=pulls_tdir, skipFit=args.skipFit)
         
         for name,histo in pulls_mean_summary.iteritems():
            canvas = plotting.Canvas()
            histo.Draw()
            canvas.Update()
            line = ROOT.TLine(histo.GetBinLowEdge(1),0,histo.GetBinLowEdge(histo.GetNbinsX()+1),0)
            line.SetLineColor(2)
            line.Draw("same")
            canvas.Update()
            canvas.SaveAs('%s/%s.png' % (pulls_dir,histo.GetName()))
            canvas.SaveAs('%s/%s.pdf' % (pulls_dir,histo.GetName()))
            pulls_tdir.WriteObject(histo, histo.GetName())
         for name,histo in pulls_sigma_summary.iteritems():
            canvas = plotting.Canvas()
            histo.Draw()
            canvas.Update()
            line = ROOT.TLine(histo.GetBinLowEdge(1),1,histo.GetBinLowEdge(histo.GetNbinsX()+1),1)
            line.SetLineColor(2)
            line.Draw("same")
            canvas.Update()
            canvas.SaveAs('%s/%s.png' % (pulls_dir,histo.GetName()))
            canvas.SaveAs('%s/%s.pdf' % (pulls_dir,histo.GetName()))
            pulls_tdir.WriteObject(histo, histo.GetName())
         
         for name,histo in deltas_mean_summary.iteritems():
            canvas = plotting.Canvas()
            histo.Draw()
            canvas.Update()
            line = ROOT.TLine(histo.GetBinLowEdge(1),0,histo.GetBinLowEdge(histo.GetNbinsX()+1),0)
            line.SetLineColor(2)
            line.Draw("same")
            canvas.Update()
            canvas.SaveAs('%s/%s.png' % (pulls_dir,histo.GetName()))
            canvas.SaveAs('%s/%s.pdf' % (pulls_dir,histo.GetName()))
            pulls_tdir.WriteObject(histo, histo.GetName())
         for name,histo in deltas_sigma_summary.iteritems():
            histo.Draw()
            canvas.Update()
            #line = ROOT.TLine(histo.GetBinLowEdge(1),1,histo.GetBinLowEdge(histo.GetNbinsX()+1),1)
            #line.Draw("same")
            canvas.Update()
            canvas.SaveAs('%s/%s.png' % (pulls_dir,histo.GetName()))
            canvas.SaveAs('%s/%s.pdf' % (pulls_dir,histo.GetName()))
            pulls_tdir.WriteObject(histo, histo.GetName())


   if not args.noshapes:
      #Overlays the prefit values of the different shapes with the envelope of 
      #what is fitted by the toys
      out = os.path.join(args.out, 'shapes')
      mkdir(out)
      biased_shapes={}
      if args.biasFile:
         with io.root_open(args.biasFile) as biased:
            biased_dir= biased.prefit \
               if hasattr(biased, 'prefit') else \
               None
            ROOT.TH1.AddDirectory(False)
            for key in biased_dir.keys():
               biased_shapes[key.name] = asrootpy(key.ReadObj().Clone())

      with io.root_open(args.harvested) as harvest:
         has_prefit = hasattr(harvest, 'prefit')
         prefit = harvest.prefit if has_prefit else None
         toys = EnvelopeView(
            *[harvest.get(i.GetName()).get(args.variable) 
              for i in harvest.keys() 
              if i.GetName().startswith('toy_')
              and (i.GetName() not in failed_fits) ]
             )
         #shapes = [i.GetName() for i in prefit.keys()] #FIXME! should not depend on prefit!
         first_toy = [i.GetName() for i in harvest.keys() if i.GetName().startswith('toy_')][0]
         not_shapes = set('correlation_matrix')
         shapes = [i.GetName() for i in harvest.get(first_toy).get(args.variable).keys() if i.GetName() not in not_shapes]

         for shape in shapes:
            legend = plotting.Legend(5, rightmargin=0.07, topmargin=0.05, leftmargin=0.45)
            legend.SetBorderSize(0)
            biased_shape = biased_shapes.get(shape, None)
            
            if has_prefit:
               pre_shape = prefit.Get(shape)
               pre_shape.title = 'input shape'
               pre_shape.legendstyle = 'p'
               pre_shape.drawstyle = 'p'
               if biased_shape:
                  pre_shape.legendstyle = 'l'
                  pre_shape.drawstyle = 'hist'
                  pre_shape.linecolor = 'blue'
                  pre_shape.fillstyle = 0
            toy_shape = toys.Get(shape)
            
            toy_shape.Draw()
            if has_prefit:
               pre_shape.Draw('same')
               
            if biased_shape:
               biased_shape.title = 'true shape'
               biased_shape.legendstyle = 'p'
               biased_shape.inlegend = True               
               biased_shape.drawstyle = 'p'
               biased_shape.Draw('same')

            legend.AddEntry(toy_shape.two_sigma)
            legend.AddEntry(toy_shape.one_sigma)
            legend.AddEntry(toy_shape.median)
            if has_prefit:
               legend.AddEntry(pre_shape)
            if biased_shape:
               legend.AddEntry(biased_shape)
            legend.Draw()

            canvas.Update()
            canvas.SaveAs('%s/%s.png' % (out, shape))
            canvas.SaveAs('%s/%s.pdf' % (out, shape))
            with open(os.path.join(out, '%s.json' % shape), 'w') as jfile:
               jfile.write(toy_shape.json())

   output_file.Close()


if __name__ == '__main__':
   parser = ArgumentParser()
   parser.add_argument('variable')
   parser.add_argument('harvested')
   parser.add_argument('mlfit')
   parser.add_argument(
      '-o', dest='out', type=str, default='./',
      help='output directory'
      )
   parser.add_argument('--noshapes', action='store_true')
   parser.add_argument('--nopars', action='store_true')
   parser.add_argument('--no-post-pulls', dest='postpulls', action='store_true')
   parser.add_argument('--use-prefit', dest='useprefit', action='store_true')
   parser.add_argument('--filter-in-pars', dest='pars_regex', type=str)
   parser.add_argument('--filter-in-sample', dest='sample_regex', type=str)
   parser.add_argument('--filter-out-pars', dest='pars_out_regex', type=str)
   parser.add_argument('--filter-out-sample', dest='sample_out_regex', type=str)
   parser.add_argument('--oneshot', action='store_true', help='to use to store data/asimov runs')
   parser.add_argument('--skipFit', action='store_true', help='skips pulls fit, when really nothing works')
   parser.add_argument('--biasFile', help='fit harvest of the biased toys with biased shapes.\n'
                       'Used to get the values of the REAL prefit shapes in case of bias tests')

   args = parser.parse_args()
   run_module(**dict(args._get_kwargs()))

