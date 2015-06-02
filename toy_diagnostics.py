#! /bin/env python

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

asrootpy = rootpy.asrootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()
log = rootpy.log["/toy_diagnostics"]

import URAnalysis.Utilities.prettyjson as prettyjson
from argparse import ArgumentParser

def mkdir(path):
   if not os.path.isdir(path):
      os.makedirs(path)
   

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

args = parser.parse_args()
mkdir(args.out)

if not args.noshapes:
   #Overlays the prefit values of the different shapes with the envelope of 
   #what is fitted by the toys
   out = os.path.join(args.out, 'shapes')
   mkdir(out)
   with io.root_open(args.harvested) as harvest:
      has_prefit = hasattr(harvest, 'prefit')
      prefit = harvest.prefit if has_prefit else None
      toys = EnvelopeView(
         *[harvest.get(i.GetName()).get(args.variable) 
           for i in harvest.keys() 
           if i.GetName().startswith('toy_')]
          )
      shapes = [i.GetName() for i in prefit.keys()]
      
      for shape in shapes:
         canvas = plotting.Canvas()
         legend = plotting.Legend(4, rightmargin=0.07, topmargin=0.05, leftmargin=0.45)
         legend.SetBorderSize(0)
         if has_prefit:
            pre_shape = prefit.Get(shape)
            pre_shape.title = 'input shape'
            pre_shape.legendstyle = 'p'
         toy_shape = toys.Get(shape)
         
         toy_shape.Draw()
         if has_prefit:
            pre_shape.Draw('P same')

         legend.AddEntry(toy_shape.two_sigma)
         legend.AddEntry(toy_shape.one_sigma)
         legend.AddEntry(toy_shape.median)
         if has_prefit:
            legend.AddEntry(pre_shape)
         legend.Draw()

         canvas.Update()
         canvas.SaveAs('%s/%s.png' % (out, shape))
         canvas.SaveAs('%s/%s.pdf' % (out, shape))


def fill_hist(vals):
   low = min(vals)
   low = low*0.8 if low > 0 else low*1.2
   hi = max(vals)*1.2
   hist = plotting.Hist(50, low, hi, title='')
   for v in vals:
      hist.Fill(v)
   return hist

def make_hist(key, rfit_vals, out, prefix=''):
   vals = [i.getVal() for i in rfit_vals]
   hist = fill_hist(vals)
   hist.Draw()
   log.debug('%s has %i entries, %.0f integral' % (key, hist.GetEntries(), hist.Integral()))
   canvas.SaveAs('%s/%s%s.png' % (out, prefix, key.replace('/','_')))
   #canvas.SaveAs('%s/%s.pdf' % (out, key.replace('/','_')))

def make_post_pulls(key, vals, out, prefix=''):
   canvas = plotting.Canvas()
   pulls = []
   if 'YieldSF' in key:
      pulls = [(i.getVal()-1)/i.getError() for i in vals]
   else:
      pulls = [i.getVal()/i.getError() for i in vals]

   hist = fill_hist(pulls)
   hist.Fit('gaus', 'IMES')
   hist.Draw()
   canvas.SaveAs('%s/%s%s.png' % (out, prefix, key.replace('/','_')))

if not args.nopars and not args.postpulls:
   #Plots the post-fit distribution of the POI and nuisances
   out = os.path.join(args.out, 'floating_parameters')
   mkdir(out)

   #Plots the post-fit pulls (nuisance(post) - nuisance(pre))/unc(post)
   pulls_dir = os.path.join(args.out, 'postfit_pulls')
   mkdir(pulls_dir)

   with io.root_open(args.mlfit) as mlfit:
      pars = {}
      yields = {}
      first = True
      toys = [i.GetName() for i in mlfit.keys() if i.GetName().startswith('toy_')]
      log.info('examining %i toys' % len(toys))
      for toy in toys:
         norms = ArgSet(
            mlfit.Get(
               join(
                  toy,
                  'norm_fit_s'
                  )
               )
            )
         norms = [i for i in norms]

         fit_result = mlfit.Get(
            join(
               toy,
               'fit_s'
               )
            )
         fit_pars = ArgList(fit_result.floatParsFinal())
         
         if first:
            first = False
            pars.update( 
               dict(
                  (i.GetName(), []) for i in fit_pars
                  )
               )
            yields.update(
               dict(
                  (i.GetName(), []) for i in norms
                  )
               )

         for i in norms:
            yields[i.GetName()].append(i)
            
         for i in fit_pars:
            pars[i.GetName()].append(i)

      if not args.nopars:
         for i, j in yields.iteritems():
            make_hist(i, j, out, prefix='yield_')
            
         for i, j in pars.iteritems():
            make_hist(i, j, out, prefix='par_')

      if not args.postpulls:
         ROOT.gStyle.SetOptFit(11111)
         for i, j in pars.iteritems():
            make_post_pulls(i, j, pulls_dir)
