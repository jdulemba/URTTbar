#! /bin/env python
# coding: utf-8
from rootpy.io import root_open
from URAnalysis.Utilities.roottools import Envelope
from pdb import set_trace
from argparse import ArgumentParser
from glob import glob
from rootpy.plotting import Hist, Hist2D
import re
import os 
from rootpy import asrootpy
from URAnalysis.PlotTools.BasePlotter import BasePlotter
from uncertainties import ufloat
from rootpy.plotting import Canvas, Graph
from rootpy import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

parser = ArgumentParser()
parser.add_argument('wp')
parser.add_argument('--asimov', '-a', action='store_true', help='run on asimov scan')
parser.add_argument('--usemean', '-m', action='store_true', help='run on asimov scan')
parser.add_argument('--dir', '-d', default='toys_simple', help='run on asimov scan')
parser.add_argument('--precomputed', '-p', action='store_true', help='run on asimov scan')
args = parser.parse_args()

toydir= args.dir
basedir = '%s/plots/%s/ctageff/mass_discriminant/%s' % (os.environ['URA_PROJECT'], os.environ['jobid'], args.wp)
indir=os.path.join(basedir, toydir)
toy_regex=re.compile('toys_res\.\d+\.(?P<csf>[01]\.\d)\.(?P<lsf>[01]\.\d)\.root$')
def get_toy_point(fname):
   base=os.path.basename(fname)
   m = toy_regex.match(base)
   if not m:
      raise RuntimeError('problem matching %s' % base)
   return float(m.group('csf')), float(m.group('lsf'))

def covers(realvar, val):
   delta = realvar.value - val   
   err_idx = 1 if delta >= 0 else 0
   delta = abs(delta)
   if hasattr(realvar.error, '__getitem__'): #two sided error        
      return (delta < abs(realvar.error[err_idx]))
   else: #one sided error
      return (delta < abs(realvar.error))
   
if not args.precomputed:
   toys_files = glob('%s/toys_res*.root' % indir)
   csf_points = set()
   lsf_points = set()
   points = {}
   for toy in toys_files:
       key = get_toy_point(toy)
       if key in points: points[key].append(toy)
       else: points[key] = [toy]
       csf, lsf = key
       csf_points.add(csf)
       lsf_points.add(lsf)
       
   csf_points = sorted(list(csf_points))
   lsf_points = sorted(list(lsf_points))

   def make_binning(lst):
       lst = [(2*lst[0] - lst[1])] + lst
       lst.append(
           (2*lst[-1] - lst[-2])
           )
       ll = len(lst)
       ret = []
       for idx in range(ll-1):
           point = lst[idx]
           next  = lst[idx+1]
           ret.append((next+point)/2)
       return ret

   csf_points = make_binning(csf_points)
   lsf_points = make_binning(lsf_points)

   cbias   = Hist2D(csf_points, lsf_points)
   lbias   = Hist2D(csf_points, lsf_points)
   ccover  = Hist2D(csf_points, lsf_points)
   lcover  = Hist2D(csf_points, lsf_points)
   failing = Hist2D(csf_points, lsf_points)

   out = root_open('%s/toys_summary.root' % indir, 'w')
   hdir = out.mkdir('histograms')

   it=0
   mtot =0;
   for key, files in points.iteritems():
      #print key
      it +=1
      cval, lval = key
      csf = Hist(200, 0, 2) if not args.asimov else None
      lsf = Hist(200, 0, 2) if not args.asimov else None
      csf_prefit = Hist(200, 0, 2) if not args.asimov else None
      lsf_prefit = Hist(200, 0, 2) if not args.asimov else None
      csf2d = Hist2D(100, 0, 2, 100, 0, 2) if not args.asimov else None
      lsf2d = Hist2D(100, 0, 2, 100, 0, 2) if not args.asimov else None
      failed = 0
      tot = 0
      cover_tot = 0
      lcover_ok = 0
      ccover_ok = 0
      for toy_file in files:
         with root_open(toy_file) as mlfit:
            toys = [''] if args.asimov else [i.name for i in mlfit.keys() if i.name.startswith('toy_')]         
            tot += len(toys)
            for i in toys:            
               path = os.path.join(i, 'fit_s')
               toy = mlfit.Get(path)
               pars = asrootpy(toy.floatParsFinal())
               prefit = asrootpy(toy.floatParsInit())
               cover_tot += 1
               if args.asimov:
                  csf = float(pars['charmSF'].value)
                  lsf = float(pars['lightSF'].value)
               else:
                  csf.Fill(float(pars['charmSF'].value))
                  lsf.Fill(float(pars['lightSF'].value))
               if covers(pars['lightSF'], lval):
                  lcover_ok += 1
               if covers(pars['charmSF'], cval):
                  ccover_ok += 1
               if not args.asimov:
                  csf_prefit.Fill(float(prefit['charmSF'].value))
                  lsf_prefit.Fill(float(prefit['lightSF'].value))
                  csf2d.Fill(prefit['charmSF'].value, pars['charmSF'].value)
                  lsf2d.Fill(prefit['lightSF'].value, pars['lightSF'].value)
               toy.IsA().Destructor(toy)


      ccover.Fill(cval, lval, float(ccover_ok)/cover_tot)
      lcover.Fill(cval, lval, float(lcover_ok)/cover_tot)

      mtot = max(tot, mtot) if not args.asimov else 1
      fcn = BasePlotter.parse_formula(
         'amplitude*TMath::Gaus(x, mean, width)', 
         'amplitude[100, 0, 1000], mean[1, 0, 2], width[0.5, 0, 2]',
         [0,2]
         )
      failing.Fill(cval, lval, failed)
      xb = cbias.xaxis.find_bin(cval)
      yb = cbias.yaxis.find_bin(lval)
      
      if args.usemean:
         bias = ufloat(csf.GetMean(), 0)
      elif args.asimov:
         bias = ufloat(csf, 0)
      else:
         csf.Fit(fcn, 'LQN0')
         bias = ufloat(fcn['mean'].value, fcn['mean'].error)
      bias -= cval
      bias /= cval
      cbias[xb,yb].value = bias.nominal_value
      cbias[xb,yb].error = bias.std_dev
      
      if args.usemean:
         bias = ufloat(lsf.GetMean(), 0)
      elif args.asimov:
         bias = ufloat(lsf, 0)
      else:
         lsf.Fit(fcn, 'LQN0')
         bias = ufloat(fcn['mean'].value, fcn['mean'].error)
      bias -= lval
      bias /= lval
      lbias[xb,yb].value = bias.nominal_value
      lbias[xb,yb].error = bias.std_dev      
      if not args.asimov:
         hdir.WriteTObject(csf, 'csf_%.1f_%.1f' % (cval, lval))
         hdir.WriteTObject(lsf, 'lsf_%.1f_%.1f' % (cval, lval))
         hdir.WriteTObject(csf_prefit, 'csf_prefit_%.1f_%.1f' % (cval, lval))
         hdir.WriteTObject(lsf_prefit, 'lsf_prefit_%.1f_%.1f' % (cval, lval))
         hdir.WriteTObject(csf2d, 'csf2D_%.1f_%.1f' % (cval, lval))
         hdir.WriteTObject(lsf2d, 'lsf2D_%.1f_%.1f' % (cval, lval))

   for i in failing:
      new_val = i.value/mtot
      i.value = new_val
   out.WriteTObject(cbias.Clone(), 'cbias')
   out.WriteTObject(lbias.Clone(), 'lbias')
   out.WriteTObject(failing.Clone(), 'failing')
   out.WriteTObject(ccover.Clone(), 'ccover')
   out.WriteTObject(lcover.Clone(), 'lcover')
else:
   infile = root_open('%s/toys_summary.root' % indir)
   cbias = infile.cbias
   lbias = infile.lbias
   failing = infile.failing
   ccover = infile.ccover
   lcover = infile.lcover

#
# Draw plots
#
ROOT.gStyle.SetOptStat(0)

fitfile = root_open('%s/MaxLikeFit.root' % basedir)
fit = fitfile.fit_s
pars = asrootpy(fit.floatParsFinal())
fitted = Graph(1)
fitted.SetPoint(0, pars['charmSF'].value, pars['lightSF'].value)
fitted.markerstyle = 20
fitted.markersize  = 3
fitted.markercolor = '#009600'

def addlabels(h):
    h.xaxis.title = 'charm SF'
    h.yaxis.title = 'light SF'
    h.xaxis.SetTitleOffset(0.9)
    h.yaxis.SetTitleOffset(1.)

addlabels(lbias)
addlabels(cbias)
addlabels(failing)
addlabels(ccover)
addlabels(lcover)

mask = failing.Clone()
mask.Reset()

for lb, cb, mb, fb in zip(lbias, cbias, mask, failing):
   if fb.value > 0.2:
      lb.value = -100.
      cb.value = -100.
      mb.value = 1.
   else:
      mb.value = -1.

ROOT.gStyle.SetPaintTextFormat('.2f')
cbias.zaxis.range_user = (-1, 1)
lbias.zaxis.range_user = (-1, 1)
canvas = Canvas(800, 800)
canvas.margin = (0.15, 0.15, 0.15, 0.1)
from array import array
ROOT.TColor.CreateGradientColorTable(
    5, array('d', [0., 0.3, .5, 0.7, 1.]), 
    array('d', [0., 0., .95, 1., 1.]), 
    array('d', [0., 0., .95, 0., 0.]), 
    array('d', [1., 1., .95, 0., 0.]), 255
    )
ROOT.gStyle.SetNumberContours(255);

mask.fillstyle = 3354
mask.fillcolor = 1
mask.zaxis.range_user = (-0.9, 0.9)

lbias.Draw('colz text')
#mask.Draw('box same')
fitted.Draw('P SAME')
canvas.SaveAs('%s/lbias.png' % indir)

cbias.Draw('colz text')
#mask.Draw('box same')
fitted.Draw('P SAME')
canvas.SaveAs('%s/cbias.png' % indir)

ROOT.TColor.CreateGradientColorTable(
    2, array('d', [0., 1.]), 
    array('d', [0., 1.]), 
    array('d', [0., 0.]), 
    array('d', [1., 0.]), 256
    )
failing.Draw('colz')
failing.zaxis.range_user = (-0.1, 1)
fitted.Draw('P SAME')
canvas.SaveAs('%s/fail.png' % indir)

ROOT.TColor.CreateGradientColorTable(
    3, array('d', [0., 0.682, 1.]), 
    array('d', [0., .95, 1.]), 
    array('d', [0., .95, 0.]), 
    array('d', [1., .95, 0.]), 255
    )
ccover.Draw('colz text')
ccover.zaxis.range_user = (0., 1)
fitted.Draw('P SAME')
canvas.SaveAs('%s/ccover.png' % indir)

lcover.Draw('colz text')
lcover.zaxis.range_user = (0., 1)
fitted.Draw('P SAME')
canvas.SaveAs('%s/lcover.png' % indir)
