#! /bin/env python

from rootpy.io import root_open as ropen
from rootpy import asrootpy as rpy
import ROOT
from argparse import ArgumentParser
import os
from glob import glob
from math import sqrt
from pdb import set_trace
#load CSV interfacing tools
ROOT.gROOT.ProcessLine('.L ../../URAnalysis/AnalysisFW/python/BTagCalibrationStandalone.cc+')
ROOT.TH1.AddDirectory(False)

def make_csv_entry(value, wpoint, flav, shift='central'):
   btagpars = ROOT.BTagEntry.Parameters(
      wpoint,  # 0 means loose
      "ttbar",  # type of measurement / ask B-Tag conveners
      shift,  # systematic
      flav,  # 0 means b-flavor
      -2.4, 2.4,   # eta range
      25, 200,   # pt range
      0, 1  # discr. range (actually not needed for loose op)
      )
   entry = ROOT.BTagEntry('%f' % value, btagpars)
   return entry

parser = ArgumentParser()
parser.add_argument('algo', help='algo to dump (CSV, CTAG etc..)')
## parser.add_argument('--wps', default='*',
##                     help='choose the working points to use')
args = parser.parse_args()
jobid = os.environ['jobid']

calibration = ROOT.BTagCalibration(args.algo)
def add_entry(value, wpoint, flav, shift='central'):
   calibration.addEntry(
      make_csv_entry(
         value,
         wpoint,
         flav,
         shift
         )
      )

class CorrelationMatrix(object):
   def __init__(self, matrix):
      self.matrix = matrix
      self.xmapping = {matrix.xaxis.GetBinLabel(idx) : idx for idx in range(1, matrix.nbins()+1)}
      self.ymapping = {matrix.yaxis.GetBinLabel(idx) : idx for idx in range(1, matrix.nbins()+1)}
      self.nbins = matrix.nbins()

   def __contains__(self, val):
      return val in self.xmapping

   def __getitem__(self, names):
      v1, v2 = names
      idx1 = self.xmapping[v1]
      idx2 = self.ymapping[v2]
      return self.matrix[idx1,idx2].value

wps = glob('plots/%s/ctageff/mass_discriminant/%s*' % (jobid, args.algo.lower()))
wp_mapping = { i:j for j, i in enumerate(['Loose', 'Medium', 'Tight'])}
sfs = ['bottomSF', 'charmSF', 'lightSF']

for wp in wps:
   wp_name = os.path.basename(wp)[len(args.algo):]
   wp_idx = wp_mapping[wp_name]
   shifts = ['%s/MaxLikeFit.root' % wp, '%s/MaxLikeFitStatistic.root' % wp] + glob('%s/sys_breakdown/*.root' % wp)
   first = True
   stat_errs = {}
   correlations = CorrelationMatrix(
      rpy(ropen('%s/MaxLikeFit.root' % wp).fit_s.correlationHist())
      )
   for shift in shifts:
      tfile = ropen(shift)
      pars = rpy(tfile.fit_s.floatParsFinal())
      for sf_idx, sf_name in enumerate(sfs):
         if sf_name in pars:
            var = pars[sf_name]
            err = max(abs(i) for i in var.error) if hasattr(var.error, '__len__') else var.error
            shift_name = ''
            sign = +1
            if 'sys_breakdown' in shift:
               base = os.path.basename(shift)
               shift_name = '_%s' % base.replace('.root', '')
               sign = -1 if shift_name in correlations and correlations[sf_name, shift_name] < 0 else +1
            elif shift.endswith('Statistic.root'):
               shift_name = '_statistic'
               
            if shift.endswith('MaxLikeFit.root'):
               #dump central
               add_entry(var.value, wp_idx, sf_idx)
            elif shift.endswith('Statistic.root'): 
               stat_errs[sf_name] = err
            else:
               if stat_errs[sf_name] > err:
                  print wp_name, shift_name, 'with negative error! (%f)' % (stat_errs[sf_name] - err)
               err = sqrt(abs(err**2 - stat_errs[sf_name]**2))
            #dump shifted up/down            
            add_entry(var.value+sign*err, wp_idx, 
               sf_idx, 'up%s' % shift_name)
            add_entry(var.value-sign*err, wp_idx, 
               sf_idx, 'down%s' % shift_name)               
      first = False

with open('%s.csv' % args.algo, 'w') as f:
   f.write(calibration.makeCSV())

