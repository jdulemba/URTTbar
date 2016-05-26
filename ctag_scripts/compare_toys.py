#! /bin/env python

from rootpy.io import root_open
from glob import glob
from argparse import ArgumentParser
from rootpy.plotting import Hist, Hist2D
import os
from URAnalysis.PlotTools.BasePlotter import BasePlotter
import URAnalysis.Utilities.prettyjson as prettyjson
from pdb import set_trace

parser = ArgumentParser()
parser.add_argument('wp')
parser.add_argument('--filter', default='*')
parser.add_argument('--allbins', action='store_true')
args = parser.parse_args()
jobid=os.environ['jobid']
urap=os.environ['URA_PROJECT']

plotter = BasePlotter('%s/plots/%s/ctageff/mass_discriminant/%s/toys_cmp' % (urap, jobid, args.wp), defaults={'save' :{'pdf' : False}})
asimovs = glob('%s/plots/%s/ctageff/mass_discriminant/%s/asimovs/toys_src.%s.root' % (urap, jobid, args.wp, args.filter))
sorting = lambda x: x[0]

def flush(lst, **decorators):
   ret = Hist(len(lst), 0, len(lst))
   for idx, info in enumerate(lst):
      ret[idx+1].value = info[1]
   ret.decorate(**decorators)
   return ret

jmap = {}

for fname in asimovs:
   total = Hist(50, 0, 100000)
   savename = os.path.basename(fname).replace('toys_src.', '').replace('.root', '')
   simples = fname.replace('toys_src', 'toys_src.*').replace('asimovs', 'toys_simple')   
   simples = glob(simples)
   toys = {}
   first = True
   for sfile in simples:
      #print "  ", sfile
      with root_open(sfile) as simple:
         for i in range(1, 201):
            toy = simple.Get('toys/toy_%d' % i)
            toy_nfo = {}
            tot=0
            for bin in toy:
               if args.allbins:
                  key = bin.leaves['CMS_channel'].index, bin.leaves['CMS_th1x'].value
               else:
                  key = bin.leaves['CMS_channel'].index
               if key not in toy_nfo:
                  #toys[key] = []
                  toy_nfo[key] = 0
               #toys[key].append(bin.weight)
               toy_nfo[key] += bin.weight
               tot += bin.weight
            total.Fill(tot)
            for key, v in toy_nfo.iteritems():
               if key not in toys:
                  toys[key] = [] 
               toys[key].append(v)
            

   for val in toys.itervalues():
      val.sort()
      
   avgs = [(k, sum(v)/len(v)) for k, v in toys.iteritems()]
   meds = [(k, v[len(v)/2]  ) for k, v in toys.iteritems()]
   
   avgs.sort(key=sorting)
   meds.sort(key=sorting)
   
   height = total.max()
   asimov = {}
   tot = 0
   with root_open(fname) as asifile:
      toy = asifile.toys.toy_asimov
      for bin in toy:
         if args.allbins:
            key = bin.leaves['CMS_channel'].index, bin.leaves['CMS_th1x'].value
         else:
            key = bin.leaves['CMS_channel'].index
            
         tot += bin.weight
         if key in asimov:
            asimov[key] += bin.weight
         else:
            asimov[key] = bin.weight

   asitot = total.Clone()
   asitot.Reset()
   asitot.Fill(tot, height)
   asitot.markerstyle = 20
   asitot.markercolor = 'red'

   f1 = plotter.parse_formula('height*TMath::Poisson(x, mean)', 'height[1], mean[%f]' % (tot), [0, 100000])
   f1.linecolor = 'red'
   f1.linewidth = 2
   peak = f1(tot)
   f1 = plotter.parse_formula('height*TMath::Poisson(x, mean)', 'height[%f], mean[%f]' % (height/peak, tot), [0, 100000])
   

   canlog = all(i > 0 for i in asimov.itervalues())
   #jmap[savename] = asimov
   asimov = asimov.items()
   asimov.sort(key=sorting)
   
   havg = flush(avgs, linewidth=2, linecolor='red')
   hmed = flush(meds, linewidth=2, linecolor='blue')
   hasi = flush(asimov, linewidth=2, linecolor='black')
   plotter.overlay_and_compare([havg, hmed], hasi, method='ratio')
   ##if canlog:
   ##   plotter.pad.SetLogy(True)
   plotter.save(savename)
   
   plotter.overlay([total, asitot, f1])
   plotter.save('tot'+savename)

#with open('%s/asimovs.json' % plotter.outputdir, 'w') as json:
#   json.write(prettyjson.dumps(jmap))
