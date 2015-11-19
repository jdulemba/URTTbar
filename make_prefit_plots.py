import rootpy.io as io
import rootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
rootpy.log["/data_views"].setLevel(rootpy.log.WARNING)
log = rootpy.log["/make_prefit_plots"]
log.setLevel(rootpy.log.INFO)
from rootpy.plotting import Hist
from argparse import ArgumentParser
import os
import re
import ROOT
from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.Utilities.prettyjson as prettyjson
from URAnalysis.Utilities.quad import quad
from styles import styles
from pdb import set_trace
from labels import set_pretty_label
import URAnalysis.Utilities.datacard as datacard

def fix_binning(postfit, correct_bin):
   ret = correct_bin.Clone()
   ret.Reset()
   ret.title = postfit.title
   ret.decorate(**postfit.decorators)
   assert(ret.nbins() == postfit.nbins())
   for rbin, pbin in zip(ret, postfit):
      rbin.value = pbin.value
      rbin.error = pbin.error
   return ret

parser = ArgumentParser()
parser.add_argument('varname')

args = parser.parse_args()
input_dir = 'plots/%s/ttxsec/%s' % (os.environ['jobid'], args.varname)

datacard = datacard.load('%s/%s.txt' % (input_dir, args.varname))
shapes = io.root_open(
   '%s/%s.root' % (input_dir, args.varname)
) #contains data

sig_yields = prettyjson.loads(
   open('%s/%s.json' % (input_dir, args.varname)).read()
   )

## binning_file = prettyjson.loads(
##    open('%s/%s.binning.json' % (input_dir, args.varname)).read()
## )

categories = datacard.bins
signals = set(datacard.signals)
regex = re.compile('^(?P<base_category>[A-Za-z0-9]+)_(?P<njets>\d+)Jets$')
plotter = BasePlotter(
   '%s/prefit' % input_dir,
   defaults = {'save' : {'png' : True, 'pdf' : False}},
   styles = {
      'only_thad_right' : {
         'legendstyle' : 'f',
         'drawstyle' : 'hist',
         'fillcolor' : '#d5aaaa',
         'linecolor' : '#d5aaaa',
         'title' : "tt, right t_{h}",
         'fillstyle': 'solid',
         },
      'qcd' : styles['QCD*'],
      'single_top' : styles['single*'],
      'tt_right' : styles['tt*'],
      'tt_wrong' : {
         'legendstyle' : 'f',
         'drawstyle' : 'hist',
         'fillcolor' : '#ab5555',
         'linecolor' : '#ab5555',
         'title' : "tt, wrong cmb",
         'fillstyle': 'solid',
         },
      'vjets' : styles['[WZ]Jets*'],
      'data_obs' : styles['data*'],
      'Total signal+background *' : {
         'legendstyle' : 'f',
         'drawstyle' : 'PE2',
         'linewidth' : 0,
         'title' : "postfit unc.",
         'markersize' : 0,
         'fillcolor' : 1,
         'fillstyle' : 3013
         }
      }
)
xlabel = set_pretty_label(args.varname)
groups = {}
for category in categories:
   m = regex.match(category)
   if not m:
      raise ValueError('Category name %s did not match the regex!' % category)
   base = m.group('base_category')
   njets = int(m.group('njets'))
   if base not in groups:
      groups[base] = []
   groups[base].append(category)

for base, categories in groups.items():
   sample_sums = {}
   plotter.set_subdir(base)
   for cat_name in categories:
      log.info("making plots for %s" % cat_name)
      shape_dir = shapes.Get(cat_name)
      data = shape_dir.data_obs
      data.title = 'data_obs'

      if data.title not in sample_sums:
      	sample_sums[data.title] = data.Clone()
      else:
      	sample_sums[data.title] += data

      samples = []
      for sample in datacard.processes:
         histo = shape_dir.Get(sample)
         if sample in signals:
            histo.Scale(sig_yields[cat_name][sample])
         if sample not in sample_sums:
            sample_sums[sample] = histo.Clone()
         else:
            sample_sums[sample] += histo
         samples.append(histo)

      stack = plotter.create_stack(*samples)
      legend = LegendDefinition(position='NE')
      plotter.overlay_and_compare( 
      	[stack], data,
      	writeTo=cat_name,
      	legend_def = legend,
        xtitle='discriminant',
        ytitle='Events',
        method='ratio'
      	)
   samples = [j for i, j in sample_sums.iteritems() if i <> 'data_obs' 
              if i <> 'postfit S+B']
   data = sample_sums['data_obs']

   plotter.overlay_and_compare(
   	[plotter.create_stack(*samples)], data,
   	writeTo=base,
   	legend_def = LegendDefinition(position='NE'),
    xtitle='discriminant',
    ytitle='Events',
    method='ratio'
   	)
	
      
