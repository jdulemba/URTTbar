import rootpy.io as io
import rootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
rootpy.log["/data_views"].setLevel(rootpy.log.WARNING)
log = rootpy.log["/make_postfit_plots"]
log.setLevel(rootpy.log.WARNING)
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

mlfit_file = io.root_open(
   '%s/%s.mlfit.root' % (input_dir, args.varname)
   )
datacard_file = io.root_open(
   '%s/%s.root' % (input_dir, args.varname)
) #contains data
binning_file = prettyjson.loads(
   open('%s/%s.binning.json' % (input_dir, args.varname)).read()
)

postfit_shapes = mlfit_file.shapes_fit_s
categories = [i.name for i in postfit_shapes.keys()]
regex = re.compile('^(?P<base_category>[A-Za-z0-9]+)_(?P<njets>\d+)Jets$')
plotter = BasePlotter(
   '%s/postfit' % input_dir,
   defaults = {'save' : {'png' : True, 'pdf' : False}},
   styles = {
      'only_thad_right *' : {
         'legendstyle' : 'f',
         'drawstyle' : 'hist',
         'fillcolor' : '#d5aaaa',
         'linecolor' : '#d5aaaa',
         'title' : "tt, right t_{h}",
         'fillstyle': 'solid',
         },
      'qcd *' : styles['QCD*'],
      'single_top *' : styles['single*'],
      'tt_right *' : styles['tt*'],
      'tt_wrong *' : {
         'legendstyle' : 'f',
         'drawstyle' : 'hist',
         'fillcolor' : '#ab5555',
         'linecolor' : '#ab5555',
         'title' : "tt, wrong cmb",
         'fillstyle': 'solid',
         },
      'vjets *' : styles['[WZ]Jets*'],
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

for base, categories in groups.items()[:1]:
   first = True
   sample_sums = {}
   plotter.set_subdir(base)
   for cat_name in categories:
      data_dir = datacard_file.Get(cat_name)
      cat_dir = postfit_shapes.Get(cat_name)
      data = data_dir.data_obs
      data.title = 'data_obs'
      hsum = fix_binning(cat_dir.total, data)
      if first:
      	sample_sums[data.title] = data.Clone()
        sample_sums['postfit S+B'] = hsum.Clone()
      else:
      	sample_sums[data.title] += data
        sample_sums['postfit S+B'] += hsum

      sample_names = [i.name for i in cat_dir.keys() if not i.name.startswith('total')]
      samples = [fix_binning(cat_dir.Get(i), data) for i in sample_names]
      for i, j in zip(sample_names, samples):
         if first:
            sample_sums[i] = j.Clone()
         else:
            sample_sums[i] += j

      first = False
      
      stack = plotter.create_stack(*samples)
      legend = LegendDefinition(position='NE')
      plotter.overlay( 
      	[stack, hsum, data],
      	writeTo=cat_name,
      	legend_def = legend,
        xtitle='discriminant',
        ytitle='Events'
      	)
   samples = [j for i, j in sample_sums.iteritems() if i <> 'data_obs' 
              if i <> 'postfit S+B']
   data = sample_sums['data_obs']
   hsum = sample_sums['postfit S+B']

   plotter.overlay(
   	[plotter.create_stack(*samples), hsum, data],
   	writeTo=base,
   	legend_def = LegendDefinition(position='NE'),
    xtitle='discriminant',
    ytitle='Events'
   	)
	
      
#global summary plots are best gathered from the samples normalizations
#to take into account at least a bit of nuisance correlations
norms = mlfit_file.norm_fit_s
edges = set()
binning = binning_file[args.varname]
for i in binning.itervalues():
   edges.add(i['up_edge'])
   edges.add(i['low_edge'])
edges = sorted(list(edges))
template_hist = Hist(edges)
samples = set(i.name.split('/')[-1] for i in norms)
fit_histos = {}
for name in samples:
   fit_histos[name] = template_hist.Clone()
   fit_histos[name].title = '%s ' % name
data = template_hist.Clone()
data.title = 'data_obs'

for categories in groups.itervalues():
   #check that all categories have identical bin index
   assert(len(set(binning[i]['idx'] for i in categories)) <= 1)
   bin_idx = binning[categories[0]]['idx']+1
   for sample in samples:
      val = sum(
         norms['%s/%s' % (i, sample)].value 
         for i in categories
         if '%s/%s' % (i, sample) in norms
         )
      err = quad(*[
            norms['%s/%s' % (i, sample)].error
            for i in categories
            if '%s/%s' % (i, sample) in norms
            ])
      log.debug("Setting value of %s bin %i to %.2f +/- %.2f" % (sample, bin_idx, val, err))
      fit_histos[sample][bin_idx].value = val
      fit_histos[sample][bin_idx].error = err
   #WARNING! this will break if the jet categories have 
   #different disciminator binning!
   data_of_cat = sum(
      datacard_file.Get('%s/data_obs' % i) 
      for i in categories      
      )
   integral_err = ROOT.Double()
   data[bin_idx].value = data_of_cat.IntegralAndError(1, data_of_cat.nbins(), integral_err)
   data[bin_idx].error = integral_err

samples_sum = sum(i for i in fit_histos.itervalues())
samples_sum.title = 'Total signal+background '
plotter.set_subdir('')
plotter.overlay_and_compare(
   [plotter.create_stack(*fit_histos.values()), samples_sum],
   data,
   writeTo='%s_postfit_compare' % args.varname,
   legend_def = LegendDefinition(position='NE'),
   xtitle=xlabel,
   ytitle='Events'
   )
plotter.overlay(
   [plotter.create_stack(*fit_histos.values()), samples_sum, data],
   writeTo='%s_postfit' % args.varname,
   legend_def = LegendDefinition(position='NE'),
   xtitle=xlabel,
   ytitle='Events'
   )
