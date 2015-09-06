import rootpy.io as io
import rootpy
rootpy.log["/"].setLevel(rootpy.log.ERROR)
log = rootpy.log["/make_postfit_plots"]
log.setLevel(rootpy.log.INFO)
from argparse import ArgumentParser
import os
import re
from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
from styles import styles

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
         'name' : "tt, right t_{h}",
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
         'name' : "tt, wrong cmb",
         'fillstyle': 'solid',
         },
      'vjets *' : styles['[WZ]Jets*'],
      'data_obs' : styles['data*'],
      }
)

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

for base, categories in groups.iteritems():
   first = True
   sample_sums = {}
   plotter.set_subdir(base)
   for cat_name in categories:
      cat_dir = postfit_shapes.Get(cat_name)
      data_dir = datacard_file.Get(cat_name)
      sample_names = [i.name for i in cat_dir.keys() if not i.name.startswith('total')]
      samples = [cat_dir.Get(i) for i in sample_names]
      for i, j in zip(sample_names, samples):
         if first:
            sample_sums[i] = j.Clone()
         else:
            sample_sums[i] += j

      data = data_dir.data_obs
      data.title = 'data_obs'
      if first:
      	sample_sums[data.title] = data.Clone()
      else:
      	sample_sums[data.title] += data
      first = False
      	
      plotter.overlay( 
      	[plotter.create_stack(*samples), data],
      	writeTo=cat_name,
      	legend_def = LegendDefinition(position='NE')
      	)
   samples = [j for i, j in sample_sums.iteritems() if i <> data.title]
   data = sample_sums[data.title]
   plotter.overlay(
   	[plotter.create_stack(*samples), data],
   	writeTo=cat_name,
   	legend_def = LegendDefinition(position='NE'),
   	)
	
      

