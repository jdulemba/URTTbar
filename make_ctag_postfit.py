import rootpy.io as io
import rootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
rootpy.log["/data_views"].setLevel(rootpy.log.WARNING)
log = rootpy.log["/make_postfit_plots"]
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

def fix_binning(postfit, correct_bin):
   ret = correct_bin.Clone()
   ret.Reset()
   ret.title = postfit.title
   ret.decorate(**postfit.decorators)
   assert(ret.nbins() <= postfit.nbins())
   for rbin, pbin in zip(ret, postfit):
      rbin.value = pbin.value
      rbin.error = pbin.error
   return ret

parser = ArgumentParser()
parser.add_argument('working_point')

args = parser.parse_args()
input_dir = 'plots/%s/ctageff/mass_discriminant/%s' % (os.environ['jobid'], args.working_point)

mlfit_file = io.root_open(
   '%s/MaxLikeFit.root' % input_dir
   )
datacard_file = io.root_open(
   '%s/datacard.root' % input_dir
) #contains data

def ordering(histo):
   multiplier = 1.
   if histo.title.startswith('right_whad'):
      multiplier = 100000.
   elif histo.title.startswith('wrong_whad'):
      multiplier = 100.
   return multiplier*histo.Integral()

plotter = BasePlotter(
   '%s/postfit' % input_dir,
   defaults = {'save' : {'png' : True, 'pdf' : False}},
   styles = {
      'right_whad *' : {
         'legendstyle' : 'f',
         'drawstyle' : 'hist',
         'fillcolor' : '#6666b3',
         'linecolor' : '#6666b3',
         'title' : "tt, right W_{h}",
         'fillstyle': 'solid',
         },
      'qcd *' : styles['QCD*'],
      'single_top *' : styles['single*'],
      'wrong_whad *' : {
         'legendstyle' : 'f',
         'drawstyle' : 'hist',
         'fillcolor' : '#ab5555',
         'linecolor' : '#ab5555',
         'title' : "tt, wrong W_{h}",
         'fillstyle': 'solid',
         },
      'nonsemi_tt *' : {
         'legendstyle' : 'f',
         'drawstyle' : 'hist',
         'fillcolor' : '#668db3',
         'linecolor' : '#668db3',
         'title' : "tt, non semi-lep.",
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

for pfit, tdir in [('postfit', 'shapes_fit_s'), ('prefit', 'shapes_prefit')]:
   shapes = mlfit_file.Get(tdir)
   plotter.set_outdir('%s/%s' % (input_dir,pfit))
   categories = [i.name for i in shapes.keys()]

   for cat_name in categories:
      log.info("making plots for %s" % cat_name)
      data_dir = datacard_file.Get(cat_name)
      cat_dir = shapes.Get(cat_name)
      data = data_dir.data_obs
      data.title = 'data_obs'
      hsum = fix_binning(cat_dir.total, data)

      sample_names = [i.name for i in cat_dir.keys() if not i.name.startswith('total')]
      samples = [fix_binning(cat_dir.Get(i), data) for i in sample_names]
      samples.sort(key=ordering)

      stack = plotter.create_stack(*samples, sort=False)
      legend = LegendDefinition(position='NE')
      plotter.overlay( 
      	[stack, hsum, data],
      	writeTo=cat_name,
      	legend_def = legend,
        xtitle='discriminant',
        ytitle='Events'
        #,method='ratio'
      	)
