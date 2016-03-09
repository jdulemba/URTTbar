'''
Top Spin Gen Plotting macro
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import os
from rootpy.io import root_open
from rootpy import asrootpy

jobid = os.environ['jobid']
plotter = BasePlotter(
   'plots/%s/tsgen' % jobid,
   defaults = {'save' : {'png' : True, 'pdf' : False}}
)

styles = {
   'lep' : {
      'title' : 'l',
      'linecolor' : 'blue',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'nu' : {
      'title' : '#nu',
      'linecolor' : 'cyan',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'utype' : {
      'title' : 'up-type',
      'linecolor' : 'red',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'dtype' : {
      'title' : 'down-type',
      'linecolor' : 'violet',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'wrongd' : {
      'title' : 'wrong (d-type)',
      'linecolor' : 'black',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'wrongu' : {
      'title' : 'wrong (u-type)',
      'linecolor' : '#00960a',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      }
}

tfile = root_open('results/%s/topspin_gen/ttJets.root' % jobid)

variables = [
   (1, 'NE', 'helframe_costheta_%s', 'cos #theta*', ['dtype', 'utype', 'lep', 'nu']),
   (2, 'SE', 'helframe_cosdelta_lep_%s'   , 'cos #delta(l, q)', ['dtype', 'utype']),
   (1, 'NE', 'helframe_prodcosth_lep_%s'  , 'cos #theta*_{l} cos #theta*_{q}', ['dtype', 'utype']),
   (1, 'NE', 'labframe_cosdeltaphi_lep_%s', 'cos #Delta#varphi(l, q)', ['dtype', 'utype']),
]

for tdir in ['matched', 'part_acceptance', 'parton']:
   path_dir = tdir
   plotter.set_subdir(path_dir)
   for rbin, lpos, var, xaxis, types in variables:
      to_draw = []
      add = ['wrongu', 'wrongd'] if tdir == 'matched' else []         
      base_path = '%s/%s' % (path_dir, var)            
      for ptype in types+add:
         if ptype.startswith('wrong'):
            path = base_path % ('dtype' if ptype.endswith('d') else 'utype')
            path = path.replace(tdir, 'wrong')
         else:
            path = base_path % ptype
         histo = asrootpy(tfile.Get(path).ProjectionX()).Clone()
         if histo.Integral():
            histo.Scale(1./histo.Integral())
         else:
            print '%s has not integral!' % path
         histo.Rebin(rbin)
         plotter.set_histo_style(histo, **styles[ptype])
         to_draw.append(histo)
      legend = LegendDefinition(position=lpos)
      plotter.overlay(to_draw, legend_def=legend, xtitle=xaxis, ytitle='a.u.')
      plotter.save(var.replace('_%s', ''))
