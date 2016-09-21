#! /bin/env python
'''
Top Spin Gen Plotting macro
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
from rootpy.plotting import views, Graph, Hist
from argparse import ArgumentParser
from URAnalysis.Utilities.roottools import slice_hist
import ROOT

project = os.environ['URA_PROJECT']
##parser = ArgumentParser()
##parser.add_argument('todo', nargs='+', help='things to do', choices=['shapes', 'comparisons', 'effs'])
##args = parser.parse_args()

jobid = os.environ['jobid']
plotter = BasePlotter(
   '%s/plots/%s/tsgen' % (project, jobid),
   defaults = {'clone' : False, 'save' : {'png' : True, 'pdf' : False}}
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
   '' : {
      'title' : 'right perm',
      'linecolor' : 'blue',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'wrong' : {
      'title' : 'wrong perm',
      'linecolor' : 'red',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'wrongdtype' : {
      'title' : 'wrong (d-type)',
      'linecolor' : 'black',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'wrongutype' : {
      'title' : 'wrong (u-type)',
      'linecolor' : '#00960a',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      }
}

tfile = root_open('%s/results/%s/permProbComputer/ttJets.root' % (project, jobid))
def make_view(tf, subdir, title, colour):
	return urviews.NormalizedView(
		views.TitleView(
			views.StyleView(
				views.SubdirectoryView(
					tf,
					subdir
					),
				linecolor = colour,
				drawstyle = 'hist',
				linewidth = 3,
				legendstyle = 'l'
				),
			title
			)
		)

samples = [
	make_view(tfile, 'semilep_visible_right/nosys', 'right cmb', '#6666b3'),
	make_view(tfile, 'semilep_right_thad/nosys'   , 'wrong tl' , '#FFCC66'),
	make_view(tfile, 'semilep_right_tlep/nosys'   , 'wrong th' , '#2aa198'),
	make_view(tfile, 'semilep_wrong/nosys'        , 'wrong'    , '#0055ff')
]

variables = [
	 (1, 'NE', 'lb_ratio', 'p_{T}(l)/p_{T}(b_{l})'),
	 (1, 'NE', 'w1b_ratio', 'p_{T}(j1)/p_{T}(b_{h})'),
	 (1, 'NE', 'w2b_ratio', 'p_{T}(j2)/p_{T}(b_{h})'),
]

vars2D = [
	('mWhad_vs_mtophad', 'm(W_{h})', 'm(t_{h})'),
	('wjets_qgt', 'qgtag', 'qgtag'),
	('bjets_qgt', 'qgtag', 'qgtag'),
	('lb_w2b_ratio', 'p_{T}(l)/p_{T}(b_{l})', 'p_{T}(j2)/p_{T}(b_{h})'),
## 	('wjets_cMVA', 'cMVA^{5}', 'cMVA^{5}'),
## 	('bjets_cMVA', 'cMVA^{5}', 'cMVA^{5}'),
## 	('wjets_cMVA_WP', 'cMVA WP', 'cMVA WP'),
## 	('bjets_cMVA_WP', 'cMVA WP', 'cMVA WP'),
	]

projections = [
	(1, 1, 'NE', 'wjets_cMVA', 'cMVA^{5}', 'cMVA^{5}'),
	(1, 1, 'NE', 'bjets_cMVA', 'cMVA^{5}', 'cMVA^{5}'),
	(1, 1, 'NE', 'wjets_cMVA_WP', 'cMVA WP', 'cMVA WP'),
	(1, 1, 'NE', 'bjets_cMVA_WP', 'cMVA WP', 'cMVA WP'),
	(1, 1, 'NE', 'wjets_qgt', 'qgtag', 'qgtag'),
	(1, 1, 'NE', 'bjets_qgt', 'qgtag', 'qgtag'),
]

		
for _, lpos, name, xaxis in variables:
	histos = [i.Get(name) for i in samples]	
	legend = LegendDefinition(position=lpos)
	plotter.overlay(histos, legend_def=legend, xtitle=xaxis, ytitle='a.u.', y_range='shape', ignore_style=True)
	#plotter.add_legend(histos, False)
	plotter.save(name)

for _, _, lpos, name, xaxis1, xaxis2 in projections:
	histos = [i.Get(name) for i in samples]	
	histos_x = [asrootpy(i.projection_x()) for i in histos]
	histos_y = [asrootpy(i.projection_y()) for i in histos]
	for i, j, k in zip(histos_x, histos_y, histos):
		i.decorate(**k.decorators)
		j.decorate(**k.decorators)
		
	legend = LegendDefinition(position=lpos)
	plotter.overlay(histos_x, legend_def=legend, xtitle=xaxis1, ytitle='a.u.', y_range='shape', ignore_style=True)
	#plotter.add_legend(histos_x, False)
	plotter.save('%s_X' % name)
	plotter.overlay(histos_y, legend_def=legend, xtitle=xaxis2, ytitle='a.u.', y_range='shape', ignore_style=True)
	#plotter.add_legend(histos_y, False)
	plotter.save('%s_Y' % name)

for name, xaxis, yaxis in vars2D:
	histos = [i.Get(name) for i in samples]
	for h in histos:
		h.drawstyle = 'colz'
		plotter.plot(h, xtitle=xaxis, ytitle=yaxis, ignore_style=True)
		plotter.save(('%s_%s' % (name, h.title)).replace(' ', '_'))
