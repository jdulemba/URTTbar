#! /bin/env python

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
from rootpy.plotting import views, Graph, Hist
from argparse import ArgumentParser
import ROOT

jobid = os.environ['jobid']
plotter = BasePlotter(
   'plots/%s/permFeatures' % jobid,
   defaults = {'save' : {'png' : True, 'pdf' : False}}
)

ttjets = urviews.NormalizedView(
	root_open('results/%s/permProbComputer/ttJets.root' % jobid)
	)
subs_and_col = [
	('semilep_visible_right','violet'),
	('other_tt_decay','blue'),
	('semilep_right_thad','cyan'),
	('semilep_right_tlep','red'),
	('semilep_wrong','#00960a'),
	]
samples = [
	views.TitleView(
		views.StyleView(
			views.SubdirectoryView(ttjets, i),
			linecolor = j,
      drawstyle = 'hist',
      legendstyle = 'l'
			), i)
	for i, j in subs_and_col
	]

plots = [
	'nusolver_chi2',
	'thfr',
	'tlfr',
	'wfr_jcosth_delta'
]
#plotter.set_subdir('shapes')
#for var in plots:
#	shapes = [i.Get('nosys/%s' % var) for i in samples]
#	legend = LegendDefinition(position='NE')
#	plotter.overlay(shapes, legend_def=legend, ytitle='a.u.', y_range='shape', ignore_style=True)
#	plotter.save(var)
	
systematics = [
	('nosys'   , 'violet'),
	('jes_down', 'blue'),
	('jes_up'  , 'cyan'),
	('met_down', 'red'),
	('met_up'  ,'#00960a'),
	]

sys_samples = [
	views.TitleView(
		views.StyleView(
			views.SubdirectoryView(samples[0], i),
			linecolor = j,
      drawstyle = 'hist',
      legendstyle = 'l'
			), i)
	for i, j in systematics
	]

#plotter.set_subdir('systematics')
#for var in plots:
#	shapes = [i.Get(var) for i in sys_samples]
#	legend = LegendDefinition(position='NE')
#	plotter.overlay(shapes, legend_def=legend, ytitle='a.u.', y_range='shape', ignore_style=True)
#	plotter.save(var)

from glob import glob
signals = glob('results/%s/permProbComputer/[AH]toTT*.root' % jobid)
Hs, As, HIs, AIs = [], [], [], []
for i in signals:
	base = os.path.basename(i)
	if base.startswith('HtoTT'):
		if '_Interf' in i:
			HIs.append(i)
		else:
			Hs.append(i)
	else:
		if '_Interf' in i:
			AIs.append(i)
		else:
			As.append(i)

higgses = views.SumView(*[
		root_open(i) for i in Hs
	])
pscals = views.SumView(*[
		root_open(i) for i in As
	])
hinter = views.SumView(*[
		root_open(i) for i in HIs
	])
ainter = views.SumView(*[
		root_open(i) for i in AIs
	])

pys_samples = [
	(higgses, 'H', 'blue'),
	(pscals , 'A', 'cyan'),
#	(hinter , 'HI', 'red'),
#	(ainter , 'AI', '#00960a'),
]

pys_samples = [
	urviews.NormalizedView(
		views.TitleView(
			views.StyleView(
				views.SubdirectoryView(i, 'semilep_visible_right'),
				linecolor = k,
				drawstyle = 'hist',
				legendstyle = 'l'
				), j
			) )
	for i, j, k in pys_samples
]
pys_samples.append(
	views.TitleView(
		samples[0],
		'ttjets'
		))

plotter.set_subdir('samples')
from pdb import set_trace
#set_trace()
for var in plots:
  shapes = [i.Get('nosys/%s' % var) for i in pys_samples]
  legend = LegendDefinition(position='NE')
  plotter.overlay(shapes, legend_def=legend, ytitle='a.u.', y_range='shape', ignore_style=True)
  plotter.save(var)
