#! /bin/env python

from URAnalysis.PlotTools.Plotter import Plotter, BasePlotter
import URAnalysis.PlotTools.views as urviews
import os
import glob
from styles import styles
import sys
import logging
#import rootpy.plotting.views as views
import rootpy.plotting as plotting
views = plotting.views
from array import array
import ROOT
from pdb import set_trace
import rootpy
import rootpy.io as io
from fnmatch import fnmatch
import URAnalysis.Utilities.prettyjson as prettyjson
import URAnalysis.Utilities.quad as quad
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
from argparse import ArgumentParser
import math
from uncertainties import ufloat
from URAnalysis.Utilities.datacard import DataCard
from URAnalysis.Utilities.tables import latex_table
from URAnalysis.Utilities.latex import t2latex
import re
import itertools

parser = ArgumentParser()
parser.add_argument('mode', choices=['electrons', 'muons'], help='choose leptonic decay type')
#parser.add_argument('--pdfs', action='store_true', help='make plots for the PDF uncertainties')
args = parser.parse_args()

def syscheck(cmd):
	out = os.system(cmd)
	if out == 0:
		return 0
	else:
		raise RuntimeError("command %s failed executing" % cmd)

class HTTPlotter(Plotter):
	def __init__(self, mode, lumi=None):
		'Inits the plotter, mode is either "electrons" or "muons" and identifies what will be plotted'
		lumi = lumi if lumi > 0 else None
		filtering = lambda x: not os.path.basename(x).startswith('data') or \
			 (os.path.basename(x).startswith('data') and mode[:-1] in os.path.basename(x).lower())
		self.tt_to_use = 'ttJets'
		self.tt_shifted = {
			'mtop_up' : 'ttJets_mtopup',
			'mtop_down' : 'ttJets_mtopdown',
			'hadscale_up' : 'ttJets_scaleup', 
			'hadscale_down' : 'ttJets_scaledown',
			}
		jobid = os.environ['jobid']
		files = glob.glob('results/%s/htt_simple/*.root' % jobid)
		files = filter(filtering, files)
		logging.debug('files found %s' % files.__repr__())
		lumis = glob.glob('inputs/%s/*.lumi' % jobid)
		lumis = filter(filtering, lumis)
		logging.debug('lumi files found %s' % lumis.__repr__())
		
		outdir= 'plots/%s/htt/%s' % (jobid, mode)
		super(HTTPlotter, self).__init__(
			files, lumis, outdir, styles, None, lumi
			)

		#select only mode subdir
		for info in self.views.itervalues():
			info['view'] = views.SubdirectoryView(info['view'], mode)
			info['unweighted_view'] = views.SubdirectoryView(info['unweighted_view'], mode)

		self.defaults = {
			'blurb' : [13, self.views['data']['intlumi']]
			}
		self.jobid = jobid

		self.views['ttJets_preselection'] = self.views['ttJets']
		#{
		#	'view' : self.create_tt_subsample(
		#		['semilep_visible_right'],
		#		't#bar{t}',
		#		'#6666b3'
		#		)
		#	}
		self.views['ttJets_right'] = {
			'view' : self.create_tt_subsample(
				['semilep_visible_right'],
				't#bar{t}, right matching',
				'#6666b3'
				)
			}
		self.views['ttJets_partial'] = {
			'view' : self.create_tt_subsample(
				['semilep_right_tlep',  'semilep_right_thad', 'semilep_right_whad', 'semilep_wrong'],
				't#bar{t}, partial matching',
				'#ab5555'
				)
			}

		self.views['ttJets_other'] = {
			'view' : self.create_tt_subsample(
				['other'], 
				'Other tt decay',
				'#668db3',
				)
			}

		self.mc_samples = [
			'QCD*',
			'[WZ]Jets',
			#'WJets',
			#'ZJets',
			'single*',
			'ttJets',
			#'ttJets_other',
			#'ttJets_partial',
			#'ttJets_right',
			]

		self.card_names = {
			'qcd' : ['QCD*'],
			'vjets'		: ['[WZ]Jets'],
			'right_whad' : ['ttJets_sig'],
			'wrong_whad' : ['ttJets_bkg'],
			'nonsemi_tt' : ['ttJets_other'],
			'single_top' : ['single*'],
			'data_obs'	: ['data']
			}
		self.card_by_title = {
			'QCD' : 'qcd' ,
			'V + jets' : 'vjets' ,
			'single top' : 'single_top' ,
			'Other tt decay' : 'nonsemi_tt' ,
			't#bar{t}, wrong W_{h}' : 'wrong_whad', 
			't#bar{t}, right W_{h}' : 'right_whad', 
			'Observed' : 'data_obs'
			}
		self.signal = 'right_whad'
		self.signal_yields = {}

		self.systematics = {
			}
		self.card = None
		self.binning = {
			}

	def create_tt_subsample(self, subdirs, title, color='#9999CC'):
		dirmap = {
			'' : views.SumView(
				*[views.SubdirectoryView(self.views[self.tt_to_use]['view'], i) for i in subdirs]
				 )
			}
		for shift, view in self.tt_shifted.iteritems():
			dirmap[shift] = views.SumView(
				*[views.SubdirectoryView(self.views[view]['view'], '%s/nosys' % i) for i in subdirs]
				)
		
		return views.StyleView(
			views.TitleView(
				urviews.MultifileView(**dirmap),
				title
				),
			fillcolor = color,
			linecolor = color
			)

	def cut_flow(self):
		BasePlotter.set_canvas_style(self.canvas)
		BasePlotter.set_canvas_style(self.pad)
		lab_f1, _ = self.dual_pad_format()
		self.label_factor = lab_f1
		views_to_flow = filter(lambda x: 'ttJets' not in x and 'QCD' not in x, self.mc_samples)
		views_to_flow.append(self.tt_to_use)
		qcd_samples = [i for i in self.views if 'QCD' in i]
		samples = []

		for vtf in views_to_flow:
			histo = self.get_view(vtf).Get('cut_flow')
			print vtf, len(histo)
			self.keep.append(histo)
			samples.append(histo)

		#QCD may not have all the bins filled, needs special care
		qcd_histo = histo.Clone()			
		qcd_histo.Reset()
		for sample in qcd_samples:
			qcd_flow = self.get_view(sample).Get('cut_flow')
			qcd_histo = qcd_histo.decorate(
				**qcd_flow.decorators
				)
			qcd_histo.title = qcd_flow.title
			for sbin, qbin in zip(qcd_histo, qcd_flow):
				if sbin.overflow: continue
				sbin.value = qbin.value+sbin.value
				sbin.error = quad.quad(sbin.error, qbin.error)
		samples.append(qcd_histo)
		self.keep.append(qcd_histo)
		samples.sort(key=lambda x: x[-2].value)
		stack = plotting.HistStack()
		self.keep.append(stack)
		for i in samples:			
			stack.Add(i)

		cflow = ['']
		for idx in range(1,len(samples[0])):
			vals = {}
			for s in samples:
				vals['name'] = s.xaxis.GetBinLabel(idx)
				vals[s.title] = s[idx].value
			cflow.append(vals)

		self.style_histo(stack)
		self.style_histo(histo, **histo.decorators)

		histo.Draw() #set the proper axis labels
		histo.yaxis.title = 'Events'
		data = self.get_view('data').Get('cut_flow')
		for idx in range(1,len(samples[0])):
			cflow[idx][data.title] = data[idx].value
		smin = min(stack.min(), data.min(), 1.2)
		smax = max(stack.max(), data.max())
		histo.yaxis.range_user = smin*0.8, smax*1.2
		stack.Draw('same')
		data.Draw('same')
		self.keep.append(data)
		self.add_legend([stack, data], False, entries=len(views_to_flow)+1)
		self.pad.SetLogy()
		self.add_ratio_plot(data, stack, ratio_range=0.4)
		self.lower_pad.SetLogy(False)
		return cflow
		#cut_flow.GetYaxis().SetRangeUser(1, 10**7)

	def make_preselection_plot(self, *args, **kwargs):
		systematics = None
		if 'sys_effs' in kwargs:
			systematics = kwargs['sys_effs']
			del kwargs['sys_effs']
		mc_default = self.mc_samples
		self.mc_samples = ['QCD*', '[WZ]Jets', 'single*', 'ttJets_preselection']
		self.plot_mc_vs_data(*args, **kwargs)
		if systematics is not None:
			data = [i for i in self.keep if isinstance(i, ROOT.TH1)][0]
			mc_stack = [i for i in self.keep if isinstance(i, ROOT.THStack)][0]
			stack_sum = sum(mc_stack.hists)
			self.reset()
			dirname = args[0].split('/')[0]
			path = args[0]
			args = list(args)
			args[0] = path.replace(dirname, '%s_up' % systematics)
			kwargs['nodata'] = True
			self.plot_mc_vs_data(*args, **kwargs)
			stack_up = [i for i in self.keep if isinstance(i, ROOT.THStack)][0]
			self.reset()
			s_up = sum(stack_up.hists)
			for ibin, jbin in zip(stack_sum, s_up):
				ibin.error = quad.quad(ibin.error, abs(ibin.value - jbin.value))
			stack_sum.fillcolor = 'black'
			stack_sum.fillstyle = 3013
			stack_sum.title = 'uncertainty'
			stack_sum.drawstyle = 'pe2'
			stack_sum.markerstyle = 0
			plotter.overlay_and_compare(
				[mc_stack, stack_sum], data,
				xtitle = kwargs.get('xaxis',''),
				ytitle='Events', ignore_style=True,				
				method='datamc'
				)
			# Add legend
			self.pad.cd()
			self.add_legend(
				[mc_stack, stack_sum, data], kwargs.get('leftside', True), 
				entries=len(mc_stack.hists)+2
				)

		self.mc_samples = mc_default

	def make_sample_table(self):
		stack = [i for i in self.keep if isinstance(i, plotting.HistStack)][0]
		names = [i.title for i in stack.hists]
		yields = [i.Integral() for i in stack.hists]
		tot = sum(yields)
		yields = [100*i/tot for i in yields]
		mlen = max(len(i) for i in names)
		format = '%'+str(mlen)+'s     %.1f%%\n'
		fullname = os.path.join(self.outputdir, 'yields.raw_txt')
		with open(fullname, 'w') as f:
			for i in zip(names, yields):
				f.write(format % i)

plotter = HTTPlotter(args.mode)

variables = [
	("lep_b_pt" , "p_{T}(bl) (GeV)", 10, None, False),
	("had_b_pt" , "p_{T}(bh) (GeV)", 10, None, False),
	("wjets_pt" , "p_{T}(j) (GeV)", 10, None, False),
	("Whad_DR"  , "#DeltaR(jj)", 10, None, False),
	("Whad_pt"  , "p_{T}(W_{had}) (GeV)", 10, None, False),
  ("njets"	 , "# of selected jets", range(13), None, False),
  ("lep_pt"	, "p_{T}(l) (GeV)", 20, None, False),		
]

preselection = [
  (False, "njets"	 , "# of selected jets", range(13), None, False),
  (False, "jets_eta", "#eta(jet)", 10, None, False),
  (False, "jets_pt", "p_{T}(jet)", 10, None, False),
  (False, "lep_eta", "#eta(l)", 10, None, False),
  (False, "lep_pt", "p_{T}(l)", 10, None, False),
  (False, "nvtx", "# of reconstructed vertices", range(41), None, False),
  (False, "rho", "#rho", range(40), None, False),
	(True , "max_jets_CSV", 'Max CSV', 1, [0, 1], False),
	(False, "lep_iso", 'l rel Iso', 1, [0,1], False),
	(False, "lep_wp" , "electron wp", 1, None, False),
	(False, "MT" , "MT",  10, None, False),
]

permutations = [
	("comb_discriminant", "#lambda_{C}", 1, (8, 30), True),
	("nu_chi2", "#chi^{2}_{#nu}", 1, (2,6), True),
	("mass_discriminant", "#lambda_{M}", 1, None, True),
	("Wmasshad", "M(W_{h})", 2, (0,400), False),
	("tmasshad", "M(t_{h})", 2, None, False),
]

jet_categories = ["3jets", "4jets", "5Pjets"]

#cut flow
#flow = plotter.cut_flow()
#plotter.save('cut_flow')

plotter.set_subdir('preselection')
for logy, var, axis, rebin, x_range, leftside in preselection:
	plotter.make_preselection_plot(
		'nosys/preselection', var, sort=True,
		xaxis=axis, leftside=leftside, rebin=rebin, 
		show_ratio=True, ratio_range=0.5, xrange=x_range, logy=logy)		
	plotter.save(var)

for dirid in itertools.product(['looseNOTTight', 'tight'], ['inverted', 'pass']):
	tdir = '%s/%s' % dirid
	plotter.set_subdir(tdir)
	first = True
	for logy, var, axis, rebin, x_range, leftside in preselection:
		plotter.plot_mc_vs_data(
			'nosys/%s' % tdir, var, sort=True,
			xaxis=axis, leftside=leftside, rebin=rebin,
			show_ratio=True, ratio_range=0.5, xrange=x_range, logy=logy)
		if first:
			first = False
			plotter.make_sample_table()
		plotter.save(var)

#for var, axis, rebin, x_range, leftside in permutations:
#		plotter.make_preselection_plot(
#			'nosys/preselection', var, sort=True,
#			xaxis=axis, leftside=leftside, rebin=rebin, 
#			show_ratio=True, ratio_range=0.5)
#		plotter.save(var)
#
#for jdir in jet_categories:
#	plotter.set_subdir(jdir)
#	folder = 'nosys/%s' % jdir
#	for var, axis, rebin, x_range, leftside in variables+permutations:
#		plotter.plot_mc_vs_data(
#			folder, var, rebin, #sort=True,
#			xaxis=axis, leftside=False,
#			xrange=x_range, show_ratio=True, 
#			ratio_range=0.5)		
#		plotter.save(var)
		
