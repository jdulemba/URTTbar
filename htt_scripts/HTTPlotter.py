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
parser.add_argument('--preselection', action='store_true', help='')
parser.add_argument('--plots', action='store_true', help='')
parser.add_argument('--flow', action='store_true', help='')
parser.add_argument('--shapes', action='store_true', help='')
parser.add_argument('--all', action='store_true', help='')
parser.add_argument('--btag', action='store_true', help='')
parser.add_argument('--card', action='store_true', help='')
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
			#defaults = {'save' : {'png' : True, 'pdf' : False}}
			)

		#select only mode subdir
		for info in self.views.itervalues():
			info['view'] = views.SubdirectoryView(info['view'], mode)
			info['unweighted_view'] = views.SubdirectoryView(info['unweighted_view'], mode)
		self.views['data']['view'] = urviews.BlindView(
			self.views['data']['view'], 
			'\w+/tight/MTHigh/(:?(:?m_tt)|(:?.+_ctstar)|(:?cdelta_ld)|(:?hframe_ctheta_d))'
			)

		self.defaults = {
			'blurb' : [13, self.views['data']['intlumi']]
			}
		self.jobid = jobid

		self.views['ttJets_preselection'] = self.views['ttJets']

		self.views['ttJets_right'] = {
			'view' : self.create_tt_subsample(
				['right'],
				't#bar{t} matched',
				'#6666b3'
				)
			}
		self.views['ttJets_matchable'] = {
			'view' : self.create_tt_subsample(
				['matchable'], #,  'right_th', 'wrong'],
				't#bar{t} matchable',
				'#dddddd'
				)
			}
		self.views['ttJets_unmatchable'] = {
			'view' : self.create_tt_subsample(
				['unmatchable'],#  'right_th', 'wrong'],
				't#bar{t} unmatchable',
				'#ab5555'
				)
			}

		self.views['ttJets_other'] = {
			'view' : self.create_tt_subsample(
				['noslep'], 
				'Other t#bar{t}',
				'#668db3',
				)
			}

		##
		## signal sub-samples
		## 
		added_samples = []
		for sample in self.views.keys():
			if sample.startswith('AtoTT_'): 
				raise ValueError("I did not implement it yet remember to swap the H and A")
			if sample.startswith('HtoTT_'):# or sample.startswith('AtoTT_'):
				_, mval, width, pI = tuple(sample.split('_'))				
				new_name = ['ggA']
				interference = False
				if pI == 'Int':
					interference = True
					new_name.append('int')
				new_name.extend(['neg', mval[1:], width])
				new_name = '_'.join(new_name)
				self.views[new_name] = {
					'view' : self.create_subsample(sample, ['negative'], '%s negative' % sample, color='#9999CC')
					}
				added_samples.append(new_name)
				new_name = new_name.replace('neg_' if interference else '_neg_', '')
				self.views[new_name.replace('neg_' if interference else '_neg_', '')] = {
					'view' : self.create_subsample(sample, ['positive'], '%s positive' % sample, color='#9999CC')
					}
				added_samples.append(new_name)

		self.mc_samples = [
			'QCD*',
			'[WZ][WZ]',
			'[WZ]Jets',
			'tt[WZ]*',
			## #'WJets',
			#'ZJets',
			'single*',
			#'ttJets',
			'ttJets_other',
			'ttJets_unmatchable',
			'ttJets_matchable',
			'ttJets_right',
			]

		self.card_names = {
			#synced
			'QCDmujets' if mode == 'muons' else 'QCDejets' : ['QCD*'],
			'TT' : ['ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right'],
			'VV' : ['[WZ][WZ]'],
			#contested
			'TTV': ['tt[WZ]*'],
			'single_top' : ['single*'],
			'vjets'		: ['[WZ]Jets'],
			'data_obs'	: ['data']
			}
		self.card_names.update({i : [i] for i in added_samples})

		self.systematics = {
			#synched
			'lumi_13TeV' : {
				'type' : 'lnN',
				'samples' : ['.*'],
				'categories' : ['.*'],
				'value' : 1.027,
				},
			'CMS_scale_j_13TeV' : {
				'samples' : ['.*_whad', 'nonsemi_tt', 'single_top'],
				'categories' : ['.*'],
				'type' : 'shape',
				'+' : lambda x: x.replace('nosys', 'jes_up'),
				'-' : lambda x: x.replace('nosys', 'jes_down'),
				'constants' : ('jes_down', 'jes_up'),
				'value' : 1.00,
				},
			}
		{
			#to sync
			'lepEff' : {
				'type' : 'lnN',
				'samples' : ['.*'],
				'categories' : ['.*'],
				'value' : 1.02,
				},
			'ttxsec' : {
				'type' : 'lnN',
				'samples' : ['TT'],
				'categories' : ['.*'],
				'value' : 1.16,
				},			
			'singletxsec' : {
				'type' : 'lnN',
				'samples' : ['single_top'],
				'categories' : ['.*'],
				'value' : 1.50,
				},			
			'qcdxsec' : {
				'type' : 'lnN',
				'samples' : ['qcd'],
				'categories' : ['.*'],
				'value' : 1.70,
				},			
			'vjetxsec' : {
				'type' : 'lnN',
				'samples' : ['vjets'],
				'categories' : ['.*'],
				'value' : 1.10,
				},			
			'pu' : {
				'samples' : ['.*'],
				'categories' : ['.*'],
				'type' : 'lnN',
				'+' : lambda x: x.replace('nosys', 'pu_up'),
				'-' : lambda x: x.replace('nosys', 'pu_down'),
				},
			'JER' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad', 'single_top'],
				'categories' : ['.*'],
				'type' : 'shape',
				'+' : lambda x: x.replace('nosys', 'jer_up'),
				'-' : lambda x: x.replace('nosys', 'jer_down'),
				'constants' : ('jer_down', 'jer_up'),
				'value' : 1.00,				
				},
			'MTOP' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
				'categories' : ['.*'],
				'type' : 'shape',
				'+' : lambda x: x.replace('nosys', 'mtop_up'),
				#'-' : lambda x: x.replace('nosys', 'mtop_down'),
				'value' : 1.00,				
				},
			'HScale' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
				'categories' : ['.*'],
				'type' : 'shape',
				'+' : lambda x: x.replace('nosys', 'hadscale_up'),
				#'-' : lambda x: x.replace('nosys', 'hadscale_down'),
				'constants' : ('hadscale_down', 'hadscale_up'),
				'value' : 1.00,				
				},
			'BTAG' : {
				'samples' : ['.*'],
				'categories' : ['.*'],
				'type' : 'lnN',
				'+' : lambda x: x.replace('nosys', 'btag_up'),
				'-' : lambda x: x.replace('nosys', 'btag_down'),
				},
			'CTAGB' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'single_top', 'vjets', 'qcd'],
				'categories' : ['.*'],
				'type' : 'lnN',
				'+' : lambda x: x.replace('nosys', 'btagb_up'),
				'-' : lambda x: x.replace('nosys', 'btagb_down'),
				},
			'CTAGL' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'single_top', 'qcd', 'vjets'],#, 'right_whad'],
				'categories' : ['.*'],
				'type' : 'lnN',
				'+' : lambda x: x.replace('nosys', 'btagl_up'),
				'-' : lambda x: x.replace('nosys', 'btagl_down'),
				},
			'CTAGC' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'single_top', 'qcd', 'vjets'],#, 'right_whad'],
				'categories' : ['.*'],
				'type' : 'lnN',
				'+' : lambda x: x.replace('nosys', 'btagc_up'),
				'-' : lambda x: x.replace('nosys', 'btagc_down'),
				},
			'PDF' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
				'categories' : ['.*'],
				'type' : 'pdf',
				},
			}
		self.card = None
		self.binning = {
			}

	def save_card(self, name):
		if not self.card:
			raise RuntimeError('There is no card to save!')
		for sys_name, info in self.systematics.iteritems():
			if 'value' not in info: continue
			for category in info['categories']:
				for sample in info['samples']:
					plotter.card.add_systematic(
						sys_name, 
						info['type'], 
						category, 
						sample, 
						info['value']
						)
		self.card.save(name, self.outputdir)
		self.card = None
	
	def write_shapes(self, category_name, folder, variable,
									 rebin=1, preprocess=None):
		'''
		should use get_shape above in some memoized form in the future
		'''
		if not self.card: self.card = DataCard('TT')
		self.card.add_category(category_name)
		category = self.card[category_name]

		path = os.path.join(folder, variable)
		card_views = {}
		for name, samples in self.card_names.iteritems():
			card_views[name] = self.rebin_view(
				views.SumView(
					*[self.get_view(i) for i in samples]
					),
				rebin
				)

		if preprocess:
			for name in card_views:
				card_views[name] = preprocess(card_views[name])

		for name, view in card_views.iteritems():
			try:
				histo = view.Get(path)
			except:
				set_trace()
			integral = histo.Integral()
			
			category[name] = histo
			if name == 'data_obs': continue #skip systematics for data!
			
			#
			# shape and dynamically assigned systematics
			#
			for sys_name, info in self.systematics.iteritems():
				if not any(re.match(i, category_name) for i in info['categories']): continue
				if not any(re.match(i, name) for i in info['samples']): continue
				shift = 1.0
				if info['type'] == 'pdf': #pdf special case
					samples = self.card_names[name]
					histos = [self.pdf_unc_histo(i, j, 'mass_discriminant', rebin, 'nnpdf') for i in samples for j in folders]
					hpdf = sum(histos) if len(histos) > 1 else histos[0]
					hu = category[name].Clone()
					hd = category[name].Clone()
					for i,j,k in zip(hpdf, hu, hd):
						err = (i.error/i.value)*j.value if i.value else 0.
						j.value += err
						k.value -= err
					category['%s_%sUp'	% (name, sys_name)] = hu
					category['%s_%sDown' % (name, sys_name)] = hd					
					plotter.card.add_systematic(sys_name, 'shape', category_name, name, 1.00)
				elif info['type'] == 'shape' or 'value' not in info:
					path_up = info['+'](path) if '+' in info else None
					path_dw = info['-'](path) if '-' in info else None
					hup = view.Get(path_up) if path_up else None
					hdw = view.Get(path_dw) if path_dw else None

					if hup is None and hdw is None:
						raise RuntimeError('%s systematic does not define neither "+" nor "-" values' % sys_name)
					elif hup is None or hdw is None: 
						mirrored = histo.Clone()
						multiplier = integral 
						mirrored.Reset()
						src = hup if hup is not None else hdw
						for mbin, sbin, dbin in zip(mirrored, src, histo):
							delta = sbin.value - multiplier*dbin.value
							mbin.value = multiplier*dbin.value - delta
						if hup is None:
							hup = mirrored
						else:
							hdw = mirrored
					
					if info['type'] == 'lnN':
						if integral <= 0: continue
						upi = hup.Integral()
						dwi = hdw.Integral()
						rel_u = abs(upi-integral)/integral if integral else 0
						rel_d = abs(dwi-integral)/integral if integral else 0
						#rootpy.log["/"].info("Sys: %s %s/%s: %f -- %f" % (sys_name, category_name, name, rel_u, rel_d))
						sign =  1 if (integral - dwi) >= 0 else -1.
						delta = max(min(max(rel_u, rel_d), 1), -1)
						if (integral - dwi) >= 0:
							value = 1.00+delta
						else:
							value = 1. / (1.+delta)

						plotter.card.add_systematic(
							sys_name, info['type'],
							category_name, name, value
						)

					#shapes only: store shape in root file					
					if info['type'] == 'shape':
						if hup.Integral() == 0.:
							rootpy.log["/"].warning('%s Up for %s/%s has normalization == 0, forcing it to 10**-6' %(sys_name, category_name, name))
							mbin = category[name].GetMaximumBin()
							hup[mbin].value = 10**-6
						if hdw.Integral() == 0.:
							rootpy.log["/"].warning('%s Down for %s/%s has normalization == 0, forcing it to 10**-6' %(sys_name, category_name, name))
							mbin = category[name].GetMaximumBin()
							hdw[mbin].value = 10**-6
						category['%s_%sUp'	% (name, sys_name)] = hup
						category['%s_%sDown' % (name, sys_name)] = hdw					


	def create_subsample(self, baseview, subdirs, title, color='#9999CC'):
		dirmap = {
			'' : views.SumView(
				*[views.SubdirectoryView(self.views[baseview]['view'], i) for i in subdirs]
				 )
			}
		for shift, view in self.tt_shifted.iteritems():
			dirmap[shift] = views.SumView(
				*[views.SubdirectoryView(self.views[view]['view'], '%s/nosys' % i) for i in subdirs]
				)
		
		return views.StyleView(
			views.TitleView(
				dirmap[''],#				urviews.MultifileView(**dirmap),
				title
				),
			fillcolor = color
			)

	def create_tt_subsample(self, subdirs, title, color='#9999CC'):
		return self.create_subsample('ttJets', subdirs, title, color=color)

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
		self.mc_samples = [
			'[WZ][WZ]', 'QCD*', '[WZ]Jets', 
			'single*', 'ttJets_preselection'
			]
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

	def make_sample_table(self, threshold=None, absolute=False, fname='yields.raw_txt'):
		stack = [i for i in self.keep if isinstance(i, plotting.HistStack)][0]
		names = [i.title for i in stack.hists]
		yields = [i.Integral() for i in stack.hists]
		def ratio(lst):
			tot = sum(lst)/100. if not absolute else 1.
			return [i/tot for i in lst] if tot else [0. for _ in lst]
		yields = ratio(yields)
		mlen = max(len(i) for i in names)
		format = '%'+str(mlen)+'s     %.1f%%\n'
		header = ('%'+str(mlen)+'s    frac') % 'sample'
		if threshold is not None:
			binid = stack.hists[0].xaxis.FindBin(threshold)
			fbin = stack.hists[0].nbins()+1
			less = ratio([i.Integral(1, binid-1) for i in stack.hists])
			above = ratio([i.Integral(binid, fbin) for i in stack.hists])
			yields = [i for i in zip(yields, less, above)]
			format = format.replace('\n', '    %.1f%%    %.1f%%\n')
			header += '    <%.0f %s     >%.0f %s' % (threshold, stack.hists[0].xaxis.title, threshold, stack.hists[0].xaxis.title)
		header += '\n'
		fullname = os.path.join(self.outputdir, fname)
		with open(fullname, 'w') as f:
			f.write(header)
			for i in zip(names, yields):
				info = i
				if threshold is not None:
					#repack info
					a, b = i
					c, d, e = b
					info = (a, c, d, e)
				f.write(format % info)

	def get_yields(self, thr):
		stack = [i for i in self.keep if isinstance(i, plotting.HistStack)][0]
		obs = self.keep[1]
		assert(obs.title == 'Observed')
		hists = stack.hists
		hists.append(obs)
		ret = []
		for hist in hists:
			name = hist.title
			binid = hist.xaxis.FindBin(thr)
			fbin = hist.nbins()+1
			loe = ROOT.Double()
			lov = hist.IntegralAndError(1, binid-1, loe)
			hie = ROOT.Double()
			hiv = hist.IntegralAndError(binid, fbin, hie)
			ret.append(
				(name, (lov, loe), (hiv, hie))
				)
		return ret

plotter = HTTPlotter(args.mode)

variables = [
  (False, "m_tt"	, "m(tt) (GeV)", 2, None, False),		
  (False, "tlep_ctstar"	, "cos #theta^{*}(tlep)", 2, None, False),		
  (False, "thad_ctstar"	, "cos #theta^{*}(thad)", 2, None, False),		
  (False, "cdelta_ld", "cos #delta(ld)", 2, None, False),		
  (False, "hframe_ctheta_d", "cos #theta(d-jet)", 2, None, False),		

  (False, "pt_thad"	, "p_{T}(t_{had}) (GeV)", 5, (0, 1000), False),		
  (False, "eta_thad"	, "#eta_{T}(t_{had}) (GeV)", 5, None, False),		
  (False, "pt_tlep"	, "p_{T}(t_{lep}) (GeV)", 5, (0, 1000), False),		
  (False, "eta_tlep"	, "#eta_{T}(t_{lep}) (GeV)", 5, None, False),		
  (False, "pt_tt"	, "p_{T}(tt) (GeV)", 5, (0, 1000), False),		
  (False, "eta_tt"	, "#eta_{T}(tt) (GeV)", 2, None, False),		
  (False, "full_discriminant_4j"	, "#lambda_{comb} 4 jets",   2, (5, 30), False),		
  (False, "full_discriminant_5j"	, "#lambda_{comb} 5 jets",   2, (5, 30), False),		
  (False, "full_discriminant_6Pj"	, "#lambda_{comb} > 6 jets", 2, (5, 30), False),		
]

preselection = [
	(False, "MT" , "MT",  10, (0, 300), False),
  (False, "njets"	 , "# of selected jets", range(13), None, False),
  (False, "jets_eta", "#eta(jet)", 5, None, False),
  (False, "jets_pt", "p_{T}(jet)", 10, (0, 300), False),
  (False, "lep_eta", "#eta(l)", 5, None, False),
  (False, "lep_pt"	, "p_{T}(l) (GeV)", 20, (0, 300), False),		
  (False, "nvtx", "# of reconstructed vertices", range(41), None, False),
  (False, "rho", "#rho", range(40), None, False),
	(False, "lep_iso", 'l rel Iso', 1, [0,1], False),
	(False, "lep_wp" , "electron wp", 1, None, False),
	(True	, "cMVA"    , "cMVA",  1, None, False),
	(True	, "cMVA_p11", "cMVA^{11}", 1, None, False),
	(True	, "qgtag"   , "quark-gluon tag",  1, None, False),
]

permutations = [
	(False, "bjets_pt", "p_{T} (b-jets)", 10, (0, 300), False),
	(False, "wjets_pt", "p_{T} (W-jets)", 10, (0, 200), False),
	(False, "full_discriminant", "#lambda_{C}", 1, (8, 30), False),
	(False, "nu_discriminant", "#chi^{2}_{#nu}", 1, None, False),
	(False, "mass_discriminant", "#lambda_{M}", 1, None, False),
	(False, "Wmasshad", "M(W_{h})", 2, (0,400), False),
	(False, "tmasshad", "M(t_{h})", 2, None, False),
	(False, "lbratio" , "p_{T}(l)/p_{T}(b)", 1, (0, 6), False),
	(False, "j2bratio", "p_{T}(j)/p_{T}(b)", 1, (0, 6), False),
]

jet_categories = ["3jets", "4jets", "5Pjets"]

#cut flow
if args.flow or args.all:
	flow = plotter.cut_flow()
	plotter.save('cut_flow')

if args.preselection or args.all:
	plotter.set_subdir('preselection')
	for logy, var, axis, rebin, x_range, leftside in preselection + permutations:
		plotter.make_preselection_plot(
			'nosys/preselection', var, sort=True,
			xaxis=axis, leftside=leftside, rebin=rebin, 
			show_ratio=True, ratio_range=0.5, xrange=x_range, logy=logy)		
		plotter.save(var)

if args.plots or args.all:
	vals = []
	for dirid in itertools.product(['looseNOTTight', 'tight'], ['MTLow', 'MTHigh']):
		tdir = '%s/%s' % dirid
		plotter.set_subdir(tdir)
		first = True
		for logy, var, axis, rebin, x_range, leftside in preselection+variables+permutations:
			plotter.plot_mc_vs_data(
				'nosys/%s' % tdir, var, sort=True,
				xaxis=axis, leftside=leftside, rebin=rebin,
				show_ratio=True, ratio_range=0.5, xrange=x_range,
				logy=logy)
			if first:
				first = False
				plotter.make_sample_table(threshold=50)
				plotter.make_sample_table(threshold=50, absolute=True, fname='yields_abs.raw_txt')
				vals.append(
					(tdir, plotter.get_yields(50))
					)
			plotter.save(var)

	jmap = {}
	for tdir, sams in vals:
		for sam, lo, hi in sams:
			if sam not in jmap:
				jmap[sam] = {}
			jmap[sam]['%s/%s' % (tdir, 'hi')] = hi
			jmap[sam]['%s/%s' % (tdir, 'lo')] = lo
	with open('yields.json', 'w') as out:
		out.write(prettyjson.dumps(jmap))

if args.shapes or args.all:
	for peak, dname in [('*', 'shapes'), ('Peak', 'shapes_peak'), ('Int', 'shapes_interference')]:
		plotter.set_subdir(dname)
		for mass in [400, 500, 600, 750]:
			histos = []
			for width, color in zip([5, 10, 25, 50], ['#f9a505', '#2aa198', '#0055ff', '#6666b3']):
				htt_view = plotter.get_view('HtoTT_M%d_%dpc_%s' % (mass, width, peak))
				histos.append(
					sum(
						htt_view.Get('%s/nosys/tight/MTHigh/m_tt' % i) \
							for i in ['right', 'matchable', 'unmatchable', 'noslep']
						)
					)
				histos[-1].Rebin(2)
				histos[-1].linecolor = color
				if histos[-1].Integral():
					histos[-1].Scale(1/abs(histos[-1].Integral()))
			plotter.overlay(histos, xtitle='m(tt) (GeV)', ytitle='a.u.', x_range=(0, 1300))
			plotter.add_legend(histos)
			plotter.save('M%s' % mass)

if args.btag:
	mc_default = plotter.mc_samples
	plotter.mc_samples = [
		'[WZ][WZ]', 'QCD*', '[WZ]Jets', 
		'single*', 'ttJets_preselection'
		]
	hists = [i.Get('jets_cMVA_WP') for i in plotter.mc_views(1, None, 'nosys/preselection')]
	plotter.mc_samples = mc_default
	ttb = hists[-1]
	bkg = sum(hists[:-1])
	qcd = [i for i in hists if i.title == 'QCD'][0]
	rest = sum(i for i in hists if i.title != 'QCD')
	sig_sqrt_bkg = bkg.Clone()
	sig_sqrt_sb  = bkg.Clone()
	qcd_part = bkg.Clone()
	ttfrac = bkg.Clone()
	labels = [None, 'None', 'Loose', 'Medium', 'Tight']
	for xstrt, ystrt in itertools.product(range(1,5), range(1,5)):
		s, b = 0, 0
		q, e = 0, 0
		sig_sqrt_sb.xaxis.SetBinLabel(xstrt, labels[xstrt])
		sig_sqrt_sb.yaxis.SetBinLabel(ystrt, labels[ystrt])
		sig_sqrt_bkg.xaxis.SetBinLabel(xstrt, labels[xstrt])
		sig_sqrt_bkg.yaxis.SetBinLabel(ystrt, labels[ystrt])
		qcd_part.xaxis.SetBinLabel(xstrt, labels[xstrt])
		qcd_part.yaxis.SetBinLabel(ystrt, labels[ystrt])
		if not ttb[xstrt, ystrt].value and not bkg[xstrt, ystrt].value: continue
		for x, y in itertools.product(range(xstrt,5), range(ystrt,5)):
			s += ttb[x,y].value
			b += bkg[x,y].value
			q += qcd[x,y].value
			e += rest[x,y].value
		print labels[xstrt], labels[ystrt], q, e, q / (q+e)
		sig_sqrt_bkg[xstrt, ystrt].value = s/math.sqrt(b)
		sig_sqrt_sb[ xstrt, ystrt].value = s/math.sqrt(s+b)
		ttfrac[xstrt, ystrt].value = s
		qcd_part[xstrt, ystrt].value = q / (q+e)
	ttfrac.Scale(1/ttb.Integral())
	plotter.set_subdir('preselection')
	plotter.plot(sig_sqrt_bkg, drawstyle='colz')
	plotter.save('sig_sqrt_bkg')
	plotter.plot(sig_sqrt_sb, drawstyle='colz')
	plotter.save('sig_sqrt_sb')
	plotter.plot(qcd_part, drawstyle='colz')
	plotter.save('qcd_contamination')
	plotter.plot(ttfrac, drawstyle='colz')
	plotter.save('fraction_tt')

if args.card:
	plotter.write_shapes(		
		'mujets' if args.mode == 'muons' else 'ejets',
		'nosys/tight/MTHigh', 'mtt_tlep_ctstar',
		rebin=1, preprocess=None)
	plotter.set_subdir()
	plotter.save_card(args.mode)

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
		
