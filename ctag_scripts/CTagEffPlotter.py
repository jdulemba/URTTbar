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

parser = ArgumentParser()
parser.add_argument('--wps', default='*',
                    help='choose the working points to use')
parser.add_argument('--plots', dest='plots', action='store_true',
                    help='make control plots')
parser.add_argument('--mceffs', dest='mceffs', action='store_true',
                    help='Dump MC-based efficiencies')
parser.add_argument('--systematics', dest='systematics', action='store_true',
                    help='make systematic plots')
parser.add_argument('--shapes', dest='shapes', action='store_true',
                    help='make shapes root files')
parser.add_argument('--noLightFit', dest='noLightFit', action='store_true',
                    help='set up fitting without light fitting')
parser.add_argument('--noBBB', action='store_true',
                    help='do not run bin-by-bin uncertainties')
parser.add_argument('--inclusive', action='store_true',
                    help='use inclusive categories')
parser.add_argument('--lumi', type=float, default=-1.,
                    help='force luminosity')
parser.add_argument('--pdfs', action='store_true', help='make plots for the PDF uncertainties')
parser.add_argument('--noPOIpropagation', action='store_true')
args = parser.parse_args()

def syscheck(cmd):
	out = os.system(cmd)
	if out == 0:
		return 0
	else:
		raise RuntimeError("command %s failed executing" % cmd)

class CTagPlotter(Plotter):
	def __init__(self, lumi=None):
		lumi = lumi if lumi > 0 else None
		#self.tt_to_use = 'ttJets' #original
		self.flavour_info = 'hadronflav'
		self.tt_shifted = {
			'hdamp_up'  : ['ttJetsDiLep_hdampUP',  'ttJetsHad_hdampUP',  'ttJetsSL_hdampUP'],
			'hdamp_down': ['ttJetsDiLep_hdampDOWN','ttJetsHad_hdampDOWN','ttJetsSL_hdampDOWN'],
			'mtop_up'   : ['ttJetsDiLep_mtopUP',   'ttJetsHad_mtopUP',   'ttJetsSL_mtopUP'],
			'mtop_down' : ['ttJetsDiLep_mtopDOWN', 'ttJetsHad_mtopDOWN', 'ttJetsSL_mtopDOWN'],
			'ps_up'     : ['ttJetsDiLep_psUP',     'ttJetsHad_psUP',     'ttJetsSL_psUP'],
			'ps_down'   : ['ttJetsDiLep_psDOWN',   'ttJetsHad_psDOWN',   'ttJetsSL_psDOWN'],
			#'isr_up'   : 'ttJets_isrup',
			#'isr_down' : 'ttJets_isrdown',
			#'fsr_up'   : 'ttJets_fsrup',
			#'fsr_down' : 'ttJets_fsrdown',
			}
		jobid = os.environ['jobid']

		files = filter(lambda x: 'SingleElectron' not in x, glob.glob('results/%s/ctag_eff/*.root' % jobid))

		logging.debug('files found %s' % files.__repr__())
		lumis = glob.glob('inputs/%s/*.lumi' % jobid)
		logging.debug('lumi files found %s' % lumis.__repr__())
		
		outdir= 'plots/%s/ctageff' % jobid

		self.tt_lhe_weights = prettyjson.loads(
			open('inputs/%s/ttJets.weights.json' % jobid).read()
		)

		super(CTagPlotter, self).__init__(
			files, lumis, outdir, styles, None, lumi,
			)
		self.defaults = {
			'watermark' : ['(13 TeV, 25ns)', True, self.views['data']['intlumi']]
			}
		self.jobid = jobid

		self.tt_to_use = 'ttJetsAll'
		self.views['ttJets_preselection'] = {
			'view' : self.create_tt_subsample(
				['semilep_visible_right'],
				't#bar{t}',
				'r'
				#'#6666b3' #purplish
				)
			}
		self.views['ttJets_sig'] = {
			'view' : self.create_tt_subsample(
				['semilep_visible_right', 'semilep_right_thad', 'semilep_right_whad'], 
				't#bar{t}, right W_{h}',
				'r'
				#'#6666b3' #purplish
				)
			}
		self.views['ttJets_bkg'] = {
			'view' : self.create_tt_subsample(
				['semilep_right_tlep', 'semilep_wrong'], 
				't#bar{t}, wrong W_{h}',
				'b'
				#'#ab5555' #burgundy
				)
			}
		self.views['ttJets_other'] = {
			'view' : self.create_tt_subsample(
				['other'], 
				'Other tt decay',
				'm',
				#'#668db3', #light bluish
				)
			}

		self.views['VJets'] = {
			'view' : views.SumView(
				self.get_view('WJets'),
				#self.get_view('W[1-4]Jets'),
				#self.get_view('W[1-5]Jets'),
				self.get_view('ZJets'),
				)
			}
		##self.views['SingleTOP'] = {
    ##  'view' : views.ScaleView(
    ##    self.get_view('singlet*'),
		##		2
    ##    )
    ##  }

		self.mc_samples = [
			'QCD*',
			#'tt[WZ]*',
			'VJets',
			#'WJets',
			#'ZJets',
			'singlet*',
			'ttJets_other',
			'ttJets_bkg',
			'ttJets_sig',		 
			]

		self.card_names = {
			'qcd' : ['QCD*'],
			#'vjets'		: ['ZJets', 'W[1-4]Jets'],
			#'vjets'		: ['ZJets', 'W[1-5]Jets'],
			'vjets'		: ['[WZ]Jets'],#, 'Z[1-3]Jets'],
			'right_whad' : ['ttJets_sig'],
			'wrong_whad' : ['ttJets_bkg'],
			'nonsemi_tt' : ['ttJets_other'],
			'single_top' : ['singlet*'],
			#'single_top' : ['SingleTOP'],
			'data_obs'	: ['data']
			}

		self.card_by_title = {
			'QCD' : 'qcd' ,
			'V + jets' : 'vjets' ,
			'single top' : 'single_top' ,
			'Other tt decay' : 'nonsemi_tt' ,
			't#bar{t}, wrong W_{h}' : 'wrong_whad', 
			't#bar{t}, right W_{h}' : 'right_whad', 
			#'ttV' : 'ttV',
			'Observed' : 'data_obs'
			}

		#set_trace()
		self.signal = 'right_whad'
		self.signal_yields = {}
		self.right_whad_yields = {}

		self.systematics = {
			'lumi' : {
				'type' : 'lnN',
				'samples' : ['.*'],
				'categories' : ['.*'],
				'value' : 1.062,
				},
			'lepEff' : {
				'type' : 'lnN',
				'samples' : ['.*'],
				'categories' : ['.*'],
				'value' : 1.02,
				},
			'ttxsec' : {
				'type' : 'lnN',
				'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
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
				'value' : 3.00,
				#'value' : 1.70,
				},	
			'vjetxsec' : {
				'type' : 'lnN',
				'samples' : ['vjets'],
				'categories' : ['.*'],
				'value' : 1.20,
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
			'QCDscaleMERenorm_TT' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
				'categories' : ['.*'],
				'type' : 'shape',
				'+' : lambda x: x.replace('nosys', 'renorm_up'),
				'-' : lambda x: x.replace('nosys', 'renorm_down'),
				'value' : 1.00,				
				'scales' : (1./self.tt_lhe_weights['3'], 1./self.tt_lhe_weights['6']),
				},
			'QCDscaleMEFactor_TT' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
				'categories' : ['.*'],
				'type' : 'shape',
				'+' : lambda x: x.replace('nosys', 'factor_up'),
				'-' : lambda x: x.replace('nosys', 'factor_down'),
				'value' : 1.00,				
				'scales' : (1./self.tt_lhe_weights['1'], 1./self.tt_lhe_weights['2']),
				},
			'MTOP' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
				'categories' : ['.*'],
				'type' : 'shape',
				'+' : lambda x: x.replace('nosys', 'mtop_up'),
				'-' : lambda x: x.replace('nosys', 'mtop_down'),
				'value' : 1.00,				
				},
			'HDAMP' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
				'categories' : ['.*'],
				'type' : 'shape',
				'+' : lambda x: x.replace('nosys', 'hdamp_up'),
				'-' : lambda x: x.replace('nosys', 'hdamp_down'),
				'value' : 1.00,				
				},
			'PS' : {
				'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
				'categories' : ['.*'],
				'type' : 'shape',
				'+' : lambda x: x.replace('nosys', 'ps_up'),
				'-' : lambda x: x.replace('nosys', 'ps_down'),
				'value' : 1.00,				
				},
			#'ISR' : {
			#	'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
			#	'categories' : ['.*'],
			#	'type' : 'shape',
			#	'+' : lambda x: x.replace('nosys', 'isr_up'),
			#	'-' : lambda x: x.replace('nosys', 'isr_down'),
			#	'constants' : (
			#		'isr_down', 
			#		'isr_up'
			#		),
			#	'value' : 1.00,				
			#	},
			#'FSR' : {
			#	'samples' : ['wrong_whad', 'nonsemi_tt', 'right_whad'],
			#	'categories' : ['.*'],
			#	'type' : 'shape',
			#	'+' : lambda x: x.replace('nosys', 'fsr_up'),
			#	'-' : lambda x: x.replace('nosys', 'fsr_down'),
			#	'constants' : (
			#		'fsr_down', 
			#		'fsr_up'
			#		),
			#	'value' : 1.00,				
			#	},
			'JES' : {
				'samples' : ['.*_whad', 'nonsemi_tt', 'single_top'],
				'categories' : ['.*'],
				'type' : 'shape',
				'+' : lambda x: x.replace('nosys', 'jes_up'),
				'-' : lambda x: x.replace('nosys', 'jes_down'),
				'constants' : ('jes_down', 'jes_up'),
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
				'samples' : ['wrong_whad', 'nonsemi_tt', 'vjets', 'single_top', 'qcd'],
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
				'samples' : ['wrong_whad', 'nonsemi_tt', 'vjets', 'single_top', 'qcd'],#, 'right_whad'],
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
			'notag' : {
				'notag'	: [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
				},
			'DeepCSVLoose' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag' 	: [6, 8, 10, 12, 20],	
				},
			'DeepCSVMedium' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag' 	: [6, 8, 10, 12, 20],
				},
			'DeepCSVTight' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag'  : [6, 8, 10, 12, 20],
				},			
			'DeepJetLoose' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag' 	: [6, 8, 10, 12, 20],	
				},
			'DeepJetMedium' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag' 	: [6, 8, 10, 12, 20],
				},
			'DeepJetTight' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag'  : [6, 8, 10, 12, 20],
				},			
			'DeepCSVctagLoose' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag' 	: [6, 8, 10, 12, 20],
				},
			'DeepCSVctagMedium' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag' 	: [6, 8, 10, 12, 20],
				},
			'DeepCSVctagTight' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag' 	: [6, 8, 10, 12, 20],
				},
			'DeepJetctagLoose' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag' 	: [6, 8, 10, 12, 20],
				},
			'DeepJetctagMedium' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag' 	: [6, 8, 10, 12, 20],
				},
			'DeepJetctagTight' : {
				'notag' 	: [6, 8, 10, 12, 20],
				'leadtag' : [6, 8, 10, 12, 20],
				'subtag'  : [6, 8, 10, 12, 20],
				'ditag' 	: [6, 8, 10, 12, 20],
				},
			}


	def create_tt_subsample(self, subdirs, title, color='#9999CC'):
		dirmap = {
			'' : views.SumView(
				*[views.SubdirectoryView(self.views[self.tt_to_use]['view'], i) for i in subdirs]
				 )
			}
		for shift, view in self.tt_shifted.iteritems():
			#set_trace()
			dirmap[shift] = views.SumView(
				    *[views.SumView(
				        *[views.SubdirectoryView(self.views[v]['view'], '%s/nosys' % i) for i in subdirs]
				    ) for v in view]
				)
		#set_trace()
		
		return views.StyleView(
			views.TitleView(
				urviews.MultifileView(**dirmap),
				title
				),
			fillcolor = color,
			)

	def add_systematics(self, nobbb=False):
		if not plotter.card:
			raise ValueError('The card is not defined!')
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
		bbbs = []
		if not nobbb:
			print 'adding bbb'
			#plotter.card.add_bbb_systematics('.*', '^(?!right)')
			bbbs.extend(plotter.card.add_bbb_systematics('.*', 'wrong_whad'))
			bbbs.extend(plotter.card.add_bbb_systematics('.*', 'nonsemi_tt'))
			bbbs.extend(plotter.card.add_bbb_systematics('.*', 'single_top')) 
			bbbs.extend(plotter.card.add_bbb_systematics('.*', 'vjets')) 
			for cname, integral in self.signal_yields.iteritems():
				bbbs.extend(plotter.card.add_bbb_systematics(cname, 'right_whad', multiplier=integral))
			#bbbs.extend(plotter.card.add_bbb_systematics('.*', 'qcd'))

			#plotter.card.add_bbb_systematics('.*', 'right_whad', 0.05, relative=False)
			#plotter.card.add_bbb_systematics('.*', '.*')  
		return bbbs

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
		#set_trace()
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
		#stack.Draw('text')
		#data.Draw('text')
		return cflow
		#cut_flow.GetYaxis().SetRangeUser(1, 10**7)

	def save_card(self, name):
		if not self.card:
			raise RuntimeError('There is no card to save!')
		self.card.save(name, self.outputdir)
		self.card = None

	def get_shape(self, folders, process, rebin=1, systematic='', preprocess=None):
		'''get_shape(self, folders, process, rebin=1, systematic='', preprocess=None)
		get a shape fit style systematic return the sys-shifted shape, the name should 
		be passed as 'NAME+'/'NAME-' to indicate the shift up/down
		'''
		paths = [os.path.join(i, 'mass_discriminant') for i in folders]
		view = self.rebin_view(
			views.SumView(*[self.get_view(i) for i in self.card_names[process]]),
			rebin
			)

		if preprocess:
			view = preprocess(view)

		if systematic:
			sys_name = systematic[:-1]
			shift = systematic[-1]
			paths = [self.systematics[sys_name][shift](path) for path in paths] 
			
		return sum(view.Get(path) for path in paths)
		

	def write_mass_discriminant_shapes(self, category_name, folders, rebin=1, 
												  preprocess=None):
		'''
		should use get_shape above in some memoized form in the future
		'''
		if not self.card: self.card = DataCard(self.signal)
		self.card.add_category(category_name)
		category = self.card[category_name]

		paths = [os.path.join(i, 'mass_discriminant') for i in folders]
		card_views = {}
		for name, samples  in self.card_names.iteritems():
			card_views[name] = self.rebin_view(
				views.SumView(
					*[self.get_view(i) for i in samples]
					),
				rebin
				)

		if preprocess:
			for name in card_views:
				card_views[name] = preprocess(card_views[name])

		data_is_fake = False
		if 'data_obs' not in card_views:
			logging.error('Using fake data given real one are missing!')
			card_views['data_obs'] = views.SumView(*card_views.values())
			data_is_fake = True

		for name, view in card_views.iteritems():
			histo = sum(view.Get(path) for path in paths)
			integral = histo.Integral()
			#if args.eras=='B':
			#	era_title = 'Run B'
			#elif args.eras=='CtoE':
			#	era_title = 'Runs C-E'
			#elif args.eras=='EtoF':
			#	era_title = 'Runs E-F'
			#elif args.eras == 'All':
			#	era_title = 'All 2017 (B-F)'
			#self.views['data']['view'] = views.TitleView(self.views['data']['view'], era_title)
			##set_trace()
			if name == self.signal:
				histo.Scale(1./integral)
				self.signal_yields[category_name] = integral
				#set_trace()
			elif name == 'data_obs' and data_is_fake:
				int_int = float(int(integral))
				histo.Scale(int_int/integral)
			
			category[name] = histo
			if name == 'data_obs': continue #skip systematics for data!
			
			if name == 'right_whad':
				self.right_whad_yields[category_name] = integral
				#set_trace()
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
					paths_up = [info['+'](path) for path in paths] if '+' in info else None
					paths_dw = [info['-'](path) for path in paths] if '-' in info else None
					hup = sum(view.Get(i) for i in paths_up) if paths_up else None
					hdw = sum(view.Get(i) for i in paths_dw) if paths_dw else None
					if 'scales' in info:
						sup, sdw = info['scales']
						hup.Scale(sup)
						hdw.Scale(sdw)

					if hup is None and hdw is None:
						raise RuntimeError('%s systematic does not define neither "+" nor "-" values' % sys_name)
					elif hup is None or hdw is None: 
						mirrored = histo.Clone()
						multiplier = integral if name == self.signal else 1.
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
						#if name == self.signal: set_trace()

					#shapes only: store shape in root file					
					if info['type'] == 'shape':
						if name == self.signal:
							hup.Scale(1./integral)
							hdw.Scale(1./integral)
						if hup.Integral() == 0.:
							rootpy.log["/"].warning('%s Up for %s/%s has normalization == 0, forcing it to 10**-6' %(sys_name, category_name, name))
							mbin = category[name].GetMaximumBin()
							hup[mbin].value = 10**-6
						if hdw.Integral() == 0.:
							rootpy.log["/"].warning('%s Down for %s/%s has normalization == 0, forcing it to 10**-6' %(sys_name, category_name, name))
							mbin = category[name].GetMaximumBin()
							hdw[mbin].value = 10**-6
						#set_trace()
						category['%s_%sUp'	% (name, sys_name)] = hup
						category['%s_%sDown' % (name, sys_name)] = hdw					
					
		return self.right_whad_yields
		#set_trace()

	def write_summary_table(self, wpoints, shift='nosys'):
		right_wpoints = [i for i in wpoints if i <> 'notag']
		rightW_view = self.get_view('ttJets_sig') 
		#mc_weight = self.views[self.tt_to_use]['weight'] ## original

		info = {}
		pflav_path = '{systematic}/mass_discriminant/{wpoint}/{jtag}/{jrank}/%s' % self.flavour_info
		def get_info(sys, wpoint, jtag, jrank):
			hpath = pflav_path.format(
				systematic=sys, 
				wpoint=wpoint,
				jtag=jtag,
				jrank=jrank
				)
			histo = rightW_view.Get(hpath)
			ret = {}
			err = ROOT.Double()
			tot = histo.IntegralAndError(1, histo.GetNbinsX(), err)
			ret['total'] = ufloat(tot, err)
			c_bin = histo.FindBin(4) #should be 5
			assert(c_bin == 5)
			ret['charm'] = ufloat(histo[c_bin].value, histo[c_bin].error)
			light, lerr = 0, 0
			for idx in range(1, c_bin):
				light += histo[idx].value
				lerr = quad.quad(lerr, histo[idx].error)
			ret['light'] = ufloat(light,lerr)
			return ret		

		def cjet_fraction(sys, jrank):
			info = get_info(sys, 'notag', 'both_untagged', jrank)
			return info['charm']/info['total']
			
		def eff_val(n_passing, n_total):
			#print 'passing nominal: ', n_passing.nominal_value
			#print 'total nominal: ', n_total.nominal_value
			#print 'sq err: ', n_passing.nominal_value * (1 - n_passing.nominal_value / n_total.nominal_value)
			val = n_passing.nominal_value / n_total.nominal_value
			err = math.sqrt(n_passing.nominal_value * (1 - n_passing.nominal_value / n_total.nominal_value)) / n_total.nominal_value			
			return ufloat(val, err)

		def effs(sys, wp, jrank):
			notag = get_info(sys, wp, 'both_untagged', jrank)
			ltag = get_info(sys, wp, 'lead_tagged', jrank)
			stag = get_info(sys, wp, 'sublead_tagged', jrank)
			dtag = get_info(sys, wp, 'both_tagged', jrank)

			def pass_and_tot(flav, jrank):
				ctot = sum(i[flav] for i in [notag, ltag, stag, dtag])
				cpass = dtag[flav]
				if jrank == 'leading': cpass += ltag[flav]
				else: cpass += stag[flav]				
				return cpass, ctot
			c_eff = eff_val(*pass_and_tot('charm', jrank))
			l_eff = eff_val(*pass_and_tot('light', jrank))
			return c_eff, l_eff

		info['lcfrac'] = cjet_fraction(shift, 'leading')
		info['scfrac'] = cjet_fraction(shift, 'subleading')
		info['total'] = get_info(shift, 'notag', 'both_untagged', 'leading')['total']
		for wp in right_wpoints:
			lceff, lleff = effs(shift, wp, 'leading')
			sceff, sleff = effs(shift, wp, 'subleading')
			info[wp] = {
				'lead' : {
					'leff' : lleff,
					'ceff' : lceff,
					},
				'sub' : {
					'leff' : sleff,
					'ceff' : sceff
					}
				}
		return info 

	def plot_mc_shapes(self, folder, variable, rebin=1, xaxis='',
							 leftside=True, xrange=None, preprocess=None,
							 normalize=False, logy=False, filt=None, show_err=False):
		mc_views = None
		if not filt:
			mc_views = self.mc_views(rebin, preprocess, folder)
		else:
			mc_views = [views.SubdirectoryView(
					self.rebin_view(
						self.get_view(i), 
						rebin
						), 
					folder
					) for i in filt]
			if preprocess:
				mc_views = [preprocess(i) for i in mc_views]

		histos = [i.Get(variable) for i in mc_views]

		if normalize:
			for i in histos:
				if i.Integral():
					i.Scale(1./i.Integral())
		m = 0.
		first = True
		for hist in histos:
			hist.fillstyle = 'hollow'
			hist.legendstyle = 'l'
			hist.linewidth = 2
			hist.markercolor = hist.linecolor
			if show_err:
				hist.drawstyle = 'pe'
				hist.markerstyle = 20
				hist.legendstyle = 'pe'

			hist.GetXaxis().SetTitle(xaxis)
			if first:
				hist.Draw()
				first = False
			else:
				hist.Draw('same')
			
			self.keep.append(hist)
			m = max(m,max(bin.value for bin in hist))

		histos[0].GetYaxis().SetRangeUser(0., m*1.2)
		if xrange:
			histos[0].GetXaxis().SetRangeUser(*tuple(xrange))
		self.add_legend(histos, leftside)

		samples = filt if filt else self.mc_samples
		try:
			ref_idx = samples.index('ttJets_allRight') 
			ytitle='bkg / signal'
		except ValueError:
			ref_idx = 0
			ytitle='other / %s' % histos[0].GetTitle()
		try:
			ref = histos[ref_idx]
		except:
			set_trace()
		tests = histos[:ref_idx]+histos[ref_idx+1:]

		self.add_ratio_plot(
			ref, *tests, x_range=xrange,
			ratio_range=3, ytitle=ytitle)

	def make_flavor_table(self, histo, fname='', to_json=False, to_txt=True, hadron_flavour=False):
		ids = {
			0 : 'pile-up',
			1 : 'down',
			2 : 'up',
			3 : 'strange',
			4 : 'charm',
			5 : 'beauty',
			6 : 'top',
			21: 'gluon',				
			}
		yield_names = {
			'light' : set([1,2,3]), 
			'charm' : set([4]), 
			'other' : set([0, 5, 6, 21]),
			'strange' : set([3]),
			}
		if hadron_flavour:
			ids = {
				0 : 'light',
				4 : 'charm',
				5 : 'beauty',
				}
			yield_names = {
				'light' : set([0]),
				'charm' : set([4]),
				'other' : set([5]),
				}
		format = '%10s %10s +/- %10s %10s\n'
		yields = {
			'light' : 0, 'charm' : 0, 'other' : 0, 'strange' : 0,
			'light_err' : 0, 'charm_err' : 0, 'other_err' : 0, 'strange_err' : 0,
					 }
		quark_yields = {}
		quark_errors = {}

		total = histo.Integral()
		for idx, name in ids.iteritems():
			entries  = histo.GetBinContent(histo.FindBin(idx))
			#entries += histo.GetBinContent(histo.FindBin(-idx))
			error = histo.GetBinError(histo.FindBin(idx))
			quark_yields[name] = entries
			quark_errors[name] = error
			for name, idset in yield_names.iteritems():
				if idx in idset: 
					yields[name] += entries

		if to_txt:
			fullname = os.path.join(self.outputdir, fname)
			with open(fullname, 'w') as f:
				f.write('Total events: %.0f\n' % total)
				header = format % ('quarks', 'yield', 'error', 'relative')
				f.write(header)
				f.write('-'*len(header)+'\n')
				for _, name in ids.iteritems():
					entries = quark_yields[name]
					error = quark_errors[name]
					ratio = 0 if total == 0 else (entries/total)*100
					f.write(format % (name, '%.1f' % entries, '%.1f' % error, '%.f%%' % (ratio)))
			logging.info('Wrote to %s' % fullname)
		if to_json:
			with open(os.path.join(self.outputdir, fname.replace('.raw_txt','.json')), 'w') as f:
				f.write(prettyjson.dumps(yields))

		return (quark_yields, quark_errors), yields

	def bake_pie(self, path):
		var = '%s/%s' % (path, self.flavour_info)
		mc_views = self.mc_views()
		mc_histos = [i.Get(var) for i in mc_views]
		integrals = [i.Integral() for i in mc_histos]
		data = self.get_view('data').Get(var).Integral()		
		total = sum(integrals)
		ratios = [i/total for i in integrals]
		names = [i.GetTitle() for i in mc_histos]
		colors = [i.GetFillColor('root') for i in mc_histos]		
		format = '%20s %10s %10s\n'
		entries = [('observed', '%d' % data)] + [(i, '%d' % j) for i, j in zip(names, integrals)]
		with open(os.path.join(self.outputdir, 'yields.tex'), 'w') as f:
			f.write(latex_table(['process', 'yield'], entries))
		with open(os.path.join(self.outputdir, 'yields.raw_txt'), 'w') as f:
			f.write('Total events: %.0f\n' % total)
			header = format % ('process', 'yield', 'relative')
			f.write(header)
			f.write('-'*len(header)+'\n')
			for name, entries, ratio in zip(names, integrals, ratios):
				f.write(format % (name, '%.0f' % entries, '%.0f%%' % (ratio*100)))

		hsum = sum(mc_histos)
		self.make_flavor_table(hsum, 'flavors.raw_txt', to_json=True)

		right_cmb = self.get_view('ttJets_sig').Get(var)
		self.make_flavor_table(right_cmb, 'flavors_sig.raw_txt')
		## 
		## wright_view = views.SumView(
		##	 *[self.get_view(i) for i in ['ttJets_rightHad', 'ttJets_rightWHad']]
		##	  )
		## wright = wright_view.Get(var) + right_cmb
		## self.make_flavor_table(wright, 'flavors_rightw.raw_txt', to_json=True)

	def make_preselection_plot(self, *args, **kwargs):
		systematics = None
		if 'sys_effs' in kwargs:
			systematics = kwargs['sys_effs']
			del kwargs['sys_effs']
		mc_default = self.mc_samples
		self.mc_samples = ['QCD*', 'VJets*', 'single*', 'ttJets_preselection']
		#self.mc_samples = ['QCD*', 'ZJets*', 'W[1-4]Jets*', 'single*', 'ttJets_preselection']
		#self.mc_samples = ['QCD*', '[WZ]Jets', 'single*', 'ttJets_preselection']
		##self.mc_samples = ['[WZ]Jets*', 'single*', 'ttJets_preselection']
		#set_trace()
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
			#set_trace()
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
			self.keep[0].SetMinimum(10**-3)
			self.canvas.Update()
			# Add legend
			self.pad.cd()
			self.add_legend(
				[mc_stack, stack_sum, data], kwargs.get('leftside', True), 
				entries=len(mc_stack.hists)+2
				)

		self.mc_samples = mc_default

	def pdf_unc_histo(self, sample, folder, var, rebin, pdf):
		view = self.rebin_view(
			views.SubdirectoryView(
				views.StyleView(
					self.get_view(sample), 
					drawstyle='E2',
					legendstyle='f'
					),
				folder
				),
			rebin
			)
		nnpdf = range(0, 107) ## changed for 108 pdf hists kept from Fall17 MC
		#nnpdf = [1] + range(10, 110)
		ct10 = range(112, 165)
		mmht = range(167, 218)
		vartemplate = '%s_mcws_pdf_%d'

		def rms_envelope(view, var, idrange):
			central = view.Get(vartemplate % (var, idrange[0]))
			sqsum = central.Clone()
			sqsum.Reset()
			for idx in idrange[1:]:
				histo = view.Get(vartemplate % (var, idx))
				for sbin, cbin, hbin in zip(sqsum, central, histo):
					sbin.value += (cbin.value-hbin.value)**2

			nreplicas = len(idrange)-1
			for sbin, cbin in zip(sqsum, central):
				err = math.sqrt(sbin.value/nreplicas)
				cbin.error = err
			return central
					
		def quad_envelope(view, var, idrange):
			central = view.Get(vartemplate % (var, idrange[0]))
			sqsum = central.Clone()
			sqsum.Reset()
			for pos, idx in enumerate(idrange[1:]):
				if pos % 2: continue
				h1 = view.Get(vartemplate % (var, idx))
				h2 = view.Get(vartemplate % (var, idrange[pos+1]))
				for sbin, cbin, h1bin, h2bin in zip(sqsum, central, h1, h2):
					d1 = abs(cbin.value-h1bin.value)
					d2 = abs(cbin.value-h2bin.value)					
					sbin.value = quad.quad(sbin.value, max(d1,d2))

			for sbin, cbin in zip(sqsum, central):				
				cbin.error = sbin.value
			return central

		if pdf == 'nnpdf':
			return rms_envelope( view, var, nnpdf)
		elif pdf == 'ct10':
			return quad_envelope(view, var, ct10)
		elif pdf == 'mmht':
			return quad_envelope(view, var, mmht)
		else:
			raise RuntimeError('invalid pdf type')

	def pdf_unc_plots(self, sample, folder, var, rebin):
		h_nnpdf = self.pdf_unc_histo(sample, folder, var, rebin, 'nnpdf')
		h_nnpdf.xaxis.title = '#lambda_{M}'
		h_nnpdf.yaxis.title = 'events'
		h_nnpdf.title = 'nnpdf'
		h_nnpdf.markercolor = '#22528b'
		h_nnpdf.linecolor = '#22528b'
		h_nnpdf.fillcolor = '#22528b'
		h_nnpdf.fillstyle = '//'
		h_ct10  = self.pdf_unc_histo(sample, folder, var, rebin, 'ct10')
		h_ct10.xaxis.title = '#lambda_{M}'
		h_ct10.yaxis.title = 'events'
		h_ct10.title = 'ct10'
		h_ct10.markercolor = '#a32020'
		h_ct10.linecolor = '#a32020'
		h_ct10.fillcolor = '#a32020'
		h_ct10.fillstyle = '\\\\'
		h_mmht  = self.pdf_unc_histo(sample, folder, var, rebin, 'mmht')
		h_mmht.xaxis.title = '#lambda_{M}'
		h_mmht.yaxis.title = 'events'
		h_mmht.title = 'mmht'
		h_mmht.markercolor = '#006d67'
		h_mmht.linecolor = '#006d67'
		h_mmht.fillcolor = '#006d67'
		h_mmht.fillstyle = '|'
		self.overlay([h_nnpdf, h_ct10, h_mmht])
		self.add_legend([h_nnpdf, h_ct10, h_mmht], False)


	def draw_flavor_shapes(self, dirname, basename, leftleg, preselection):
		if preselection:
			mc_default = self.mc_samples
			#self.mc_samples = ['QCD*', 'ZJets*', 'W[1-4]Jets*', 'single*', 'ttJets_preselection']
			self.mc_samples = ['QCD*', 'ZJets*', 'WJets*', 'single*', 'ttJets_preselection']
			#self.mc_samples = ['QCD*', '[WZ]Jets', 'single*', 'ttJets_preselection']
			#self.mc_samples = ['[WZ]Jets*', 'single*', 'ttJets_preselection']
			##self.mc_samples = ['[WZ]Jets*', 'Z[1-3]Jets', 'single*', 'ttJets_preselection']
		mc_views = plotter.mc_views(folder=dirname)	
		if preselection: self.mc_samples = mc_default

		b_hists = [i.Get('%s_B' % basename) for i in mc_views]
		c_hists = [i.Get('%s_C' % basename) for i in mc_views]
		l_hists = [i.Get('%s_L' % basename) for i in mc_views]
		
		hb = sum(b_hists)
		hl = sum(l_hists)
		hc = sum(c_hists)
		
		if not hasattr(self, 'flav_styles'):
			self.flav_styles = {
				#'b' : {'fillcolor' : '#0055ff', 'linecolor' : '#0055ff'},
				'b' : {'fillcolor' : 'r', 'linecolor' : 'k', 'markercolor' : 'k'},
				#'c' : {'fillcolor' : '#9999CC', 'linecolor' : '#9999CC'},
				'c' : {'fillcolor' : '#008900', 'linecolor' : 'k', 'markercolor' : 'k'},
				#'l' : {'fillcolor' : '#ab5555', 'linecolor' : '#ab5555'},
				'l' : {'fillcolor' : 'b', 'linecolor' : 'k', 'markercolor' : 'k'},
				#'s' : {'fillcolor' : '#FFCC66', 'linecolor' : '#FFCC66'}
				's' : {'fillcolor' : '#FFCC66', 'linecolor' : 'k', 'markercolor' : 'k'}
				}
		hb.decorate(**self.flav_styles['b']) 
		hb.title = 'b-jets'
		hl.decorate(**self.flav_styles['l'])
		hl.title = 'l-jets'
		hc.decorate(**self.flav_styles['c'])
		hc.title = 'c-jets'
		stack = plotter.create_stack(hb, hl, hc, sort=False)
		#data = plotter.get_view('data').Get(dirname+'%s_' % basename)
		data = sum([ plotter.get_view('data').Get(dirname+'%s_%s' % (basename, quark)) for quark in ['B', 'C', 'L'] ])
		plotter.overlay_and_compare(
            [stack], data, 
            method='datamc',
            ignore_style=False,
            lower_y_range=0.5,
            xtitle='b-tag SF',
            ytitle='events',
            #ytitle='data/MC',
            )
		self.canvas.Update()
		self.pad.cd()
		plotter.add_legend([stack, data], leftleg, 5)
		plotter.save('%s_flavours' % basename)
		
		#hb.Scale(1/hb.Integral() if hb.Integral() else 0.)
		#hl.Scale(1/hl.Integral() if hl.Integral() else 0.)
		#hc.Scale(1/hc.Integral() if hc.Integral() else 0.)
		#plotter.overlay([hb,hl,hc], fillstyle='hollow', linewidth=3, linecolor=fillcolor, xtitle='b-tag SF', ytitle='Events')
		##plotter.overlay([hb,hl,hc], fillstyle='hollow', linewidth=3, xtitle='b-tag SF', ytitle='Events')
		#plotter.add_legend([hb,hl,hc], leftleg, 4)
		#plotter.save('%s_norm' % basename)
		#set_trace()


	def draw_cvsl_shapes(self, dirname, basename, disc, leftleg, preselection):
		if preselection:
			mc_default = self.mc_samples
			#self.mc_samples = ['QCD*', 'ZJets*', 'W[1-4]Jets*', 'single*', 'ttJets_preselection']
			self.mc_samples = ['QCD*', 'ZJets*', 'WJets*', 'single*', 'ttJets_preselection']
			#self.mc_samples = ['QCD*', '[WZ]Jets', 'single*', 'ttJets_preselection']
			#self.mc_samples = ['[WZ]Jets*', 'single*', 'ttJets_preselection']
			##self.mc_samples = ['[WZ]Jets*', 'Z[1-3]Jets', 'single*', 'ttJets_preselection']
		mc_views = plotter.mc_views(folder=dirname)	
		if preselection: self.mc_samples = mc_default

		b_hists = [i.Get('%s_%s_B' % (basename, disc)) for i in mc_views]
		c_hists = [i.Get('%s_%s_C' % (basename, disc)) for i in mc_views]
		l_hists = [i.Get('%s_%s_L' % (basename, disc)) for i in mc_views]
		
		hb = sum(b_hists)
		hl = sum(l_hists)
		hc = sum(c_hists)
		
		if not hasattr(self, 'flav_styles'):
			self.flav_styles = {
				#'b' : {'fillcolor' : '#0055ff', 'linecolor' : '#0055ff'},
				'b' : {'fillcolor' : 'r', 'linecolor' : 'k', 'markercolor' : 'k'},
				#'c' : {'fillcolor' : '#9999CC', 'linecolor' : '#9999CC'},
				'c' : {'fillcolor' : '#008900', 'linecolor' : 'k', 'markercolor' : 'k'},
				#'l' : {'fillcolor' : '#ab5555', 'linecolor' : '#ab5555'},
				'l' : {'fillcolor' : 'b', 'linecolor' : 'k', 'markercolor' : 'k'},
				#'s' : {'fillcolor' : '#FFCC66', 'linecolor' : '#FFCC66'}
				's' : {'fillcolor' : '#FFCC66', 'linecolor' : 'k', 'markercolor' : 'k'}
				}
		hb.decorate(**self.flav_styles['b']) 
		hb.title = 'b-jets'
		hl.decorate(**self.flav_styles['l'])
		hl.title = 'l-jets'
		hc.decorate(**self.flav_styles['c'])
		hc.title = 'c-jets'
		stack = plotter.create_stack(hb, hl, hc, sort=False)
		#print dirname+'%s_%s' % (basename, disc.replace('hflav_','').replace('pflav_',''))
		data = plotter.get_view('data').Get(dirname+'%s_%s' % (basename, disc.replace('hflav_','').replace('pflav_','')))
		#set_trace()
		##plotter.overlay([stack, data], xtitle=disc.replace('hflav_','')+ ' (%s)' % basename, ytitle='events')
		plotter.overlay_and_compare(
            [stack], data, 
            method='datamc',
            ignore_style=False,
            lower_y_range=0.5,
            xtitle=disc.replace('hflav_','')+ ' (%s)' % basename,
            ytitle='events',
            #ytitle='data/MC',
            )
		self.canvas.Update()
		self.pad.cd()
		plotter.add_legend([stack, data], leftleg, 5)
		plotter.save('%s_%s_flavour' % (basename, disc))
		
		hb.Scale(1/hb.Integral() if hb.Integral() else 0.)
		hl.Scale(1/hl.Integral() if hl.Integral() else 0.)
		hc.Scale(1/hc.Integral() if hc.Integral() else 0.)
		plotter.overlay([hb,hl,hc], fillstyle='hollow', linewidth=3, xtitle=disc.replace('hflav_','')+ ' (%s)' % basename, ytitle='Events')
		plotter.add_legend([hb,hl,hc], leftleg, 4)
		plotter.save('%s_%s_norm' % (basename, disc))
		#print '\n\n%s_%s_norm\n\n' % (basename, disc)

			#plotter.overlay_and_compare(
			#	[mc_stack, stack_sum], data,
			#	xtitle = kwargs.get('xaxis',''),
			#	ytitle='Events', ignore_style=True,				
			#	method='datamc'
			#	)

	def get_efficiency(self, sample, wp, ptcat, flav):
		sample_view = views.SubdirectoryView(
			self.get_view(sample),
			'nosys/mass_discriminant/%s' % wp
			)
		project = 'ProjectionX' if ptcat == 'leading' else 'ProjectionY'
		ntag = rootpy.asrootpy( getattr(sample_view.Get('both_untagged/flav_map') , project)() )
		ltag = rootpy.asrootpy( getattr(sample_view.Get('lead_tagged/flav_map')	, project)() )
		stag = rootpy.asrootpy( getattr(sample_view.Get('sublead_tagged/flav_map'), project)() )
		dtag = rootpy.asrootpy( getattr(sample_view.Get('both_tagged/flav_map')	, project)() )
		all_jets = sum([ntag, ltag, stag, dtag])
		if ptcat == 'leading':
			pass_jets = ltag+dtag
		elif ptcat == 'subleading':
			pass_jets = stag+dtag
		else:
			raise ValueError('ptcat can only be (sub)leading')
		
		cname = self.card_by_title[pass_jets.title]
		eff = pass_jets[flav].value / all_jets[flav].value
		return cname, eff


plotter = CTagPlotter(args.lumi)

jet_variables = [
	('energy', 5, '%s jet energy (GeV)', None),
	('pt' ,  10, '%s jet p_{T} (GeV)', None),
	('eta',	1, '%s jet #eta', None),
	('phi',  10, '%s jet #varphi', None),
	#('abs_pflav_smart', 1, '%s jet parton flavour', None),#[-7, 22]),
	('hadronflav', 1, '%s jet parton flavour', None),#[-7, 22]),
	#("ncharged", 1, "%s jet charged constituents", None),
	#("nneutral", 1, "%s jet neutral constituents", None),
	#("ntotal"  , 1, "%s jet total constituents", None),
]

systematics_to_check = [
	'JES',
	'JER',
	'BTAG',
	'CTAGL', 
	'CTAGB', 
	'CTAGC', 
	'MTOP', 
	'FSR', 
	'ISR', 
	'pu'
]

vars2D = [
	# ('Whad_jet_pts', 'p_{T}(lead W jet) (GeV)', 'p_{T}(sub W jet) (GeV)', (1,1)),
]

variables = [
    ("njets"	 , "# of selected jets", range(13), None, False),
    #("evt_weight", "event weight", range(-3,4), None, False),
    #("btag_sf"	 , "SF applied to b-tagged jets", 1, None, False),
    #("muon_sf"	 , "SF applied to muon", 1, None, False),
    ("lep_pt"	, "p_{T}(l) (GeV)", 20, None, False),
    ("Whad_mass", "m_{W}(had) (GeV)", 1, None, False),
    ("thad_mass", "m_{t}(had) (GeV)", 1, None, False),
    ("thad_pt", "p_{T}(t_{h}) (GeV)", 1, None, False),
    #("mass_discriminant", "#lambda_{M}", 1, [0,15], False), #[5, 20]),
    ("mass_discriminant", "#lambda_{M}", 1, None, False), #[5, 20]),
    ("Wjets_DeepCSVb" , "DeepCSV Prob(b)", 1, None, False),
    ("Wjets_DeepCSVl" , "DeepCSV Prob(l) ", 1, None, False),
    ("Wjets_DeepCSVbb", "DeepCSV Prob(bb)", 1, None, False),
    ("Wjets_DeepCSVc" , "DeepCSV Prob(c)", 1, None, False),
    ("Wjets_DeepCSVbD", "DeepCSV b Disc", 1, None, False),
    ("Wjets_DeepCSVCvsB", "DeepCSV CvsB Disc (W Jet)", 1, None, False),
    ("Wjets_DeepCSVCvsL", "DeepCSV CvsL Disc (W Jet)", 1, None, False),
    #("Wlep_mass", "m_{W}(lep) (GeV)", 10, None, False),
    ("Wjets_DeepJetCvsB", "DeepJet CvsB Disc (W Jet)", 1, None, False),
    ("Wjets_DeepJetCvsL", "DeepJet CvsL Disc (W Jet)", 1, None, False),
    ("Wjets_DeepFlavourb" ,   "DeepJet Prob(b)", 1, None, False),
    ("Wjets_DeepFlavourbb",   "DeepJet Prob(bb)", 1, None, False),
    ("Wjets_DeepFlavourlepb" ,"DeepJet Prob(lepb) ", 1, None, False),
    ("Wjets_DeepFlavourc" ,   "DeepJet Prob(c)", 1, None, False),
    ("Wjets_DeepFlavouruds" , "DeepJet Prob(uds)", 1, None, False),
    ("Wjets_DeepFlavourg" ,   "DeepJet Prob(g)", 1, None, False),
    ("Wjets_DeepFlavourbD",   "DeepJet b Disc", 1, None, False),
    #("Whad_DR"  , "#DeltaR(jj) W_{had} (GeV)", 1, [0,7]),
    #("Whad_pt"  , "p_{T}(W_{had}) (GeV)", 10, None, False),
    #("Whad_leading_DR", "#DeltaR(W_{had}, leading jet)"	, 1, [0,7]),
    #("Whad_sublead_DR", "#DeltaR(W_{had}, subleading jet)", 1, [0,7]),
    #("nu_chisq"			, "nu_chisq"			, 1, [0, 20]),
    #("nu_discriminant"	, "nu_discriminant"	 , 1, None, False),
    #("btag_discriminant", "btag_discriminant", 1, [-11, -2]),
    #("nbjets"	, "# of bjets", 1, [0, 12]),
    #("lep_b_pt" , "p_{T}(b) (GeV)", 10, None, False),
    #("had_b_pt" , "p_{T}(b) (GeV)", 10, None, False),
]

preselection = [
    ("njets"	 , "# of selected jets", range(13), None, False),
    ("jets_DeepCSVb" , "DeepCSV Prob(b)" , 1, None, False),
    ("jets_DeepCSVl" , "DeepCSV Prob(l)" , 1, None, False),
    ("jets_DeepCSVbb", "DeepCSV Prob(bb)", 1, None, False),
    ("jets_DeepCSVc" , "DeepCSV Prob(c)" , 1, None, False),
    ("jets_DeepCSVbD", "DeepCSV b Disc", 1, None, False),
    ("jets_DeepCSVCvsB", "DeepCSV CvsB Disc", 1, None, False),
    ("jets_DeepCSVCvsL", "DeepCSV CvsL Disc", 1, None, False),
    ("jets_DeepJetCvsB", "DeepJet CvsB Disc", 1, None, False),
    ("jets_DeepJetCvsL", "DeepJet CvsL Disc", 1, None, False),
    ("jets_DeepFlavourb" ,   "DeepJet Prob(b)", 1, None, False),
    ("jets_DeepFlavourbb",   "DeepJet Prob(bb)", 1, None, False),
    ("jets_DeepFlavourlepb" ,"DeepJet Prob(lepb) ", 1, None, False),
    ("jets_DeepFlavourc" ,   "DeepJet Prob(c)", 1, None, False),
    ("jets_DeepFlavouruds" , "DeepJet Prob(uds)", 1, None, False),
    ("jets_DeepFlavourg" ,   "DeepJet Prob(g)", 1, None, False),
    ("jets_DeepFlavourbD",   "DeepJet b Disc", 1, None, False),
    ("jets_eta", "#eta(jet)", 10, None, False),
    ("jets_pt", "p_{T}(jet)", 10, None, False),
    ("lep_eta", "#eta(l)", 10, None, False),
    ("lep_pt", "p_{T}(l)", 10, None, False),
    ("nvtx", "# of reconstructed vertices", range(41), None, False),
    ("weight", "event weight", range(-3,4), None, False),
    ("rho", "#rho", range(40), None, False),
]

permutation_presel = [
	("mass_discriminant", "#lambda_{M}", 1, None, False),
	("Wmasshad", "M(W_{h})", 2, None, False),
	("tmasshad", "M(t_{h})", 2, None, False),
	("tpthad", "p_{T}(t_{h})", 2, None, False),
]


order = "mass_discriminant"

available_wps = [
	"notag",
	"DeepCSVLoose" ,
	"DeepCSVTight" ,
	"DeepCSVMedium",
	"DeepJetLoose" ,
	"DeepJetTight" ,
	"DeepJetMedium",
	"DeepCSVctagLoose",
	"DeepCSVctagMedium",
	"DeepCSVctagTight",
	"DeepJetctagLoose",
	"DeepJetctagMedium",
	"DeepJetctagTight",
]
working_points = [i for i in available_wps if fnmatch(i, args.wps)]

jet_categories = [
	("both_untagged" , 'notag'  ),
	("lead_tagged"	, 'leadtag'), 
	("sublead_tagged", 'subtag' ), 
	("both_tagged"	, 'ditag'  ), 
]
jet_types = ['leading', 'subleading']

def set_pdg_bins(histo):
	ids = {
		1 : 'd', #-1 : 'd',
		2 : 'u', #-2 : 'u',
		3 : 's', #-3 : 's',
		4 : 'c', #-4 : 'c',
		5 : 'b', #-5 : 'b',
		6 : 't', #-6 : 't',
		21: 'g',
		}
	
	for idx, name in ids.iteritems():
		bin_idx = histo.GetXaxis().FindBin(idx)
		histo.GetXaxis().SetBinLabel(bin_idx, name)
	return histo

additional_opts = {
	'btag_discriminant' : {
		'xrange' : [-11, -2]
		},
	'nu_chisq' : {'xrange' : [0, 20]},
}

if args.mceffs:
	json = {}
	for wp in working_points:
		json[wp] = {}
		for flav, idx in [('light', 1), ('charm', 2), ('beauty', 3)]: 
			json[wp][flav] = {}
			for ptc in ['leading', 'subleading']:
				json[wp][flav][ptc] = {}
				for sample in plotter.mc_samples:
					name, eff = plotter.get_efficiency(sample, wp, ptc, idx)
					json[wp][flav][ptc][name] = eff

	sample_view = views.SubdirectoryView(
		plotter.get_view('QCD*'),
		'nosys/mass_discriminant/ctagLoose'
		)
	## ntag = sample_view.Get('both_untagged/flav_map') 
	## ltag = sample_view.Get('lead_tagged/flav_map')	
	## stag = sample_view.Get('sublead_tagged/flav_map')
	## dtag = sample_view.Get('both_tagged/flav_map')	
	## total = ntag+ltag+stag+dtag
	## set_trace()
	with open('%s/info.json' % plotter.outputdir, 'w') as out:
		out.write(prettyjson.dumps(json))

if args.plots:
	#cut flow
	#flow = plotter.cut_flow()
	#plotter.save('cut_flow')
	##exit()

	#if args.eras == 'B':
	#    era_name = 'Run_B'
	#elif args.eras == 'CtoE':
	#    era_name = 'Runs_CtoE'
	#elif args.eras == 'EtoF':
	#    era_name = 'Runs_EtoF'
	#else:
	#    era_name = 'All_Runs'

	#set_trace()
	#plotter.set_subdir('%s/preselection' % era_name)
	plotter.set_subdir('preselection')
	for var, axis, rebin, x_range, leftside in preselection:
		plotter.make_preselection_plot(
			'nosys/preselection', var, sort=True,
			xaxis=axis, leftside=leftside, rebin=rebin, 
			show_ratio=True, ratio_range=0.5, sys_effs='pu' if var == 'nvtx' or var == 'rho' else None)		
		plotter.save(var)
		
	for var, axis, rebin, x_range, leftside in permutation_presel:
		plotter.make_preselection_plot(
			'nosys/permutations', var, sort=True,
			xaxis=axis, leftside=leftside, rebin=rebin, 
			show_ratio=True, ratio_range=0.5)
		plotter.save(var)

	##
	## Special plot
	##
	#plotter.draw_flavor_shapes('nosys/mass_discriminant/notag/both_untagged/', 'btag_sf', False, False) 
	#plotter.draw_cvsl_shapes('nosys/mass_discriminant/notag/both_untagged/', 'Wjets', 'hflav_CvsL', False, False) 
	#plotter.draw_cvsl_shapes('nosys/mass_discriminant/notag/both_untagged/', 'Wjets', 'hflav_CvsB', False, False) 
	#plotter.draw_cvsl_shapes('nosys/mass_discriminant/notag/both_untagged/', 'Wjets', 'hflav_DeepCSVCvsL', False, False) 
	#plotter.draw_cvsl_shapes('nosys/mass_discriminant/notag/both_untagged/', 'Wjets', 'hflav_DeepCSVCvsB', False, False) 
	#plotter.draw_cvsl_shapes('nosys/mass_discriminant/notag/both_untagged/', 'Wjets', 'hflav_DeepCSVbD', False, False) 
	###plotter.draw_cvsl_shapes('nosys/mass_discriminant/notag/both_untagged/', 'Wja', 'hflav_CvsL', False, False) 
	###plotter.draw_cvsl_shapes('nosys/mass_discriminant/notag/both_untagged/', 'Wjb', 'hflav_CvsL', False, False) 
	###plotter.draw_cvsl_shapes('nosys/preselection/', 'jets', 'hflav_CvsL', False, True) 


	def write_btag_muon_sf_hists():
		comparison_hists = []
        #
        #btag SFs
        #
		flavors = ['B','B_dw', 'B_up', 'C', 'C_dw', 'C_up', 'L', 'L_dw', 'L_up']
		for flavor in flavors:
			btag_hists = [i.Get('btag_sf_vs_pt_%s' % flavor) for i in plotter.mc_views(1, None, 'nosys/mass_discriminant/notag/both_untagged/')]
			#set_trace()
			sum_btag_hist = sum(btag_hists)
			sum_btag_hist.drawstyle='colz'
			sum_btag_hist.Draw()
			sum_btag_hist.SetXTitle('p_{T}(%s quarks) (GeV)' % flavor)
			sum_btag_hist.SetYTitle('b-tag SF')
			sum_btag_hist.name = 'btag_sf_vs_pt_%s' % flavor
			plotter.save('btag_sf_vs_pt_%s' % flavor)
			comparison_hists.append(sum_btag_hist)
			
		for flavor in flavors:
			btag_hists = [i.Get('btag_sf_%s' % flavor) for i in plotter.mc_views(1, None, 'nosys/mass_discriminant/notag/both_untagged/')]
			sum_btag_hist = sum(btag_hists)
			sum_btag_hist.Draw()
			sum_btag_hist.SetXTitle('b-tag SF (%s quarks)' % flavor)
			sum_btag_hist.SetYTitle('Events')
			sum_btag_hist.name = 'btag_sf_%s' % flavor
			plotter.save('btag_sf_%s' % flavor)
			comparison_hists.append(sum_btag_hist)

		systs = ['_dw', '_up', '']
		for sys in systs:
			btag_hists = [i.Get('btag_sf%s' % sys) for i in plotter.mc_views(1, None, 'nosys/mass_discriminant/notag/both_untagged/')]
			sum_btag_hist = sum(btag_hists)
			sum_btag_hist.Draw()
			sum_btag_hist.SetXTitle('b-tag SF (%s)' % sys)
			sum_btag_hist.SetYTitle('Events')
			sum_btag_hist.name = 'btag_sf%s' % sys
			plotter.save('btag_sf%s' % sys)
			comparison_hists.append(sum_btag_hist)


        #muon SFs
        #
		kinvars = {'pt' : 'p_{T}(#mu) (GeV)', 'eta' : '#eta(#mu)', 'pt_etaL0p9': 'p_{T}(#mu) (GeV)', 'pt_0p9eta1p2': 'p_{T}(#mu) (GeV)', 'pt_1p2eta2p1': 'p_{T}(#mu) (GeV)'}
		for kvar in kinvars.keys():
			musf_hists = [i.Get('muon_sf_vs_mu_%s' % kvar) for i in plotter.mc_views(1, None, 'nosys/mass_discriminant/notag/both_untagged/')]
			sum_musf_hist = sum(musf_hists)
			sum_musf_hist.drawstyle='colz'
			sum_musf_hist.Draw()
			sum_musf_hist.SetXTitle(kinvars[kvar])
			sum_musf_hist.SetYTitle('#mu SF')
			sum_musf_hist.name = 'muon_sf_vs_mu_%s' % kvar
			plotter.save('muon_sf_vs_mu_%s' % kvar)
			comparison_hists.append(sum_musf_hist)
			
        #set_trace()
		with io.root_open('%s/sf_comparison.root' % plotter.outputdir, 'w') as out:
			for i in comparison_hists:
				i.Write()



	#
	#Flavour maps
	#
	histograms = [i.Get('flav_map') for i in plotter.mc_views(1, None, 'nosys/mass_discriminant/notag/both_untagged/')]
	for i, j in enumerate(histograms):
		ROOT.gStyle.SetOptTitle(1)
		j.drawstyle='colz'
		j.Draw()
		plotter.save('flavour_map_%d' % i)
		ROOT.gStyle.SetOptTitle(0)
	with io.root_open('%s/falvour_maps.root' % plotter.outputdir, 'w') as out:
		for i in histograms:
			i.name = plotter.card_by_title[i.title]
			i.Write()

	#
	#C-Tag light/charm avg PTs
	#
	#charm
	stacks = plotter.make_stack(folder='nosys/mass_discriminant/notag/both_untagged/')
	charm = stacks.Get('Wjets_hflav_jpt_C')
	print "\n\nCharm Pts"
	print "sample\t\t\tIntegral\t\tMean pT\t\tpT RMS"
	for h in charm.hists:
		print h.title, '\t\t', h.Integral(), '\t', h.GetMean(), '\t', h.GetRMS()
	charm_signal = charm.hists[-1].Clone()
	
	charm = stacks.Get('Wjets_hflav_jpt_L')
	print "\n\nLight Pts"
	print "sample\t\t\tIntegral\t\tMean pT\t\tpT RMS"
	for h in charm.hists:
		print h.title, '\t\t', h.Integral(), '\t', h.GetMean(), '\t', h.GetRMS()
	light_signal = charm.hists[-1].Clone()
	with io.root_open('%s/jpts.root' % plotter.outputdir, 'w') as o:
		o.WriteTObject(charm_signal, 'cjets')
		o.WriteTObject(light_signal, 'ljets')

	for wpoint in working_points:
	  for cat_dir, cat_name in jet_categories:
		 if wpoint == 'notag' and cat_dir <> "both_untagged": continue
		 if wpoint != 'notag': continue
		 base = os.path.join(order, wpoint, cat_dir)
		 for jtype in jet_types:
			 folder = os.path.join(base, jtype)
			 plotter.set_subdir(os.path.join(order, wpoint, cat_name, jtype))
			 folder = os.path.join('nosys', folder)
			 logging.info(folder)
			 for var, rebin, axis, x_range in jet_variables:
				 plotter.plot_mc_vs_data(
					 folder, var, rebin, sort=True,
					 xaxis=axis % jtype, leftside=False, 
					 xrange=x_range, show_ratio=True, ratio_range=0.5)
				 if 'pflav' in var:
					 for h in plotter.keep:
						 if isinstance(h, plotting.HistStack):
							 set_pdg_bins(h)
				 plotter.save('_'.join((jtype,var)))
			 plotter.bake_pie(folder)

		 plotter.set_subdir(os.path.join(order, wpoint, cat_name))

		 #write_btag_muon_sf_hists()

		 for var, xaxis, yaxis, rebin in vars2D:
			 ROOT.gStyle.SetPalette(56)
			 path = os.path.join('nosys', base, var)
			 plotter.plot_shape('data', path, drawopt='colz', xaxis=xaxis, yaxis=yaxis)
			 plotter.keep[0].xaxis.SetTitleOffset(1.2)
			 plotter.keep[0].yaxis.SetTitleOffset(1.5)
			 plotter.save(var)
		 
		 if wpoint != 'notag': continue
		 for var, axis, rebin, x_range, leftside in variables:
			 if var == 'mass_discriminant': rebin = plotter.binning[wpoint][cat_name]
			 folder = os.path.join('nosys', base)
			 if var == 'evt_weight':
				plotter.plot_mc_vs_data(
				 folder, var, rebin,
				 xaxis=axis, leftside=leftside,
				 xrange=x_range, show_ratio=True, 
				 logy=True, #ratio_range=1.0)
				 ratio_range=0.5)
			 else:
				plotter.plot_mc_vs_data(
				 folder, var, rebin,
				 xaxis=axis, leftside=leftside,
				 xrange=x_range, show_ratio=True, 
				 #ratio_range=1.0)
				 ratio_range=0.5)
			# plotter.plot_mc_vs_data(
			#	 folder, var, rebin,
			#	 xaxis=axis, leftside=leftside,
			#	 xrange=x_range, show_ratio=True, 
			#	 ratio_range=1.0)
			#	 #ratio_range=0.5)
			 plotter.save(var)

#
# separated plots of systematics uncertainties 
#
if args.systematics:
	for wpoint in working_points:
		for cat_dir, cat_name in jet_categories:			
			if wpoint == 'notag' and cat_dir <> "both_untagged": continue
			base = os.path.join(order, wpoint, cat_dir)

			if args.pdfs:
				plotter.set_subdir(os.path.join(order, wpoint, cat_name, 'systematics', 'PDFs'))
				samples = [i for i in plotter.mc_samples if i.startswith('ttJets')]
				for sample in samples:
					plotter.pdf_unc_plots(sample, 'nosys/%s' % base, 'mass_discriminant', plotter.binning[wpoint][cat_name])
					plotter.save(sample)

			for systematic in systematics_to_check:
				plotter.set_subdir(os.path.join(order, wpoint, cat_name, 'systematics', systematic))
				samples = [re.compile(i) for i in plotter.systematics[systematic]['samples']]
				for cardname, grouping in plotter.card_names.iteritems():
					if cardname == 'data_obs': continue
					if not any(i.match(cardname) for i in samples): continue

					path	= 'nosys/%s/mass_discriminant' % base
					path_p = plotter.systematics[systematic]['+'](path)
					path_m = plotter.systematics[systematic]['-'](path) if '-' in plotter.systematics[systematic] else None
					#merge groups
					groupview = Plotter.rebin_view(
						views.StyleView(
							views.SumView(
								*[plotter.get_view(i) for i in grouping]
								 ),
							fillstyle = 'hollow',
							linestyle = 1,
							linewidth = 2,
							drawstyle = 'hist e',
							markersize = 0,
							),
						plotter.binning[wpoint][cat_name],
						)

					histo	= groupview.Get(path)
					histo_p = groupview.Get(path_p)
					histo_m = groupview.Get(path_m) if path_m else None

					to_overlay = [histo_p]
					histo.title = '%s central' % cardname
					histo.linecolor = 'black'
					histo_p.title = '%s+' % systematic
					histo_p.linecolor = 'red'
					histo_p.markercolor = 'red'
					if histo_m:
						histo_m.title = '%s-' %systematic
						histo_m.linecolor = 'blue'
						histo_m.markercolor = 'blue'
						to_overlay.append(histo_m)
					plotter.overlay_and_compare(
						to_overlay, histo, method='ratio', lower_y_range=0.55,
						legend_def=None, xtitle='mass discriminant', ytitle='events'
						)
					plotter.save('%s_%s' % (cardname, systematic))


if args.shapes:
	if args.noLightFit:
		plotter.systematics['CTAGL']['samples'].append('right_whad')
	plotter.set_subdir("")
	info = plotter.write_summary_table(working_points)
	shifts = {'lcfrac' : {}, 'scfrac' : {}}
	lfer = info['lcfrac'].std_dev
	sfer = info['scfrac'].std_dev	
	for name, shift in plotter.systematics.iteritems():
		if 'constants' in shift:			
			dw, up = shift['constants']
			su = plotter.write_summary_table(working_points, dw) if dw is not None else info
			sd = plotter.write_summary_table(working_points, up) if up is not None else info
			for attr in ['lcfrac', 'scfrac']:
				delta = max(
					abs(info[attr].nominal_value - su[attr].nominal_value),
					abs(info[attr].nominal_value - sd[attr].nominal_value)
					)
				info[attr] += ufloat(0, delta)
				## lead_cfrac
				## sub_cfrac
				#quadratic eqn factors
				## du = (su[attr].nominal_value - info[attr].nominal_value)/info[attr].nominal_value
				## dd = (sd[attr].nominal_value - info[attr].nominal_value)/info[attr].nominal_value
				## a = (du+dd)/2
				## b = (du-dd)/2
				## shifts[attr][name] = '%f*@0*@0%s%f*@0+1' % (a, '+' if b>=0 else '', b)
	
	## print 'lcfrac', info['lcfrac'], lfer
	## print 'scfrac', info['scfrac'], sfer
	cwd = os.getcwd()
	
	#
	# Common shapes, not enough stats to have one each category
	#	
	all_tag = [os.path.join('nosys', 'mass_discriminant', 'notag', 'both_untagged')]
	common_shapes = { #pointless to check effect of JES on samples with such low statitics
		'qcd$' : plotter.get_shape(all_tag, 'qcd', rebin=1) if 'qcd' in plotter.card_names else None,
		## 'qcd.+JESUp$' : plotter.get_shape(all_tag, 'qcd', rebin=mass_discriminant_binning, systematic='JES+'),
		## 'qcd.+JESDown$' : plotter.get_shape(all_tag, 'qcd', rebin=mass_discriminant_binning, systematic='JES-'),
		'vjets$' : plotter.get_shape(all_tag, 'vjets', rebin=1) if 'vjets' in plotter.card_names else None,
		## 'vjets.+JESUp$' : plotter.get_shape(all_tag, 'vjets', rebin=mass_discriminant_binning, systematic='JES+'),
		## 'vjets.+JESDown$' : plotter.get_shape(all_tag, 'vjets', rebin=mass_discriminant_binning, systematic='JES-'),
		}

	for wpoint in working_points:
		wpoint_dir = os.path.join(order, wpoint)
		plotter.set_subdir(wpoint_dir)
		## fname = os.path.join(shapedir, 'shapes.root')
		## shape_file = ROOT.TFile(fname, 'RECREATE')
		
		#plotter.write_mass_discriminant_shapes(
		#	shape_file.mkdir('notused'),
		#	os.path.join('all', order, 'notag', 'both_untagged'), 
		#	rebin=2
		#	)
		categories= ['Inc_nolead', 'Inc_nosub', 'Inc_leadtag', 'Inc_subtag'] \
			if args.inclusive else ['notag', 'leadtag', 'subtag', 'ditag']
		folders = [['both_untagged', 'sublead_tagged'], ['both_untagged', 'lead_tagged'], ['lead_tagged', 'both_tagged'], ['sublead_tagged', 'both_tagged']]\
			if args.inclusive else [['both_untagged'], ['lead_tagged'], ['sublead_tagged'], ['both_tagged']]
		if wpoint == 'notag':
			categories = ['notag']
			folders = [['both_untagged']]

		base = os.path.join('nosys', order, wpoint)
		for folders, category in zip(folders, categories):
			right_whad_yields = plotter.write_mass_discriminant_shapes(
				category,
				[os.path.join(base, i) for i in folders], 
				rebin=plotter.binning[wpoint][category]
				)

		#set_trace()
		#
		# yield formulas
		#
		yield_formulas = {
			'notag'	: '((1-{LCharmE})*(1-{sub_lightEff}*(1+x))*{lead_cfrac}+(1-{lead_lightEff}*(1+x))*(1-{SCharmE})*{sub_cfrac}+(1-{lead_lightEff}*(1+x))*(1-{sub_lightEff}*(1+x))*(1-{lead_cfrac}-{sub_cfrac}))',
			'leadtag' : '({LCharmE}*(1-{sub_lightEff}*(1+x))*{lead_cfrac}+{lead_lightEff}*(1+x)*(1-{SCharmE})*{sub_cfrac}+{lead_lightEff}*(1+x)*(1-{sub_lightEff}*(1+x))*(1-{lead_cfrac}-{sub_cfrac}))',
			'subtag'  : '((1-{LCharmE})*{sub_lightEff}*(1+x)*{lead_cfrac}+(1-{lead_lightEff}*(1+x))*{SCharmE}*{sub_cfrac}+(1-{lead_lightEff}*(1+x))*{sub_lightEff}*(1+x)*(1-{lead_cfrac}-{sub_cfrac}))',
			'ditag'	: '({LCharmE}*{sub_lightEff}*(1+x)*{lead_cfrac}+{lead_lightEff}*(1+x)*{SCharmE}*{sub_cfrac}+{lead_lightEff}*(1+x)*{sub_lightEff}*(1+x)*(1-{lead_cfrac}-{sub_cfrac}))',
			} #TODO add formulas for inclusive categories
		lightctag_nuisance_name = 'CTAGL'
		ctag_nuisance_name = 'CTAGC'
		category_constants = {
			'parshifts' : shifts, 
			'light_SF' : {
				'nuisance_name' : lightctag_nuisance_name, #have another systematic to take care of it
				}, 
			'charm_SF' : {'floating' : 0.5}, #float between 0.5/1.5
			'right_whad_yields' : right_whad_yields #stores right_whad yields for all categories
			}
		for category in categories:
			#
			# Remove CTAGL systematic from right_whad sample, get the value 
			#
			if args.noLightFit:
				value = None
				for idx in range(len(plotter.card.systematics[lightctag_nuisance_name].applies_)):
					cat_pat, sample_pat, value = plotter.card.systematics[lightctag_nuisance_name].applies_[idx]
					if cat_pat.match(category) and sample_pat.match('right_whad'):
						del plotter.card.systematics[lightctag_nuisance_name].applies_[idx]
						break
				yield_str = yield_formulas[category].format(
					LCharmE = info[wpoint]['lead']['ceff'].nominal_value,
					SCharmE = info[wpoint]['sub']['ceff'].nominal_value,
					lead_lightEff = info[wpoint]['lead']['leff'].nominal_value,
					sub_lightEff  = info[wpoint]['sub' ]['leff'].nominal_value,
					lead_cfrac = info['lcfrac'].nominal_value,
					sub_cfrac  = info['scfrac'].nominal_value,
					)
				default = eval(yield_str.replace('x', '0.'))
				formula_str = '{0} - {1}/{2}'.format(value, yield_str, default)
				formula = ROOT.TF1('blah', formula_str, -2, 2)
				#look for sign changes
				xvals = [i/100. for i in range(-200, 200)]
				vals  = [formula(i) for i in xvals]
				crossings = []
				for i in range(len(xvals)-1):
					if vals[i]*vals[i+1] < 0:
						crossings.append(xvals[i])
				crossings.sort(key=lambda x: abs(x))
				if crossings:
					constant = formula.GetX(0, crossings[0], crossings[0]+0.01)				
				else:
					constant = formula.GetX(0, -2, 2)
				category_constants['light_SF'][category] = constant
				#set_trace()
				if len(crossings) > 1:
					logging.warning(
						"multiple roots found for %s, %s storing the smallest as default and the other as alternate" % (wpoint, category)
						)
					constant = formula.GetX(0, crossings[1], crossings[1]+0.01)				
					category_constants['light_SF']['ALTERNATE %s' % category] = constant
			elif not args.noPOIpropagation:
				category_constants['light_SF'][category] = {}
				for cat_pat, sample_pat, value in plotter.card.systematics[lightctag_nuisance_name].applies_:
					if not cat_pat.match(category): continue
					sample = [i for i in plotter.card_names if sample_pat.match(i)][0]
					ang_par = (value-1)/category_constants['light_SF']['floating']
					category_constants['light_SF'][category][sample] = '{ang}*(@0-1)+1'.format(ang=ang_par)
				

			#
			# Remove all CTAGC systematics and store the effect on the constants json,
			# it will be fed to the model that will take care to propagate it proportional
			# to the POI (flat prior)
			#
			if not args.noPOIpropagation:
				category_constants['charm_SF'][category] = {}
				for cat_pat, sample_pat, value in plotter.card.systematics[ctag_nuisance_name].applies_:
					if not cat_pat.match(category): continue
					sample = [i for i in plotter.card_names if sample_pat.match(i)][0]
					ang_par = (value-1)/category_constants['charm_SF']['floating']
					category_constants['charm_SF'][category][sample] = '{ang}*(@0-1)+1'.format(ang=ang_par)
				
		if not args.noPOIpropagation:
			del plotter.card.systematics[ctag_nuisance_name]
		if not args.noLightFit and not args.noPOIpropagation: 
			del plotter.card.systematics[lightctag_nuisance_name]

		#
		# Replace V+Jets and QCD shapes with cumulative ones
		#
		for pattern, shape in common_shapes.iteritems():
			if shape is None: continue
			for category in categories:
				binned_shape = shape.Clone()
				binned_shape = urviews.RebinView.rebin(binned_shape, plotter.binning[wpoint][category])
				sys = ('Up' not in pattern and 'Down' not in pattern)
				plotter.card.replace_shape(category, pattern, binned_shape, sys)

		plotter.card.clamp_negative_bins('.*', '.*')
		plotter.set_subdir(wpoint_dir)
		if wpoint == 'notag': 
			plotter.save_card('datacard')
			continue
		#
		# Add Bin-by-bin uncertainties and plot their effect 
		#
		bbbs = plotter.add_systematics(args.noBBB)
		#set_trace()
		for category in categories:
			plotter.set_subdir('%s/%s/systematics/BBB' % (wpoint_dir, category))
			systematics = [i for i in bbbs if category in i]
			groups = {}
			for name in systematics:
				sample, bin = tuple(name.replace('%s_' % category, '').split('_bin_'))
				bin  = int(bin)
				if sample not in groups:
					groups[sample] = []
				groups[sample].append((name, bin))
			
			for sample, syss in groups.iteritems():
				central = plotter.card[category][sample].Clone()
				central.title = '%s central' % sample
				central.linecolor = 'black'
				central.fillstyle = 'hollow'
				central.drawstyle = 'hist'

				up = central.Clone()
				down = central.Clone()
				up.linecolor = 'red'
				down.linecolor = 'blue'
				for name, ibin in syss:
					up[ibin].value = plotter.card[category]['%s_%sUp' % (sample, name)][ibin].value
					down[ibin].value = plotter.card[category]['%s_%sDown' % (sample, name)][ibin].value

				plotter.overlay(
					[up, down, central], ytitle='events',
					legend_def=None, xtitle='mass discriminant', 
					)
				plotter.save('%s_BBB' % (sample))
				
		plotter.set_subdir(wpoint_dir)
		plotter.card.add_systematic('lead_cfrac', 'param', '', '',		  info['lcfrac'].nominal_value, info['lcfrac'].std_dev)
		plotter.card.add_systematic('sub_cfrac' , 'param', '', '',		  info['scfrac'].nominal_value, info['scfrac'].std_dev)
		plotter.card.add_systematic('mc_lead_charm_eff', 'param', '', '', info[wpoint]['lead']['ceff'].nominal_value, info[wpoint]['lead']['ceff'].std_dev)
		plotter.card.add_systematic('mc_lead_light_eff', 'param', '', '', info[wpoint]['lead']['leff'].nominal_value, info[wpoint]['lead']['leff'].std_dev)
		plotter.card.add_systematic('mc_sub_charm_eff' , 'param', '', '', info[wpoint]['sub']['ceff' ].nominal_value, info[wpoint]['sub']['ceff' ].std_dev)
		plotter.card.add_systematic('mc_sub_light_eff' , 'param', '', '', info[wpoint]['sub']['leff' ].nominal_value, info[wpoint]['sub']['leff' ].std_dev)
		plotter.card.add_systematic('signal_norm', 'param', '', '', info['total'].nominal_value, info['total'].std_dev)
		if args.noLightFit:
			plotter.card.add_comment("to be used without fitting the light SF -- NOLIGHTSFFIT")
		if args.noPOIpropagation:
			plotter.card.add_comment("this card was set to run without POI propagation to backgrounds -- NOPOIPROPAGATION")
		with open(os.path.join(plotter.outputdir, 'datacard.json'), 'w') as f:
			f.write(prettyjson.dumps(category_constants))

		plotter.save_card('datacard')
