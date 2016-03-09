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
                    help='use inclusive categories')
parser.add_argument('--pdfs', action='store_true')
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
      self.tt_to_use = 'ttJets'
      self.tt_shifted = {
         'mtop_up' : 'ttJets_mtopup',
         'mtop_down' : 'ttJets_mtopdown',
         'hadscale_up' : 'ttJets_scaleup', 
         'hadscale_down' : 'ttJets_scaledown',
         }
      jobid = os.environ['jobid']
      files = glob.glob('results/%s/ctag_eff/*.root' % jobid)
      logging.debug('files found %s' % files.__repr__())
      lumis = glob.glob('inputs/%s/*.lumi' % jobid)
      logging.debug('lumi files found %s' % lumis.__repr__())
      
      outdir= 'plots/%s/ctageff' % jobid
      lscaling = 0.5
      super(CTagPlotter, self).__init__(
         files, lumis, outdir, styles, None, lumi,
         lumi_scaling=lscaling
         )
      self.defaults = {
         'blurb' : [13, self.views['data']['intlumi']]
         }
      self.jobid = jobid

      self.views['ttJets_preselection'] = {
         'view' : self.create_tt_subsample(
            ['semilep_visible_right'],
            't#bar{t}',
            '#6666b3'
            )
         }
      self.views['ttJets_sig'] = {
         'view' : self.create_tt_subsample(
            ['semilep_visible_right', 'semilep_right_thad', 'semilep_right_whad'], 
            't#bar{t}, right W_{h}',
            '#6666b3'
            )
         }
      self.views['ttJets_bkg'] = {
         'view' : self.create_tt_subsample(
            ['semilep_right_tlep', 'semilep_wrong'], 
            't#bar{t}, wrong W_{h}',
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
         'single*',
         'ttJets_other',
         'ttJets_bkg',
         'ttJets_sig',
         ]

      self.card_names = {
         'qcd' : ['QCD*'],
         'vjets'      : ['[WZ]Jets'],
         'right_whad' : ['ttJets_sig'],
         'wrong_whad' : ['ttJets_bkg'],
         'nonsemi_tt' : ['ttJets_other'],
         'single_top' : ['single*'],
         'data_obs'   : ['data']
         }
      self.signal = 'right_whad'
      self.signal_yields = {}

      self.systematics = {
         'lumi' : {
            'type' : 'lnN',
            'samples' : ['.*'],
            'categories' : ['.*'],
            'value' : 1.05,
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
            'samples' : ['.*'],
            'categories' : ['.*'],
            'type' : 'shape',
            '+' : lambda x: x.replace('nosys', 'jer_up'),
            '-' : lambda x: x.replace('nosys', 'jer_down'),
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
            'value' : 1.00,            
            },
         'JES' : {
            'samples' : ['.*_whad', 'nonsemi_tt', 'single_top'],
            'categories' : ['.*'],
            'type' : 'shape',
            '+' : lambda x: x.replace('nosys', 'jes_up'),
            '-' : lambda x: x.replace('nosys', 'jes_down'),
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
            'samples' : ['.*'],
            'categories' : ['.*'],
            'type' : 'lnN',
            '+' : lambda x: x.replace('nosys', 'btagb_up'),
            '-' : lambda x: x.replace('nosys', 'btagb_down'),
            },
         'CTAGL' : {
            'samples' : ['wrong_whad', 'nonsemi_tt', 'single_top', 'qcd', 'vjets', 'right_whad'],
            'categories' : ['.*'],
            'type' : 'lnN',
            '+' : lambda x: x.replace('nosys', 'btagl_up'),
            '-' : lambda x: x.replace('nosys', 'btagl_down'),
            },
         'CTAGC' : {
            'samples' : ['wrong_whad', 'nonsemi_tt', 'single_top', 'qcd', 'vjets', 'right_whad'],
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
            'notag'   : [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20],
            },
         'csvLoose' : {
            'notag'   : [6, 8, 10, 12, 14, 16, 18, 20],
            'leadtag' : [6, 8, 10, 12, 14, 16, 18, 20],
            'subtag'  : [6, 8, 10, 12, 14, 16, 18, 20],
            'ditag'   : [6, 8, 10, 12, 14, 16, 18, 20],   
            },
         'csvMedium' : {
            'notag'   : [6, 8, 10, 12, 14, 16, 18, 20],
            'leadtag' : [6, 8, 10, 12, 14, 16, 20],
            'subtag'  : [6, 10, 14, 20],
            'ditag'   : [6, 10, 20],
            },
         'csvTight' : {
            'notag'   : [6, 8, 10, 12, 14, 16, 18, 20],
            'leadtag' : [6, 10, 14, 20],
            'subtag'  : [6, 10, 20],
            'ditag'   : [6, 20],   
            },         
         'ctagLoose' : {
            'notag'   : [6, 10, 20],
            'leadtag' : [6, 10, 14, 20],
            'subtag'  : [6, 8, 10, 12, 14, 16, 18, 20],
            'ditag'   : [6, 8, 10, 12, 14, 16, 18, 20],
            },
         'ctagMedium' : {
            'notag'   : [6, 8, 10, 12, 14, 16, 18, 20],
            'leadtag' : [6, 8, 10, 12, 14, 16, 18, 20],
            'subtag'  : [6, 8, 10, 12, 14, 16, 18, 20],
            'ditag'   : [6, 8, 10, 12, 14, 16, 18, 20],
            },
         'ctagTight' : {
            'notag'   : [6, 10, 14, 20],
            'leadtag' : [6, 8, 10, 12, 14, 16, 18, 20],
            'subtag'  : [6, 8, 10, 12, 14, 16, 18, 20],
            'ditag'   : [6, 10, 20],
            },
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
         #plotter.card.add_bbb_systematics('.*', 'qcd') 

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
            sbin.value += qbin.value
            sbin.error = quad.quad(sbin.error, qbin.error)
      samples.append(qcd_histo)
      self.keep.append(qcd_histo)
      samples.sort(key=lambda x: x[-2].value)
      stack = plotting.HistStack()
      self.keep.append(stack)
      for i in samples:         
         stack.Add(i)

      self.style_histo(stack)
      self.style_histo(histo, **histo.decorators)

      histo.Draw() #set the proper axis labels
      histo.yaxis.title = 'Events'
      data = self.get_view('data').Get('cut_flow')
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
         if name == self.signal:
            histo.Scale(1./integral)
            self.signal_yields[category_name] = integral
         elif name == 'data_obs' and data_is_fake:
            int_int = float(int(integral))
            histo.Scale(int_int/integral)
         
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
               category['%s_%sUp'   % (name, sys_name)] = hu
               category['%s_%sDown' % (name, sys_name)] = hd               
               plotter.card.add_systematic(sys_name, 'shape', category_name, name, 1.00)
            elif info['type'] == 'shape' or 'value' not in info:
               paths_up = [info['+'](path) for path in paths] if '+' in info else None
               paths_dw = [info['-'](path) for path in paths] if '-' in info else None
               hup = sum(view.Get(i) for i in paths_up) if paths_up else None
               hdw = sum(view.Get(i) for i in paths_dw) if paths_dw else None

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
                  value = 1.00+sign*delta

                  plotter.card.add_systematic(
                     sys_name, info['type'],
                     category_name, name, value
                  )

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
                  category['%s_%sUp'   % (name, sys_name)] = hup
                  category['%s_%sDown' % (name, sys_name)] = hdw               
               

   @staticmethod
   def compute_eff(jmap, jrank, qtype):
      total = ['both_tagged', 'both_untagged', 'lead_tagged', 'sublead_tagged']
      passing = ['both_tagged']
      if jrank == 'leading': passing.append('lead_tagged')
      else: passing.append('sublead_tagged')
      
      n_total = sum(jmap[i][jrank][qtype] for i in total)
      n_passing = sum(jmap[i][jrank][qtype] for i in passing)
      return n_passing / n_total, math.sqrt(n_passing * (1 - n_passing / n_total)) / n_total

   def write_summary_table(self, orders, wpoints):
      right_wpoints = [i for i in wpoints if i <> 'notag']
      ritghtW_view = self.get_view('ttJets_sig') 
      mc_weight = self.views[self.tt_to_use]['weight']

      info = {}
      pflav_path = 'nosys/{order}/{wpoint}/{jtag}/{jrank}/hadronflav' #abs_pflav_smart'
      for order in orders:
         all_lead_pflav = ritghtW_view.Get(
            pflav_path.format(
               order = order,
               wpoint = 'notag',
               jtag = 'both_untagged',
               jrank = 'leading'
               )
            )
         all_sub_pflav = ritghtW_view.Get(
            pflav_path.format(
               order = order,
               wpoint = 'notag',
               jtag = 'both_untagged',
               jrank = 'subleading'
               )
            )

         ##compute total number of events
         err  = array('d',[-1])
         nevt = all_lead_pflav.IntegralAndError(1, all_lead_pflav.GetNbinsX(), err)
         assert(abs(all_lead_pflav.Integral() - all_sub_pflav.Integral()) < 0.01)
         
         ##compute leading charm fraction 
         #all histograms are the same
         charm_bin = all_lead_pflav.FindBin(4)
         total_charms = 0
         ncharm  = all_lead_pflav.GetBinContent(
            charm_bin
            )
         total_charms += ncharm
         lead_charm_fraction = ncharm/nevt
         lead_charm_fraction_err = math.sqrt(
            (lead_charm_fraction * (1- lead_charm_fraction))*mc_weight/nevt
            )

         ##compute subleading charm fraction 
         ncharm  = all_sub_pflav.GetBinContent(
            charm_bin
            )
         total_charms += ncharm
         sub_charm_fraction = ncharm/nevt
         sub_charm_fraction_err = math.sqrt(
            (sub_charm_fraction * (1- sub_charm_fraction))*mc_weight/nevt
            )
         total_lights = 2*nevt - total_charms

         mc_effs = {}

         ##compute MC charm efficiency
         for wp in right_wpoints:
            ncharm_tagged = 0
            nlight_tagged = 0
            mc_effs[wp] = {}
            info = {}
            for jtag in ["lead_tagged", "sublead_tagged", "both_tagged", "both_untagged"]:
               info[jtag] = {}
               for jrank in ['leading', 'subleading']:
                  histo = ritghtW_view.Get(
                     pflav_path.format(
                        order  = order,
                        wpoint = wp,
                        jtag   = jtag,
                        jrank  = jrank,
                        )
                     )
                  _, yields = self.make_flavor_table(
                     histo, 
                     fname='', 
                     to_json=False, 
                     to_txt=False,
                     hadron_flavour=True
                     )
                  info[jtag][jrank] = yields
                  err = array('d',[0])
                  mc_effs[wp][jtag] = histo.IntegralAndError(
                     1, 
                     histo.GetNbinsX(), 
                     err
                     )
                  mc_effs[wp]['%s_err' %jtag] = err[0]
            
            lead_ceff, lead_ceff_err = CTagPlotter.compute_eff(info, 'leading', 'charm')
            sub_ceff , sub_ceff_err = CTagPlotter.compute_eff(info, 'subleading', 'charm')
            lead_leff, lead_leff_err = CTagPlotter.compute_eff(info, 'leading', 'light')
            sub_leff , sub_leff_err = CTagPlotter.compute_eff(info, 'subleading', 'light')
            
            mc_effs[wp].update({
               'lead_charmEff' : lead_ceff,
               'sub_charmEff'  : sub_ceff ,
               'lead_lightEff' : lead_leff,
               'sub_lightEff'  : sub_leff ,

               'lead_charmEff_err' : lead_ceff_err,
               'sub_charmEff_err'  : sub_ceff_err,
               'lead_lightEff_err' : lead_leff_err,
               'sub_lightEff_err'  : sub_leff_err,
               })
         info[order] = {
            'nevts' : nevt,
            'nevts_err' : err[0],
            'leadCFrac' : lead_charm_fraction,
            'subCFrac' : sub_charm_fraction,
            'leadCFrac_err' : lead_charm_fraction_err,
            'subCFrac_err' : sub_charm_fraction_err,
            'working_points' : mc_effs
           } 

      with open(os.path.join(self.outputdir, 'summary.json'), 'w') as f:
         f.write(prettyjson.dumps(info))
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
      var = '%s/abs_pflav_smart' % path
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

      ## right_cmb = self.get_view('ttJets_allRight').Get(var)
      ## self.make_flavor_table(right_cmb, 'flavors_rightcmb.raw_txt')
      ## 
      ## wright_view = views.SumView(
      ##    *[self.get_view(i) for i in ['ttJets_rightHad', 'ttJets_rightWHad']]
      ##     )
      ## wright = wright_view.Get(var) + right_cmb
      ## self.make_flavor_table(wright, 'flavors_rightw.raw_txt', to_json=True)

   def make_preselection_plot(self, *args, **kwargs):
      mc_default = self.mc_samples
      self.mc_samples = ['QCD*', '[WZ]Jets', 'single*', 'ttJets_preselection']
      self.plot_mc_vs_data(*args, **kwargs)
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
      nnpdf = [1] + range(10, 110)
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

plotter = CTagPlotter(args.lumi)

jet_variables = [
   ('energy', 5, '%s jet energy (GeV)', None),
   ('pt' ,  10, '%s jet p_{T} (GeV)', None),
   ('eta',   1, '%s jet #eta', None),
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
   'HScale',
   'pu'
]

vars2D = [
   ('Whad_jet_pts', 'p_{T}(lead W jet) (GeV)', 'p_{T}(sub W jet) (GeV)', (1,1)),
]

variables = [
  ("njets"    , "# of selected jets", range(13), None, False),
  ("lep_pt"   , "p_{T}(l) (GeV)", 20, None, False),
  ("Whad_mass", "m_{W}(had) (GeV)", 10, None, False),
  ("thad_mass", "m_{t}(had) (GeV)", 10, None, False),
	("mass_discriminant", "#lambda_{M}", 1, None, False), #[5, 20]),
  ("Wjets_CvsL", "CvsL (W Jet)", 1, None, False),
  ("Wjets_CvsB", "CvsB (W Jet)", 1, None, True),
  ("Bjets_CvsB", "CvsB (B Jet)", 1, None, False),
  ("Bjets_CvsL", "CvsL (B Jet)", 1, None, True),
  #("Wlep_mass", "m_{W}(lep) (GeV)", 10, None, False),
  #("Whad_DR"  , "#DeltaR(jj) W_{had} (GeV)", 1, [0,7]),
  #("Whad_pt"  , "p_{T}(W_{had}) (GeV)", 10, None, False),
  #("Whad_leading_DR", "#DeltaR(W_{had}, leading jet)"   , 1, [0,7]),
  #("Whad_sublead_DR", "#DeltaR(W_{had}, subleading jet)", 1, [0,7]),
  #("nu_chisq"         , "nu_chisq"         , 1, [0, 20]),
	#("nu_discriminant"	, "nu_discriminant"	 , 1, None, False),
	#("btag_discriminant", "btag_discriminant", 1, [-11, -2]),
  #("nbjets"   , "# of bjets", 1, [0, 12]),
  #("lep_b_pt" , "p_{T}(b) (GeV)", 10, None, False),
  #("had_b_pt" , "p_{T}(b) (GeV)", 10, None, False),
]

preselection = [
  ("njets"    , "# of selected jets", range(13), None, False),
  ("jets_CSV", "jet CSV",   2, None, False),
  ("jets_CvsB", "jet CvsB", 1, None, True),
  ("jets_CvsL", "jet CvsL", 1, None, False),
  ("jets_eta", "#eta(jet)", 10, None, False),
  ("jets_pt", "p_{T}(jet)", 10, None, False),
  ("lep_eta", "#eta(l)", 10, None, False),
  ("lep_pt", "p_{T}(l)", 10, None, False),
  ("nvtx", "# of reconstructed vertices", range(41), None, False),
  ("rho", "#rho", range(40), None, False),
]

permutation_presel = [
   ("mass_discriminant", "#lambda_{M}", 5, None, True),
   ("Wmasshad", "M(W_{h})", 2, None, False),
   ("tmasshad", "M(t_{h})", 2, None, False),
]


order = "mass_discriminant"

available_wps = [
   "notag",
   "csvTight",
   "csvMedium",
   "csvLoose",
   "ctagLoose",
   "ctagMedium",
   "ctagTight",
]
working_points = [i for i in available_wps if fnmatch(i, args.wps)]

jet_categories = [
   ("both_untagged" , 'notag'  ),
   ("lead_tagged"   , 'leadtag'), 
   ("sublead_tagged", 'subtag' ), 
   ("both_tagged"   , 'ditag'  ), 
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

if args.plots:
   #cut flow
   plotter.cut_flow()
   plotter.save('cut_flow')

   plotter.set_subdir('preselection')
   ##for var, axis, rebin, x_range, leftside in preselection:
   ##   plotter.make_preselection_plot(
   ##      'nosys/preselection', var, sort=True,
   ##      xaxis=axis, leftside=leftside, rebin=rebin, 
   ##      show_ratio=True, ratio_range=0.5)
   ##   plotter.save(var)
   ##   
   ##for var, axis, rebin, x_range, leftside in permutation_presel:
   ##   plotter.make_preselection_plot(
   ##      'nosys/permutations', var, sort=True,
   ##      xaxis=axis, leftside=leftside, rebin=rebin, 
   ##      show_ratio=True, ratio_range=0.5)
   ##   plotter.save(var)

   #
   # Special plot
   #
   for disc, leftleg in [('CvsL', False), ('CvsB', True)]:
      dirname = 'nosys/mass_discriminant/notag/both_untagged/'
      mc_views = plotter.mc_views(folder=dirname)   
      b_hists = [i.Get('Wjets_%s_B'  % disc) for i in mc_views]
      c_hists = [i.Get('Wjets_%s_C'  % disc) for i in mc_views]
      s_hists = [i.Get('Wjets_%s_S'  % disc) for i in mc_views]
      l_hists = [i.Get('Wjets_%s_UD' % disc) for i in mc_views]
      
      hb = sum(b_hists)
      hl = sum(l_hists)
      hc = sum(c_hists)
      hs = sum(s_hists)
      
      hb.decorate(**b_hists[0].decorators) 
      hb.title = 'b-jets'
      hl.decorate(**b_hists[1].decorators)
      hl.title = 'l-jets'
      hc.decorate(**b_hists[2].decorators)
      hc.title = 'c-jets'
      hs.decorate(**b_hists[-2].decorators)
      hs.title = 's-jets'
      stack = plotter.create_stack(hb, hl, hs, hc, sort=False)
      data = plotter.get_view('data').Get(dirname+'Wjets_%s' % disc)
      plotter.overlay([stack, data])
      plotter.add_legend([stack, data], leftleg, 5)
      plotter.save('Wjets_%s_flavour' % disc)   
      
      hb.Scale(1/hb.Integral())
      hl.Scale(1/hl.Integral())
      hc.Scale(1/hc.Integral())
      hs.Scale(1/hs.Integral())
      plotter.overlay([hb,hl,hc,hs], fillstyle='hollow', linewidth=3)
      plotter.add_legend([hb,hl,hc,hs], leftleg, 4)
      plotter.save('Wjets_%s_norm' % disc)   

   for wpoint in working_points:
     for cat_dir, cat_name in jet_categories:
       if wpoint == 'notag' and cat_dir <> "both_untagged": continue
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
       for var, xaxis, yaxis, rebin in vars2D:
          ROOT.gStyle.SetPalette(56)
          path = os.path.join('nosys', base, var)
          plotter.plot('data', path, drawopt='colz', xaxis=xaxis, yaxis=yaxis)
          plotter.keep[0].xaxis.SetTitleOffset(1.2)
          plotter.keep[0].yaxis.SetTitleOffset(1.5)
          plotter.save(var)

       for var, axis, rebin, x_range, leftside in variables:
          if var == 'mass_discriminant': rebin = plotter.binning[wpoint][cat_name]
          folder = os.path.join('nosys', base)
          plotter.plot_mc_vs_data(
             folder, var, rebin,
             xaxis=axis, leftside=leftside,
             xrange=x_range, show_ratio=True, 
             ratio_range=0.5)
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

##         for systematic in systematics_to_check:
##            plotter.set_subdir(os.path.join(order, wpoint, cat_name, 'systematics', systematic))
##            samples = [re.compile(i) for i in plotter.systematics[systematic]['samples']]
##            for cardname, grouping in plotter.card_names.iteritems():
##               if cardname == 'data_obs': continue
##               if not any(i.match(cardname) for i in samples): continue
##
##               path   = 'nosys/%s/mass_discriminant' % base
##               path_p = plotter.systematics[systematic]['+'](path)
##               path_m = plotter.systematics[systematic]['-'](path) if '-' in plotter.systematics[systematic] else None
##               #merge groups
##               groupview = Plotter.rebin_view(
##                  views.StyleView(
##                     views.SumView(
##                        *[plotter.get_view(i) for i in grouping]
##                         ),
##                     fillstyle = 'hollow',
##                     linestyle = 1,
##                     linewidth = 2
##                     ),
##                  plotter.binning[wpoint][cat_name],
##                  )
##
##               histo   = groupview.Get(path)
##               histo_p = groupview.Get(path_p)
##               histo_m = groupview.Get(path_m) if path_m else None
##
##               to_overlay = [histo_p]
##               histo.title = '%s central' % cardname
##               histo.linecolor = 'black'
##               histo_p.title = '%s+' % systematic
##               histo_p.linecolor = 'red'
##               histo_p.markercolor = 'red'
##               if histo_m:
##                  histo_m.title = '%s-' %systematic
##                  histo_m.linecolor = 'blue'
##                  histo_m.markercolor = 'blue'
##                  to_overlay.append(histo_m)
##               plotter.overlay_and_compare(
##                  to_overlay, histo, method='ratio', lower_y_range=0.55,
##                  legend_def=None, xtitle='mass discriminant', ytitle='events'
##                  )
##               plotter.save('%s_%s' % (cardname, systematic))


## plotter.set_subdir("discriminants")
## ## plotter.plot_mc_vs_data(
## ##    'all/discriminators', order, 1, sort=True,
## ##    xaxis=order.replace('_', ' '), leftside=False,
## ##    **additional_opts.get(order,{}))
## ## plotter.save(order, pdf=False)
## plotter.plot_mc_shapes(
##    'nosys/all/discriminators', order, 1,
##    xaxis=order.replace('_', ' '), leftside=False,
##    **additional_opts.get(order,{}))
## plotter.save('shape_%s' % order, pdf=False)
## 
## plotter.plot_mc_shapes(
##    'nosys/all/discriminators', order, 1, normalize=True,
##    xaxis=order.replace('_', ' '), leftside=False,
##    **additional_opts.get(order,{}))
## plotter.save('normshape_%s' % order, pdf=False)

if args.shapes:
   plotter.set_subdir("")
   info = plotter.write_summary_table([order], working_points)
   cwd = os.getcwd()

   #
   # Common shapes, not enough stats to have one each category
   #   
   all_tag = [os.path.join('nosys', 'mass_discriminant', 'notag', 'both_untagged')]
   common_shapes = { #pointless to check effect of JES on samples with such low statitics
      'qcd$' : plotter.get_shape(all_tag, 'qcd', rebin=1),
      ## 'qcd.+JESUp$' : plotter.get_shape(all_tag, 'qcd', rebin=mass_discriminant_binning, systematic='JES+'),
      ## 'qcd.+JESDown$' : plotter.get_shape(all_tag, 'qcd', rebin=mass_discriminant_binning, systematic='JES-'),
      'vjets$' : plotter.get_shape(all_tag, 'vjets', rebin=1),
      ## 'vjets.+JESUp$' : plotter.get_shape(all_tag, 'vjets', rebin=mass_discriminant_binning, systematic='JES+'),
      ## 'vjets.+JESDown$' : plotter.get_shape(all_tag, 'vjets', rebin=mass_discriminant_binning, systematic='JES-'),
      }

   for wpoint in working_points:
      if wpoint == "notag": continue
      wpoint_dir = os.path.join(order, wpoint)
      plotter.set_subdir(wpoint_dir)
      ## fname = os.path.join(shapedir, 'shapes.root')
      ## shape_file = ROOT.TFile(fname, 'RECREATE')
      
      #plotter.write_mass_discriminant_shapes(
      #   shape_file.mkdir('notused'),
      #   os.path.join('all', order, 'notag', 'both_untagged'), 
      #   rebin=2
      #   )
      categories= ['Inc_nolead', 'Inc_nosub', 'Inc_leadtag', 'Inc_subtag'] \
         if args.inclusive else ['notag', 'leadtag', 'subtag', 'ditag']
      folders = [['both_untagged', 'sublead_tagged'], ['both_untagged', 'lead_tagged'], ['lead_tagged', 'both_tagged'], ['sublead_tagged', 'both_tagged']]\
         if args.inclusive else [['both_untagged'], ['lead_tagged'], ['sublead_tagged'], ['both_tagged']]
      base = os.path.join('nosys', order, wpoint)
      for folders, category in zip(folders, categories):
         plotter.write_mass_discriminant_shapes(
            category,
            [os.path.join(base, i) for i in folders], 
            rebin=plotter.binning[wpoint][category]
            )

      #
      # yield formulas
      #
      yield_formulas = {
         'notag'   : '((1-{LCharmE})*(1-{sub_lightEff}*(1+x))*{lead_cfrac}+(1-{lead_lightEff}*(1+x))*(1-{SCharmE})*{sub_cfrac}+(1-{lead_lightEff}*(1+x))*(1-{sub_lightEff}*(1+x))*(1-{lead_cfrac}-{sub_cfrac}))',
         'leadtag' : '({LCharmE}*(1-{sub_lightEff}*(1+x))*{lead_cfrac}+{lead_lightEff}*(1+x)*(1-{SCharmE})*{sub_cfrac}+{lead_lightEff}*(1+x)*(1-{sub_lightEff}*(1+x))*(1-{lead_cfrac}-{sub_cfrac}))',
         'subtag'  : '((1-{LCharmE})*{sub_lightEff}*(1+x)*{lead_cfrac}+(1-{lead_lightEff}*(1+x))*{SCharmE}*{sub_cfrac}+(1-{lead_lightEff}*(1+x))*{sub_lightEff}*(1+x)*(1-{lead_cfrac}-{sub_cfrac}))',
         'ditag'   : '({LCharmE}*{sub_lightEff}*(1+x)*{lead_cfrac}+{lead_lightEff}*(1+x)*{SCharmE}*{sub_cfrac}+{lead_lightEff}*(1+x)*{sub_lightEff}*(1+x)*(1-{lead_cfrac}-{sub_cfrac}))',
         } #TODO add formulas for inclusive categories
      lightctag_nuisance_name = 'CTAGL'
      ctag_nuisance_name = 'CTAGC'
      category_constants = {'light_SF' : {'nuisance_name' : lightctag_nuisance_name, 'floating' : 0.5}, 'charm_SF' : {'floating' : 0.5}} #float between 0.5/1.5
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
               LCharmE = info[order]['working_points'][wpoint]['lead_charmEff'],
               SCharmE = info[order]['working_points'][wpoint]['sub_charmEff'],
               lead_lightEff = info[order]['working_points'][wpoint]['lead_lightEff'],
               sub_lightEff = info[order]['working_points'][wpoint]['sub_lightEff'],
               lead_cfrac = info[order]['leadCFrac'],
               sub_cfrac = info[order]['subCFrac'],
               )
            default = eval(yield_str.replace('x', '0.'))
            formula_str = '{0} - {1}/{2}'.format(value, yield_str, default)
            formula = ROOT.TF1('blah', formula_str, -2, 2)
            constant = formula.GetX(0, -2, 2)            
            category_constants['light_SF'][category] = constant
         else:
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
         category_constants['charm_SF'][category] = {}
         for cat_pat, sample_pat, value in plotter.card.systematics[ctag_nuisance_name].applies_:
            if not cat_pat.match(category): continue
            sample = [i for i in plotter.card_names if sample_pat.match(i)][0]
            ang_par = (value-1)/category_constants['charm_SF']['floating']
            category_constants['charm_SF'][category][sample] = '{ang}*(@0-1)+1'.format(ang=ang_par)
            
      del plotter.card.systematics[ctag_nuisance_name]
      if not args.noLightFit: del plotter.card.systematics[lightctag_nuisance_name]

      #
      # Replace V+Jets and QCD shapes with cumulative ones
      #
      for pattern, shape in common_shapes.iteritems():
         for category in categories:
            binned_shape = shape.Clone()
            binned_shape = urviews.RebinView.rebin(binned_shape, plotter.binning[wpoint][category])
            sys = ('Up' not in pattern and 'Down' not in pattern)
            plotter.card.replace_shape(category, pattern, binned_shape, sys)


      plotter.card.clamp_negative_bins('.*', '.*')
      #
      # Add Bin-by-bin uncertainties and plot their effect 
      #
      bbbs = plotter.add_systematics(args.noBBB)
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

      plotter.card.add_systematic('lead_cfrac', 'param', '', '', info[order]['leadCFrac'], info[order]['leadCFrac_err'])
      plotter.card.add_systematic('sub_cfrac' , 'param', '', '', info[order]['subCFrac'] , info[order]['subCFrac_err'] )
      plotter.card.add_systematic('mc_lead_charm_eff', 'param', '', '', info[order]['working_points'][wpoint]['lead_charmEff'], info[order]['working_points'][wpoint]['lead_charmEff_err'])
      plotter.card.add_systematic('mc_lead_light_eff', 'param', '', '', info[order]['working_points'][wpoint]['lead_lightEff'], info[order]['working_points'][wpoint]['lead_lightEff_err'])
      plotter.card.add_systematic('mc_sub_charm_eff' , 'param', '', '', info[order]['working_points'][wpoint]['sub_charmEff'] , info[order]['working_points'][wpoint]['sub_charmEff_err'] )
      plotter.card.add_systematic('mc_sub_light_eff' , 'param', '', '', info[order]['working_points'][wpoint]['sub_lightEff'] , info[order]['working_points'][wpoint]['sub_lightEff_err'] )
      plotter.card.add_systematic(
         'signal_norm', 'param', '', '', 
         info[order]['nevts'], 
         info[order]['nevts_err']
         )
      if args.noLightFit:
         plotter.card.add_comment("to be used without fitting the light SF -- NOLIGHTSFFIT")
      with open(os.path.join(plotter.outputdir, 'datacard.json'), 'w') as f:
         f.write(prettyjson.dumps(category_constants))

      plotter.save_card('datacard')

