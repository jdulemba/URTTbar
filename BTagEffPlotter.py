from pdb import set_trace
from URAnalysis.PlotTools.Plotter import Plotter
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
import re

parser = ArgumentParser()
parser.add_argument('--noplots', dest='noplots', action='store_true',
                    help='skip plot making')
parser.add_argument('--noshapes', dest='noshapes', action='store_true',
                    help='skip shape making')
parser.add_argument('--noLightFit', dest='noLightFit', action='store_true',
                    help='set up fitting without light fitting')
parser.add_argument('--noBBB', action='store_true',
                    help='do not run bin-by-bin uncertainties')
parser.add_argument('--inclusive', action='store_true',
                    help='use inclusive categories')
parser.add_argument('--lumi', type=float, default=-1.,
                    help='use inclusive categories')
args = parser.parse_args()

def syscheck(cmd):
   out = os.system(cmd)
   if out == 0:
      return 0
   else:
      raise RuntimeError("command %s failed executing" % cmd)

class BTagPlotter(Plotter):
   def __init__(self, lumi=None):
      lumi = lumi if lumi > 0 else None
      self.tt_to_use = 'ttJets_madgraph'
      jobid = os.environ['jobid']
      files = glob.glob('results/%s/btag_efficiency/*.root' % jobid)
      logging.debug('files found %s' % files.__repr__())
      lumis = glob.glob('inputs/%s/*.lumi' % jobid)
      logging.debug('lumi files found %s' % lumis.__repr__())
      
      outdir= 'plots/%s/btageff' % jobid
      super(BTagPlotter, self).__init__(
         files, lumis, outdir, styles, None, lumi
         )
      self.jobid = jobid

      self.views['ttJets_allRight'] = {
         'view' : self.create_tt_subsample(
            'semilep_visible_right', 
            'tt, right cmb',
            '#6666b3'
            )
         }
      self.views['ttJets_rightHad'] = {
         'view' : self.create_tt_subsample(
            'semilep_right_thad', 
            'tt, right t_{h}',
            '#aaaad5'
            ),
         }
      self.views['ttJets_rightWHad'] = {
         'view' : self.create_tt_subsample(
            'semilep_right_whad', 
            'tt, right W_{h}',
            '#cccce6'
            )
         }
      self.views['ttJets_rightWLep'] = {
         'view' : self.create_tt_subsample(
            'semilep_right_tlep', 
            'tt, right t_{l}',
            '#ab5555'
            )
         }
      self.views['ttJets_semiWrong'] = {
         'view' : self.create_tt_subsample(
            'semilep_wrong', 
            'tt, wrong cmb',
            '#d5aaaa'
            )
         }
      self.views['ttJets_other'] = {
         'view' : self.create_tt_subsample(
            'other', 
            'Other tt decay',
            '#668db3',
            )
         }

      self.mc_samples = [
         'QCD*',
         '[WZ]Jets',
         'single*',
         'ttJets_allRight',
         'ttJets_rightHad',
         'ttJets_rightWHad',
         'ttJets_rightWLep',
         'ttJets_semiWrong',
         'ttJets_other'
         ]

      self.card_names = {
         'qcd' : ['QCD*'],
         'right_whad' : ['ttJets_allRight', 'ttJets_rightHad', 'ttJets_rightWHad'],
         'wrong_whad' : ['ttJets_rightWLep', 'ttJets_semiWrong'],
         'nonsemi_tt' : ['ttJets_other'],
         'single_top' : ['single*'],
         }
      self.signal = 'right_whad'

      self.systematics = {
         'lumi' : {
            'type' : 'lnN',
            'samples' : ['.*'],
            'categories' : ['.*'],
            'value' : 1.05,
            },
         'JES' : {
            'samples' : ['.*'],
            'categories' : ['.*'],
            'type' : 'shape',
            '+' : lambda x: x.replace('nosys', 'jes_up'),
            '-' : lambda x: x.replace('nosys', 'jes_down'),
            'value' : 1.00,
            }
         }
      self.card = None
      self.binning = {}


   def create_tt_subsample(self, subdir, title, color='#9999CC'):
      return views.StyleView(
         views.TitleView(
            views.SubdirectoryView(
               self.views[self.tt_to_use]['view'],
               subdir
               ),
            title
            ),
         fillcolor = color,
         linecolor = color
         )

   def add_systematics(self, nobbb=False):
      if not plotter.card:
         raise ValueError('The card is not defined!')
      for sys_name, info in self.systematics.iteritems():
         for category in info['categories']:
            for sample in info['samples']:
               plotter.card.add_systematic(
                  sys_name, 
                  info['type'], 
                  category, 
                  sample, 
                  info['value']
                  )
      if not nobbb:
         plotter.card.add_bbb_systematics('.*', '.*')

   def save_card(self, name):
      if not self.card:
         raise RuntimeError('There is no card to save!')
      self.card.save(name, self.outputdir)
      self.card = None

   def write_mass_discriminant_shapes(self, category_name, folders, rebin=1, 
                                      preprocess=None):
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
         card_views['data_obs'] = views.SumView(*card_views.values())
         data_is_fake = True

      for name, view in card_views.iteritems():
         histo = sum(view.Get(path) for path in paths)
         integral = histo.Integral()
         if name == self.signal:
            histo.Scale(1./integral)
         elif name == 'data_obs' and data_is_fake:
            int_int = float(int(integral))
            histo.Scale(int_int/integral)
         
         category[name] = histo
         for sys_name, info in self.systematics.iteritems():
            if not any(re.match(i, category_name) for i in info['categories']): continue
            if not any(re.match(i, name) for i in info['samples']): continue
            shift = 1.0
            if info['type'] == 'shape':
               paths_up = [info['+'](path) for path in paths] 
               paths_dw = [info['-'](path) for path in paths]
               hup = sum(view.Get(i) for i in paths_up)
               hdw = sum(view.Get(i) for i in paths_dw)
               if name == 'right_whad':
                  hup.Scale(1./integral)
                  hdw.Scale(1./integral)
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
      ritghtW_view = views.SumView(
         self.get_view('ttJets_allRight'),
         self.get_view('ttJets_rightHad'),
         self.get_view('ttJets_rightWHad')
         )
      mc_weight = self.views[self.tt_to_use]['weight']

      info = {}
      pflav_path = 'nosys/{order}/{wpoint}/{jtag}/{jrank}/abs_pflav_smart'
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
         #set_trace()
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
                     to_txt=False
                     )
                  info[jtag][jrank] = yields
                  err = array('d',[0])
                  mc_effs[wp][jtag] = histo.IntegralAndError(
                     1, 
                     histo.GetNbinsX(), 
                     err
                     )
                  mc_effs[wp]['%s_err' %jtag] = err[0]
            
            lead_ceff, lead_ceff_err = BTagPlotter.compute_eff(info, 'leading', 'charm')
            sub_ceff , sub_ceff_err = BTagPlotter.compute_eff(info, 'subleading', 'charm')
            lead_leff, lead_leff_err = BTagPlotter.compute_eff(info, 'leading', 'light')
            sub_leff , sub_leff_err = BTagPlotter.compute_eff(info, 'subleading', 'light')
            
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
         #set_trace()
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

   def make_flavor_table(self, histo, fname='', to_json=False, to_txt=True):
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
      format = '%10s %10s +/- %10s %10s\n'
      yields = {
         'light' : 0, 'charm' : 0, 'other' : 0, 'strange' : 0,
         'light_err' : 0, 'charm_err' : 0, 'other_err' : 0, 'strange_err' : 0,
                }
      quark_yields = {}
      quark_errors = {}
      yield_names = {
         'light' : set([1,2,3]), 
         'charm' : set([4]), 
         'other' : set([0, 5, 6, 21]),
         'strange' : set([3]),
         }

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
      total = sum(integrals)
      ratios = [i/total for i in integrals]
      names = [i.GetTitle() for i in mc_histos]
      colors = [i.GetFillColor('root') for i in mc_histos]
      #set_trace()
      format = '%20s %10s %10s\n'
      with open(os.path.join(self.outputdir, 'yields.raw_txt'), 'w') as f:
         f.write('Total events: %.0f\n' % total)
         header = format % ('process', 'yield', 'relative')
         f.write(header)
         f.write('-'*len(header)+'\n')
         for name, entries, ratio in zip(names, integrals, ratios):
            f.write(format % (name, '%.0f' % entries, '%.0f%%' % (ratio*100)))

      hsum = sum(mc_histos)
      self.make_flavor_table(hsum, 'flavors.raw_txt', to_json=True)

      right_cmb = self.get_view('ttJets_allRight').Get(var)
      self.make_flavor_table(right_cmb, 'flavors_rightcmb.raw_txt')

      wright_view = views.SumView(
         *[self.get_view(i) for i in ['ttJets_rightHad', 'ttJets_rightWHad']]
          )
      wright = wright_view.Get(var) + right_cmb
      self.make_flavor_table(wright, 'flavors_rightw.raw_txt', to_json=True)
   
plotter = BTagPlotter(args.lumi)

jet_variables = [
   ('pt' ,  10, '%s jet p_{T} (GeV)', None),
   ('eta',   1, '%s jet #eta', None),
   ('phi',  10, '%s jet #varphi', None),
   ('abs_pflav_smart', 1, '%s jet parton flavour', None),#[-7, 22]),
   ("ncharged", 1, "%s jet charged constituents", None),
   ("nneutral", 1, "%s jet neutral constituents", None),
   ("ntotal"  , 1, "%s jet total constituents", None),
]

variables = [
  ("njets"    , "# of selected jets", 1, [0, 12]),
  ("nbjets"   , "# of bjets", 1, [0, 12]),
  ("lep_b_pt" , "p_{T}(b) (GeV)", 10, None),
  ("had_b_pt" , "p_{T}(b) (GeV)", 10, None),
  ("lep_pt"   , "p_{T}(l) (GeV)", 10, None),
  #("Wlep_mass", "m_{W}(lep) (GeV)", 10, None),
  ("Whad_mass", "m_{W}(had) (GeV)", 10, None),
  ("Whad_DR"  , "#DeltaR(jj) W_{had} (GeV)", 1, [0,7]),
  ("Whad_pt"  , "p_{T}(W_{had}) (GeV)", 10, None),
  #("Whad_leading_DR", "#DeltaR(W_{had}, leading jet)"   , 1, [0,7]),
  #("Whad_sublead_DR", "#DeltaR(W_{had}, subleading jet)", 1, [0,7]),
  ("thad_mass", "m_{t}(had) (GeV)", 10, None),
  #("nu_chisq"         , "nu_chisq"         , 1, [0, 20]),
	#("nu_discriminant"	, "nu_discriminant"	 , 1, None),
	#("btag_discriminant", "btag_discriminant", 1, [-11, -2]),
	("mass_discriminant", "mass_discriminant", 2, None), #[5, 20]),
	("full_discriminant", "full_discriminant", 1, None),
]

shapes = set([
      'ncharged', 'nneutral', 'ntotal', 'mass_discriminant', 
      'full_discriminant', 'Whad_sublead_DR', 'Whad_pt', 
      'Whad_leading_DR', 'Whad_DR'
])

orders = [
   ## "nu_chisq"         , 
   ## "nu_discriminant"  , 
   ## "btag_discriminant", 
   "mass_discriminant", 
   ## "full_discriminant",
]

working_points = [
   "notag",
   ## "csvTight",
   ## "csvMedium",
   ## "csvLoose",
   ## "rndm10",
]

both_tag_binning = {
   'csvTight' : 100,
   'csvMedium': 10,
   'csvLoose' : 4,
   'rndm10'   : 4,
}

jet_categories = [
   "lead_tagged", "sublead_tagged", 
   "both_tagged", "both_untagged"
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

if not args.noplots:
   views_to_flow = filter(lambda x: 'ttJets' not in x and 'QCD' not in x, plotter.mc_samples)
   views_to_flow.append(plotter.tt_to_use)
   stack = plotting.HistStack()
   qcd_samples = [i for i in plotter.views if 'QCD' in i]

   for vtf in views_to_flow:
      histo = plotter.get_view(vtf).Get('cut_flow')
      plotter.keep.append(histo)
      stack.Add(
         histo
         )

   #QCD may not have all the bins filled, needs special care
   qcd_histo = histo.Clone()
   qcd_histo.Reset()
   for sample in qcd_samples:
      qcd_flow = plotter.get_view(sample).Get('cut_flow')
      qcd_histo = qcd_histo.decorate(
         **qcd_flow.decorators
         )
      qcd_histo.title = qcd_flow.title
      for sbin, qbin in zip(qcd_histo, qcd_flow):
         sbin.value += qbin.value
         sbin.error = quad.quad(sbin.error, qbin.error)
   stack.Add(qcd_histo)

   histo.Draw() #set the proper axis labels
   data = plotter.get_view('data').Get('cut_flow')
   smin = min(stack.min(), data.min(), 1.2)
   smax = max(stack.max(), data.max())
   histo.yaxis.range_user = smin*0.8, smax*1.2
   stack.Draw('same')
   data.Draw('same')
   plotter.add_legend([stack, data], False, entries=len(views_to_flow)+1)
   plotter.pad.SetLogy()
   plotter.add_ratio_plot(data, stack, ratio_range=0.4)
   plotter.lower_pad.SetLogy(False)
   #cut_flow.GetYaxis().SetRangeUser(1, 10**7)
   plotter.save('cut_flow', pdf=False)

   for order in orders:
      plotter.plot_mc_vs_data(
         'nosys', 'nvtx', sort=True,
         xaxis="# vtx", leftside=False, 
         show_ratio=True, ratio_range=0.5)
      plotter.save('nvtx', pdf=False)

      with io.root_open(os.path.join(plotter.outputdir, 'vtxs.root'), 'recreate') as out:
         vtx_stack = plotter.make_stack(folder='nosys').Get('nvtx')
         vtx_data = plotter.get_view('data').Get('nosys/nvtx')
         vtx_data.Scale(1./vtx_data.Integral())
         ratio = sum(vtx_stack.hists)
         ratio.Scale(1./ratio.Integral())
         for rbin, dbin in zip(ratio, vtx_data):
            rbin.value = dbin.value/rbin.value if rbin.value else 0.
            rbin.error = 0.
         out.WriteTObject(ratio, 'vtx_reweigt')

      for wpoint in working_points:
        for jet_cat in jet_categories:
          if wpoint == 'notag' and jet_cat <> "both_untagged": continue
          base = os.path.join(order, wpoint, jet_cat)
          for jtype in jet_types:
             folder = os.path.join(base, jtype)
             plotter.set_subdir(folder)
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
                plotter.save('_'.join((jtype,var)), pdf=False)
                if var in shapes:
                   plotter.plot_mc_shapes(
                      folder, var, rebin,
                      xaxis=axis % jtype, leftside=False,
                      xrange=x_range)
                   plotter.save('shape_'+'_'.join((jtype,var)), pdf=False)
             plotter.bake_pie(folder)

          plotter.set_subdir(base)
          for var, axis, rebin, x_range in variables:
             folder = os.path.join('nosys', base)
             plotter.plot_mc_vs_data(
                folder, var, rebin, sort=True,
                xaxis=axis, leftside=False,
                xrange=x_range, show_ratio=True, 
                ratio_range=0.5)
             plotter.save(var, pdf=False)
             if var in shapes:
                if var == "mass_discriminant" and "_tagged" in jet_cat: rebin = 4
                plotter.plot_mc_shapes(
                   folder, var, rebin,
                   xaxis=axis, leftside=False,
                   xrange=x_range)
                plotter.save('shape_%s' % var, pdf=False)
                if var <> "mass_discriminant": continue
                plotter.plot_mc_shapes(
                   folder, var, rebin, normalize=True,
                   xaxis=order.replace('_', ' '), leftside=False,
                   xrange=x_range)
                plotter.save('normshape_%s' % order, pdf=False)

                plotter.plot_mc_shapes(
                   folder, var, rebin, normalize=True, show_err=True,
                   xaxis=order.replace('_', ' '), leftside=False,
                   xrange=x_range, filt=['ttJets_allRight', 'ttJets_rightHad', 'ttJets_rightWHad'])
                plotter.save('normshape_rightwhad_%s' % order, pdf=False)

                plotter.plot_mc_shapes(
                   folder, var, rebin, normalize=True, show_err=True,
                   xaxis=order.replace('_', ' '), leftside=False,
                   xrange=x_range, filt=['ttJets_rightWLep', 'ttJets_semiWrong'])
                plotter.save('normshape_wrongwhad_%s' % order, pdf=False)

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

if not args.noshapes:
   plotter.set_subdir("")
   info = plotter.write_summary_table(orders, working_points)
   cwd = os.getcwd()
   for order in orders:
      for wpoint in working_points:
         if wpoint == "notag": continue
         plotter.set_subdir(
            os.path.join(order, wpoint)
            )
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
               rebin=2
               )
         plotter.add_systematics(args.noBBB)

         plotter.card.add_systematic('lead_cfrac', 'param', '', '', info[order]['leadCFrac'], info[order]['leadCFrac_err'])
         plotter.card.add_systematic('sub_cfrac' , 'param', '', '', info[order]['subCFrac'] , info[order]['subCFrac_err'] )
         plotter.card.add_systematic('mc_lead_charm_eff', 'param', '', '', info[order]['working_points'][wpoint]['lead_charmEff'], info[order]['working_points'][wpoint]['lead_charmEff_err'])
         plotter.card.add_systematic('mc_lead_light_eff', 'param', '', '', info[order]['working_points'][wpoint]['lead_lightEff'], info[order]['working_points'][wpoint]['lead_lightEff_err'])
         plotter.card.add_systematic('mc_sub_charm_eff' , 'param', '', '', info[order]['working_points'][wpoint]['sub_charmEff'] , info[order]['working_points'][wpoint]['sub_charmEff_err'] )
         plotter.card.add_systematic('mc_sub_light_eff' , 'param', '', '', info[order]['working_points'][wpoint]['sub_lightEff'] , info[order]['working_points'][wpoint]['sub_lightEff_err'] )
         plotter.card.add_systematic('signal_norm', 'param', '', '', info[order]['nevts'], info[order]['nevts_err'])
         if args.noLightFit:
            plotter.card.add_systematic('lightSF', 'param', '', '', 1.00, 0.1)

         plotter.save_card('datacard')

         with open(os.path.join(plotter.outputdir, 'make_workspace.sh'), 'w') as f:
            #set_trace()
            f.write("sed -i 's|$MASS||g' datacard.txt\n")
            f.write(r"sed -i 's|\x1b\[?1034h||g' datacard.txt"+'\n')
            f.write("echo creating workspace\n")
            workspace_opts = '--PO fitLightEff=False' if args.noLightFit else ''
            workspace_opts +=' --PO inclusive=True' if args.inclusive else ''
            f.write(
               'text2workspace.py datacard.txt -P HiggsAnalysis.CombinedLimit'
               '.CTagEfficiencies:ctagEfficiency -o fitModel.root %s\n' % workspace_opts
               )
            f.write("echo 'running Multi-dimensional fit with Profile-Likelyhood errors on Asimov'\n")
            f.write("combine fitModel.root -M MultiDimFit --algo=singles "
                    "--setPhysicsModelParameterRanges charmSF=0,2:lightSF=0,2 "
                    "-t -1 --expectSignal=1 > asimovFit.log\n")
            f.write("mv higgsCombineTest.MultiDimFit.mH120.root asimovFit.root\n")
            f.write("echo 'running 2D Scan on Asimov'\n")
            f.write("combine -M MultiDimFit fitModel.root --algo=grid --points=400 "
                    "--setPhysicsModelParameterRanges charmSF=0.05,2.05:lightSF=0.05,2.05 "
                    "-t -1 --expectSignal=1 > asimovScan.log\n")
            f.write("mv higgsCombineTest.MultiDimFit.mH120.root asimovScan.root\n")
            f.write("echo 'running Multi-dimensional with Profile-Likelyhood errors'\n")
            f.write("combine fitModel.root -M MultiDimFit --algo=singles "
                    "--setPhysicsModelParameterRanges charmSF=0,2:lightSF=0,2 > dataFit.log\n")
            f.write("mv higgsCombineTest.MultiDimFit.mH120.root DataFit.root\n")
            ## f.write("echo 'running 2D likelyhood scan'\n")
            ## f.write("combine -M MultiDimFit fitModel.root --algo=grid --points=900 "
            ##         "--setPhysicsModelParameterRanges charmSF=0,2:lightTagScale=0,2 > dataScan.log\n")
            ## f.write("mv higgsCombineTest.MultiDimFit.mH120.root DataScan.root\n")

