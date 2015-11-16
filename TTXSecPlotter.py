from pdb import set_trace
from URAnalysis.PlotTools.Plotter import Plotter
import URAnalysis.PlotTools.views as urviews
import os, glob, sys, logging, ROOT, rootpy, itertools
from URAnalysis.Utilities.datacard import DataCard
from URAnalysis.Utilities.roottools import slice_hist
from styles import styles
import rootpy.plotting as plotting
views = plotting.views
import rootpy.io as io
from array import array
from pdb import set_trace
from fnmatch import fnmatch
import URAnalysis.Utilities.prettyjson as prettyjson
import URAnalysis.Utilities.quad as quad
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
from argparse import ArgumentParser
from URAnalysis.Utilities.struct import Struct
import re

class TTXSecPlotter(Plotter):
   def __init__(self, ttbar_to_use='ttJets', lumi=None):
      self.ttbar_to_use = ttbar_to_use
      jobid = os.environ['jobid']
      files = glob.glob('results/%s/ttxs/*.root' % jobid)
      if len(files) == 0:
         logging.info('No file found, must be Otto\'s format...')
         all_files = glob.glob('results/%s/ttxs/*/*.root' % jobid)
         files = {}
         for ifile in all_files:
            base = os.path.basename(ifile)
            if base not in files: 
               files[os.path.basename(ifile)] = []
            files[os.path.basename(ifile)].append(
               (os.path.basename(os.path.dirname(ifile)), io.root_open(ifile))
               )
         files = dict((i, urviews.MultifileView(**dict(j))) for i, j in files.iteritems())

      logging.debug('files found %s' % files.__repr__())
      lumis = glob.glob('inputs/%s/*.lumi' % jobid)
      logging.debug('lumi files found %s' % lumis.__repr__())
      
      outdir= 'plots/%s/ttxsec' % jobid
      super(TTXSecPlotter, self).__init__(
         files, lumis, outdir, styles, None, lumi, lumi_scaling=0.5
         )
      self.jobid = jobid

      self.mc_samples = []
      self.card_names = {}
      self.systematics= {}


      self.card = None
      self.binning = {}
      self.initialized = False
      self.merged_leptons=False

   def initviews(self):
      if self.initialized: return
      for sample in self.views:
         if not sample.startswith('ttJets'):
            self.views[sample]['view'] = views.SubdirectoryView(self.views[sample]['view'], 'RECO')
      self.initialized = True

   def merge_leptons(self):
      if self.merged_leptons: return
      self.initviews()
      for sample in self.views:
         self.views[sample]['view'] = views.SumView(
            views.SubdirectoryView(self.views[sample]['view'], 'muons'),
            views.SubdirectoryView(self.views[sample]['view'], 'electrons'),
            )
      self.merged_leptons = True

   def cut_flow(self):
      views_to_flow = filter(lambda x: 'ttJets' not in x and 'QCD' not in x, self.mc_samples)
      views_to_flow.append(self.ttbar_to_use)
      stack = plotting.HistStack()
      self.keep.append(stack)
      qcd_samples = [i for i in self.views if 'QCD' in i]

      for vtf in views_to_flow:
         histo = self.get_view(vtf).Get('cut_flow')
         print vtf, len(histo)
         self.keep.append(histo)
         stack.Add(
            histo
            )

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
      stack.Add(qcd_histo)
      self.keep.append(qcd_histo)

      histo.Draw() #set the proper axis labels
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

   def create_tt_subsample(self, subdirs, title, color='#9999CC', relative_scale=itertools.repeat(1.)):
      return views.StyleView(
         views.TitleView(
            views.SumView(
               *[
                  views.ScaleView(
                     views.SubdirectoryView(
                        self.views[self.ttbar_to_use]['view'],
                        subdir
                        ),
                     scale
                     )
                  for subdir, scale in zip(subdirs, relative_scale)
                  ]
               ),
            title
            ),
         fillcolor = color,
         linecolor = color
         )

   def save_card(self, name):
      if not self.card:
         raise RuntimeError('There is no card to save!')
      self.card.normalize_signals()
      self.card.save(name, self.outputdir)
      self.card = None

      binning_file = '%s/%s.binning.json' % (self.outputdir, name)
      with open(binning_file, 'w') as json:
         json.write(prettyjson.dumps(self.binning))
      self.binning = {}
      logging.info('binning saved in %s' % binning_file)

   def add_systematics(self):
      if not self.card:
         raise ValueError('The card is not defined!')
      for sys_name, info in self.systematics.iteritems():
         if info['type'] == 'stat': continue
         for category in info['categories']:
            for sample in info['samples']:
               self.card.add_systematic(
                  sys_name, 
                  info['type'], 
                  category, 
                  sample, 
                  info['value']
                  )
      #plotter.card.add_bbb_systematics('.*', '.*')

   def write_shapes(self, folder, var, variables, var_binning, disc_binning = lambda x, *args: 8, 
                    category_template='Bin%i', slice_along='X'):
      if not self.card: self.card = DataCard('tt_*')
      #keep it there for systematics
      card_views = {}
      binning = [[1],var_binning] if slice_along.lower() == 'x' else [var_binning,[1]]
      for name, samples in self.card_names.iteritems():
         card_views[name] = self.rebin_view(
            views.SumView(
               *[self.get_view(i) for i in samples]
                ),
            binning
            )

      data_is_fake = False
      if 'data' not in self.views:
         card_views['data_obs'] = views.SumView(*card_views.values())
         data_is_fake = True
      else:
         card_views['data_obs'] = self.rebin_view(self.get_view('data'), binning)

      paths = [os.path.join(folder, i) for i in variables]
      card_hists2D = {}
      for name, view in card_views.iteritems():
         card_hists2D[name] = sum(view.Get(i) for i in paths)
      sys_hists2D = {}

      category_axis = 'GetNbinsX' if slice_along == 'Y' else 'GetNbinsY'
      nbins = getattr(card_hists2D.values()[0], category_axis)()
      fake_data = sum(i for i in card_hists2D.values())
      sliced_axis = getattr(
         fake_data, 
         'GetXaxis' if slice_along == 'Y' else 'GetYaxis'
         )()
      binning = {}

      for idx in range(nbins):
         category_name = category_template % idx
         self.card.add_category(category_name)
         category = self.card[category_name]
         logging.debug(
            'slicing bin %i from %.2f to %.2f' % (
               idx+1,
               sliced_axis.GetBinLowEdge(idx+1),
               sliced_axis.GetBinUpEdge(idx+1),
               )
            )
         binning[category_name] = {
            'var' : var,
            'idx' : idx,
            'low_edge' : sliced_axis.GetBinLowEdge(idx+1),
            'up_edge' : sliced_axis.GetBinUpEdge(idx+1)
            }

         #make 1D histograms
         card_hists1D = dict(
            (name, slice_hist(hist, idx+1, axis=slice_along))
            for name, hist in card_hists2D.iteritems()
            )
         discriminator_binning = disc_binning(idx, *card_hists1D.items())
         for name, hist in card_hists1D.iteritems():
            category[name] = hist.Rebin(discriminator_binning)
            if data_is_fake and name == 'data_obs':
               integral = category[name].Integral()
               if integral != 0:
                  int_int = float(int(integral))
                  category[name].Scale(int_int/integral)
               continue

            for sys_name, info in self.systematics.iteritems():
               if not any(re.match(i, category_name+'$') for i in info['categories']): continue
               if not any(re.match(i, name) for i in info['samples']): continue
               shift = 1.0
               if info['type'] == 'shape':
                  if sys_name not in sys_hists2D:
                     sys_hists2D[sys_name] = {}
                  if name not in sys_hists2D[sys_name]:
                     paths_up = [info['+'](i) for i in paths] 
                     paths_dw = [info['-'](i) for i in paths]
                     h2Dup = sum(card_views[name].Get(i) for i in paths_up)
                     h2Ddw = sum(card_views[name].Get(i) for i in paths_dw)
                     sys_hists2D[sys_name][name] = {
                        '+' : h2Dup,
                        '-' : h2Ddw,
                        }
                  else:
                     h2Dup = sys_hists2D[sys_name][name]['+']
                     h2Ddw = sys_hists2D[sys_name][name]['-']

                  hup = slice_hist( h2Dup, idx+1, 
                     axis=slice_along).Rebin(discriminator_binning)
                  hdw = slice_hist( h2Ddw, idx+1, 
                     axis=slice_along).Rebin(discriminator_binning)
                  if 'shape_only' in info and info['shape_only']:
                     hup_integral = hup.Integral()
                     hdw_integral = hdw.Integral()
                     hct_integral = category[name].Integral()
                     
                     if hup_integral:
                        hup.Scale(hct_integral/hup_integral)
                     if hdw_integral:
                        hdw.Scale(hct_integral/hdw_integral)

                  category['%s_%sUp'   % (name, sys_name)] = hup
                  category['%s_%sDown' % (name, sys_name)] = hdw
               if info['type'] == 'stat':
                  err = ROOT.Double()
                  integral = category[name].IntegralAndError(
                     1,
                     category[name].GetNbinsX(),
                     err
                     )
                  if integral:
                     rel_err = err/integral
                     self.card.add_systematic(
                        '%s_%s_%s' % (name, category_name, sys_name), 
                        'lnN', 
                        category_name+'$', 
                        name, 
                        1.+rel_err*info['multiplier']
                        )

      if var not in self.binning: 
         self.binning[var] = {}
      self.binning[var].update(binning)

   def plot_mc_shapes(self, folder, variable, rebin=1, xaxis='',
                      leftside=True, xrange=None, preprocess=None,
                      normalize=False, logy=False, show_err=False,
                      use_only=None, ratio_range=3):
      mc_views = []
      if not use_only:
         mc_views = self.mc_views(rebin, preprocess, folder)
      else:
         mc_views = [i for i, j in zip(self.mc_views(rebin, preprocess, folder), self.mc_samples) 
                     if j in use_only]

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
         if show_err:
            hist.drawstyle = 'pe'
            hist.markerstyle = 20
            hist.legendstyle = 'pe'
            hist.markercolor = hist.linecolor

         hist.GetXaxis().SetTitle(xaxis)
         hist.GetYaxis().SetTitle('counts' if not normalize else 'A. U.')
         if first:
            hist.Draw()
            first = False
         else:
            hist.Draw('same')
         
         self.keep.append(hist)
         m = max(m, max(bin.value for bin in hist))

      histos[0].GetYaxis().SetRangeUser(0., m*1.2)
      if xrange:
         histos[0].GetXaxis().SetRangeUser(*tuple(xrange))
      self.add_legend(histos, leftside)

      ref_idx = 0
      ytitle='ratio'
      if not use_only or 'ttJets_rightAssign' in use_only:
         ref_idx = self.mc_samples.index('ttJets_rightAssign')
         ytitle='bkg / signal'
      ref = histos[ref_idx]
      tests = histos[:ref_idx]+histos[ref_idx+1:]

      self.add_ratio_plot(
         ref, *tests, x_range=xrange,
         ratio_range=ratio_range, ytitle=ytitle)
