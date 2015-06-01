from pdb import set_trace
from URAnalysis.PlotTools.Plotter import Plotter
import URAnalysis.PlotTools.views as urviews
import os, glob, sys, logging, ROOT, rootpy
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

##################
#  DEFINITIONS
##################

pt_binning = Struct(
   gen = [0., 40., 75., 105., 135., 170., 220., 300., 1000.],
   #reco = [0., 40., 60., 75., 90., 105., 120., 135., 150., 170., 195., 220., 260., 300., 500., 1000.],
   reco = [0., 120., 1000.],
   )
discriminant =  'massDiscr'
phase_space = 'fiducialtight'
full_discr_binning = [4]#range(-15, 16)
   
## eta_binning = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.8, 8.0]
## mass_binning = [250., 350., 370., 390., 410., 430., 450., 470., 490., 510., 530., 550., 575., 600., 630., 670., 720., 800., 900, 5000.]
## y_binning = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 3.]
## ttpt_binning = [0., 20., 30., 40., 50., 60., 70., 90., 110., 140., 180., 250., 1000.]

parser = ArgumentParser()
parser.add_argument('--noplots', dest='noplots', action='store_true',
                    help='skip plot making')
parser.add_argument('--noshapes', dest='noshapes', action='store_true',
                    help='skip shape making')
opts = parser.parse_args()

class TTXSecPlotter(Plotter):
   def __init__(self, ttbar_to_use='ttJets_pu30'):
      self.ttbar_to_use = ttbar_to_use
      jobid = os.environ['jobid']
      files = glob.glob('results/%s/ttbarxsec/*.root' % jobid)
      if len(files) == 0:
         logging.info('No file found, must be Otto\'s format...')
         all_files = glob.glob('results/%s/ttbarxsec/*/*.root' % jobid)
         files = {}
         for ifile in all_files:
            base = os.path.basename(ifile)
            if base not in files: 
               files[os.path.basename(ifile)] = []
            files[os.path.basename(ifile)].append(
               (os.path.basename(os.path.dirname(ifile)), io.root_open(ifile))
               )
         files = dict((i, urviews.MultifileView(*j)) for i, j in files.iteritems())

      logging.debug('files found %s' % files.__repr__())
      lumis = glob.glob('inputs/%s/*.lumi' % jobid)
      logging.debug('lumi files found %s' % lumis.__repr__())
      
      outdir= 'plots/%s/ttxsec' % jobid
      super(TTXSecPlotter, self).__init__(
         files, lumis, outdir, styles, None, 10**3
         )
      self.jobid = jobid

      self.views['ttJets_rightAssign'] = {
         'view' : self.create_tt_subsample(
            ['semilep_visible_right'], 
            'tt, right cmb',
            '#5555ab'
            )
         }
      self.views['ttJets_rightThad'] = {
         'view' : self.create_tt_subsample(
            ['semilep_right_thad'], 
            'tt, right t_{h}',
            '#aaaad5'
            ),
         }
      self.views['ttJets_rightTlep'] = {
         'view' : self.create_tt_subsample(
            ['semilep_right_tlep'], 
            'tt, right t_{l}',
            '#ab5555'
            )
         }
      self.views['ttJets_wrongAssign'] = {
         'view' : self.create_tt_subsample(
            ['semilep_wrong'], 
            'tt, wrong cmb',
            '#d5aaaa'
            )
         }
      self.views['ttJets_other'] = {
         'view' : self.create_tt_subsample(
            ['other_tt_decay'], 
            'Other tt decay',
            '#668db3',
            )
         }

      self.views['ttJets_wrong'] = {
         'view' : self.create_tt_subsample(
            ['semilep_wrong', 
            'semilep_right_thad', 
            'semilep_right_tlep',],
            'tt, wrong cmb',
            '#ab5555'
            )
         }


      for sample in self.views:
         if not sample.startswith('ttJets'):
            self.views[sample]['view'] = views.SubdirectoryView(self.views[sample]['view'], 'RECO')

      self.mc_samples = [
         '[WZ]Jets',
         'single*',
         'ttJets_rightAssign',

         'ttJets_rightThad',
         'ttJets_rightTlep',
         'ttJets_wrongAssign',

         #'ttJets_wrong',
         'ttJets_other'
         ]

      self.card_names = {
         'vjets' : ['[WZ]Jets'],
         'single_top' : ['single*'],
         'tt_right' : ['ttJets_rightAssign'],
         'tt_wrong' : ['ttJets_rightTlep', 'ttJets_wrongAssign'],
         'tt_other' : ['ttJets_other'],
         'only_thad_right' : ['ttJets_rightThad']
         }

      self.systematics = {
         'lumi' : {
            'type' : 'lnN',
            'samples' : ['.*'],
            'categories' : ['.*'],
            'value' : 1.05,
            },
         ## 'JES' : {
         ##    'samples' : ['*'],
         ##    'categories' : ['*'],
         ##    'type' : 'shape',
         ##    '+' : lambda x: x.replace('nosys', 'jes_up'),
         ##    '-' : lambda x: x.replace('nosys', 'jes_down'),
         ##    'value' : 1.00,
         ##    },
         'MCStat' : {
            'samples' : ['only_thad_right'],
            'categories' : ['.*'],
            'type' : 'stat',
            'multiplier' : 4.,
            }
         }
      self.card = None
      self.binning = {}

   def create_tt_subsample(self, subdirs, title, color='#9999CC'):
      return views.StyleView(
         views.TitleView(
            views.SumView(
               *[views.SubdirectoryView(
                  self.views[self.ttbar_to_use]['view'],
                  subdir
                  ) for subdir in subdirs]
               ),
            title
            ),
         fillcolor = color,
         linecolor = color
         )

   def save_card(self, name):
      if not self.card:
         raise RuntimeError('There is no card to save!')
      self.card.save(name, self.outputdir)
      self.card = None

      binning_file = '%s/%s.binning.json' % (self.outputdir, name)
      with open(binning_file, 'w') as json:
         json.write(prettyjson.dumps(self.binning))
      self.binning = {}
      logging.info('binning saved in %s' % binning_file)

   def add_systematics(self):
      if not plotter.card:
         raise ValueError('The card is not defined!')
      for sys_name, info in self.systematics.iteritems():
         if info['type'] == 'stat': continue
         for category in info['categories']:
            for sample in info['samples']:
               plotter.card.add_systematic(
                  sys_name, 
                  info['type'], 
                  category, 
                  sample, 
                  info['value']
                  )
      #plotter.card.add_bbb_systematics('.*', '.*')

   def write_shapes(self, folder, variable, xbinning, ybinning, 
                    category_template='Bin%i', slice_along='X'):
      if not self.card: self.card = DataCard('tt_*')
      #keep it there for systematics
      card_views = {}
      for name, samples in self.card_names.iteritems():
         card_views[name] = self.rebin_view(
            views.SumView(
               *[self.get_view(i) for i in samples]
                ),
            [xbinning, ybinning]
            )

      data_is_fake = False
      if 'data_obs' not in card_views:
         card_views['data_obs'] = views.SumView(*card_views.values())
         data_is_fake = True

      path = os.path.join(folder, variable)
      card_hists2D = dict(
         (name, view.Get(path)) for name, view in card_views.iteritems()
         )

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

         for name, hist in card_hists2D.iteritems():
            category[name] = slice_hist(hist, idx+1, axis=slice_along)
            if data_is_fake and name == 'data_obs':
               integral = category[name].Integral()
               if integral != 0:
                  int_int = float(int(integral))
                  category[name].Scale(int_int/integral)
               continue

            for sys_name, info in self.systematics.iteritems():
               if not any(re.match(i, category_name) for i in info['categories']): continue
               if not any(re.match(i, name) for i in info['samples']): continue
               shift = 1.0
               if info['type'] == 'shape':
                  raise RuntimeError('Still to implement!')
                  #this will not work, since we handle directly 2D histograms
                  #we need to change the concept
                  paths_up = [info['+'](path) for path in paths] 
                  paths_dw = [info['-'](path) for path in paths]
                  hup = sum(view.Get(i) for i in paths_up)
                  hdw = sum(view.Get(i) for i in paths_dw)
                  category['%s_%sUp'   % (name, sys_name)] = hup
                  category['%s_%sDown' % (name, sys_name)] = hdw
               if info['type'] == 'stat':
                  err = ROOT.Double()
                  integral = category[name].IntegralAndError(
                     1,
                     category[name].GetNbinsX(),
                     err
                     )
                  rel_err = err/integral
                  self.card.add_systematic(
                     '%s_%s_%s' % (name, category_name, sys_name), 
                     'lnN', 
                     category_name, 
                     name, 
                     1.+rel_err*info['multiplier']
                     )

      self.binning[var] = binning

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
         m = max(m,max(bin for bin in hist))
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

##################
#     PLOTS
##################
plotter = TTXSecPlotter()

if not opts.noplots:
   to_plot = [
      ('all_ptthad' , pt_binning.reco, 'p_{T}(t_{had})') ,
      ('all_pttlep' , pt_binning.reco, 'p_{T}(t_{lep})') ,
      ('all_lep_pt'  , 4, 'p_{T}(l)'),
      ('all_tlep_eta', 4, '#eta(t_{lep})'),
      ('all_thad_eta', 4, '#eta(t_{had})'),
      ("all_ttm"     , 4, 'm(t#bar{t})'),
      ("all_tty"     , 4, 'y(t#bar{t})'),
      ("all_ttpt"    , 4, 'p_{T}(t#bar{t})'),
      ("all_costhetastar", 4, ''),
      ("all_njet", 4, '# of jets'),
      ('all_%s' % discriminant, 2, 'discriminant'),
      ]

   for var, rebin, xaxis in to_plot:
      plotter.plot_mc_vs_data('', var, leftside=False, rebin=rebin, xaxis=xaxis)
      plotter.save(var, pdf=False)

   plotter.plot_mc_vs_data(
      '', 'all_%s_ptthad' % discriminant, leftside=False, 
      preprocess=lambda x: urviews.ProjectionView(x, 'X', [0, 1000])
      )
   plotter.save('%s_from_projection' % discriminant, pdf=False)
   
   previous = pt_binning.reco[0]
   plotter.set_subdir('ptthad/slices')
   for idx, ptbin in enumerate(pt_binning.reco[1:]):
      plotter.plot_mc_vs_data(
         '', 'all_%s_ptthad' % discriminant, leftside=False, rebin = [pt_binning.reco, full_discr_binning],
         preprocess=lambda x: urviews.ProjectionView(x, 'X', [previous, ptbin])
         )
      plotter.save('%s_slice_%i' % (discriminant, idx), pdf=False)
   
      plotter.plot_mc_shapes(
         '', 'all_%s_ptthad' % discriminant, rebin = [pt_binning.reco, full_discr_binning], 
         xaxis='discriminant', leftside=False, normalize=True, show_err=True, xrange=(1,9),
         preprocess=lambda x: urviews.ProjectionView(x, 'X', [previous, ptbin]))
      plotter.save('%s_slice%i_shape' % (discriminant, idx), pdf=False)
      
      plotter.plot_mc_shapes(
         '', 'all_%s_ptthad' % discriminant, rebin = [pt_binning.reco, full_discr_binning], 
         xaxis='discriminant', leftside=False, normalize=True, show_err=True, xrange=(1,9),
         preprocess=lambda x: urviews.ProjectionView(x, 'X', [previous, ptbin]),
         use_only=set(['ttJets_rightTlep', 'ttJets_wrongAssign', 'ttJets_other']),
         ratio_range=1)
      plotter.save('%s_bkgslice%i_shape' % (discriminant, idx), pdf=False)
      previous = ptbin

   plotter.set_subdir('')

   plotter.plot_mc_shapes(
      '', 'all_%s' % discriminant, rebin=2, xaxis='discriminant',
      leftside=False, normalize=True, show_err=True, xrange=(1,9))
   plotter.save('%s_full_shape' % (discriminant), pdf=False)

   plotter.plot_mc_shapes(
      '', 'all_%s' % discriminant, rebin=2, xaxis='discriminant',
      leftside=False, normalize=True, show_err=True, xrange=(1,9),
      use_only=set(['ttJets_rightTlep', 'ttJets_wrongAssign', 'ttJets_other']),
      ratio_range=1)
   plotter.save('%s_bkg_shape' % (discriminant), pdf=False)

##################
#     CARDS
##################
if not opts.noshapes:
   to_fit = [
      ("ptthad", 'ptthad', pt_binning),
   ##    ("pttlep", 'leptop_pt', pt_binning),
   ##    ("etathad", eta_binning),
   ##    ("etatlep", eta_binning),
   ##    ("ttm"		, mass_binning),
   ##    ("tty"		, y_binning),
   ##    ("ttpt"   , ttpt_binning),
      ]

   for var, _, binning in to_fit:
      plotter.set_subdir(var)
      plotter.write_shapes(
         '', 'all_%s_%s' % (discriminant, var), full_discr_binning, binning.reco
         ) 
      plotter.add_systematics()
      #plotter.card.add_systematic('lumi', 'lnN', '.*', '[^t]+.*', 1.05)
      plotter.save_card(var)
   
   #Migration matrices
   plotter.set_subdir('')
   fname = os.path.join(plotter.outputdir, 'migration_matrices.root')

   with io.root_open(fname, 'recreate') as mfile:
      for name, var, binning in to_fit:
         mfile.mkdir(name).cd() ##FIXME var) 
         matrix_path = 'RECO/truth_%s_matrix_%s' % (var, phase_space)
         tt_view = plotter.get_view(plotter.ttbar_to_use, 'unweighted_view')
         matrix_view_unscaled = plotter.rebin_view(tt_view, [pt_binning.gen, pt_binning.reco])
         mig_matrix_unscaled = matrix_view_unscaled.Get(matrix_path)
         mig_matrix_unscaled.SetName('migration_matrix')
         mig_matrix_unscaled.Write()

         #plotter.rebin_view(tt_view, pt_binning.gen).Get('TRUTH/truth_response_%s_truth' % var)
         tt_view = plotter.get_view(plotter.ttbar_to_use)
         matrix_view = plotter.rebin_view(tt_view, [pt_binning.gen, pt_binning.reco])
         mig_matrix = matrix_view.Get(matrix_path)
         mig_matrix.SetName('migration_matrix_scaled')
         mig_matrix.Write()

         ## thruth_distro_proj = plotter.get_view(plotter.ttbar_to_use).Get(matrix_path).ProjectionX()
         ## from array import array
         ## thruth_distro_proj = thruth_distro_proj.Rebin(
         ##    len(pt_binning.gen)-1, 
         ##    'true_distribution_fromProjection', 
         ##    array('d', pt_binning.gen)
         ##    )
         thruth_distro_proj = mig_matrix.ProjectionX()
         thruth_distro_proj.SetName('true_distribution_fromProjection')
         thruth_distro_proj.Write()

         reco_distro = mig_matrix.ProjectionY() 
         reco_distro.SetName('reco_distribution')
         reco_distro.Write()
         
         prefit_view = plotter.rebin_view(
            plotter.get_view('ttJets_rightAssign'), 
            pt_binning.reco
            )
         prefit_plot = prefit_view.Get('all_ptthad')
         prefit_plot.name = 'prefit_distribution'
         prefit_plot.Write()

         thruth_distro = plotter.rebin_view(
            tt_view,
            pt_binning.gen
            ).Get('RECO/truth_%s_gen_%s' % (var, phase_space))
         thruth_distro.SetName('true_distribution')
         thruth_distro.Write()         

   logging.info('file: %s written' % fname)
      
