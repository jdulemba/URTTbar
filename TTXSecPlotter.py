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

parser = ArgumentParser()
parser.add_argument('--noplots', dest='noplots', action='store_true',
                    help='skip plot making')
parser.add_argument('--noshapes', dest='noshapes', action='store_true',
                    help='skip shape making')
opts = parser.parse_args()

class TTXSecPlotter(Plotter):
   def __init__(self):
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
         files, lumis, outdir, styles, None, 10000
         )
      self.jobid = jobid

      self.views['ttJets_rightAssign'] = {
         'view' : self.create_tt_subsample(
            'semilep_visible_right', 
            'tt, right cmb',
            '#6666b3'
            )
         }
      self.views['ttJets_rightThad'] = {
         'view' : self.create_tt_subsample(
            'semilep_right_thad', 
            'tt, right t_{h}',
            ),
         }
      self.views['ttJets_rightTlep'] = {
         'view' : self.create_tt_subsample(
            'semilep_right_tlep', 
            'tt, right t_{l}',
            '#cccce6'
            )
         }
      self.views['ttJets_wrongAssign'] = {
         'view' : self.create_tt_subsample(
            'semilep_wrong', 
            'tt, wrong cmb',
            '#88a7c4'
            )
         }
      self.views['ttJets_other'] = {
         'view' : self.create_tt_subsample(
            'other_tt_decay', 
            'Other tt decay',
            '#668db3',
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
         'ttJets_other'
         ]

      self.card_names = [
         'vjets',
         'single_top',
         'tt_right',
         'tt_right_th',
         'tt_right_tl',
         'tt_wrong',
         'tt_other'
         ]

      self.systematics = {
         'lumi' : {
            'type' : 'lnN',
            'samples' : ['*'],
            'categories' : ['*'],
            'value' : 1.05,
            },
         'JES' : {
            'samples' : ['*'],
            'categories' : ['*'],
            'type' : 'shape',
            '+' : lambda x: x.replace('nosys', 'jes_up'),
            '-' : lambda x: x.replace('nosys', 'jes_down'),
            }
         }
      self.card = None
      self.binning = {}

   def create_tt_subsample(self, subdir, title, color='#9999CC'):
      return views.StyleView(
         views.TitleView(
            views.SubdirectoryView(
               self.views['ttJets_pu30']['view'],
               subdir
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
      pass

   def write_shapes(self, folder, variable, xbinning, ybinning, 
                    category_template='Bin%i', slice_along='X'):
      if not self.card: self.card = DataCard('tt_*')
      #keep it there for systematics
      mc_views = dict(
            (name, view) for name, view in zip(
               self.card_names,
               self.mc_views(rebin=[xbinning, ybinning])
               )
            )

      path = os.path.join(folder, variable)
      mc_hists2D = dict(
         (name, view.Get(path)) for name, view in mc_views.iteritems()
         )

      category_axis = 'GetNbinsX' if slice_along == 'Y' else 'GetNbinsY'
      nbins = getattr(mc_hists2D.values()[0], category_axis)()
      fake_data = sum(i for i in mc_hists2D.values())
      sliced_axis = getattr(
         fake_data, 
         'GetXaxis' if slice_along == 'Y' else 'GetYaxis'
         )()
      binning = {}

      for idx in range(nbins):
         self.card.add_category(category_template % idx)
         category = self.card[category_template % idx]
         logging.debug(
            'slicing bin %i from %.2f to %.2f' % (
               idx+1,
               sliced_axis.GetBinLowEdge(idx+1),
               sliced_axis.GetBinUpEdge(idx+1),
               )
            )
         binning[category_template % idx] = {
            'var' : var,
            'idx' : idx,
            'low_edge' : sliced_axis.GetBinLowEdge(idx+1),
            'up_edge' : sliced_axis.GetBinUpEdge(idx+1)
            }

         for name, hist in mc_hists2D.iteritems():
            category[name] = slice_hist(hist, idx+1, axis=slice_along)
         category['data_obs'] = slice_hist(fake_data, idx+1, axis=slice_along)
         integral = category['data_obs'].Integral()
         if integral != 0:
            int_int = float(int(integral))
            category['data_obs'].Scale(int_int/integral)
      self.binning[var] = binning

##################
#  DEFINITIONS
##################

plotter = TTXSecPlotter()

pt_binning = Struct(
   gen = [0., 40., 75., 105., 135., 170., 220., 300., 1000.],
   reco = [0., 40., 60., 75., 90., 105., 120., 135., 150., 170., 195., 220., 260., 300., 500., 1000.],
   )
   
## eta_binning = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.3, 2.8, 8.0]
## mass_binning = [250., 350., 370., 390., 410., 430., 450., 470., 490., 510., 530., 550., 575., 600., 630., 670., 720., 800., 900, 5000.]
## y_binning = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 3.]
## ttpt_binning = [0., 20., 30., 40., 50., 60., 70., 90., 110., 140., 180., 250., 1000.]

##################
#     PLOTS
##################
if not opts.noplots:
   to_plot = [
      ('all_ptthad' , pt_binning.reco) ,
      ('all_pttlep' , pt_binning.reco) ,
      ('all_lep_pt'  , 4),
      ('all_tlep_eta', 4),
      ('all_thad_eta', 4),
      ("all_ttm"     , 4),
      ("all_tty"     , 4),
      ("all_ttpt"    , 4),
      ("all_costhetastar", 4),
      ("all_njet", 4),
      ]

   for var, rebin in to_plot:
      plotter.plot_mc_vs_data('', var, leftside=False, rebin=rebin)
      plotter.save(var, pdf=False)

   plotter.plot_mc_vs_data(
      '', 'all_fullDiscr_ptthad', leftside=False, 
      preprocess=lambda x: urviews.ProjectionView(x, 'X', [0, 1000])
      )
   plotter.save('fullDiscr_from_projection', pdf=False)
   
   previous = pt_binning.reco[0]
   for idx, ptbin in enumerate(pt_binning.reco[1:]):
      plotter.plot_mc_vs_data(
         '', 'all_fullDiscr_ptthad', leftside=False, rebin = [pt_binning.reco, [4]],
         preprocess=lambda x: urviews.ProjectionView(x, 'X', [previous, ptbin])
         )
      previous = ptbin
      plotter.save('fullDiscr_slice_%i' % idx, pdf=False)
   
   #TEST CASE FOR MARIO
   plotter.plot_mc_vs_data(
         '', 'all_fullDiscr_ptthad', leftside=False, rebin = [pt_binning.reco, [160]],
         preprocess=lambda x: urviews.ProjectionView(x, 'X', [previous, ptbin])
         )
   plotter.save('mario_test', pdf=False, dotroot=True)

##################
#     CARDS
##################
if not opts.noshapes:
   full_discr_binning = [5]#range(-15, 16)
   to_fit = [
      ("ptthad"	, pt_binning),
      ("pttlep"	, pt_binning),
   ##    ("etathad", eta_binning),
   ##    ("etatlep", eta_binning),
   ##    ("ttm"		, mass_binning),
   ##    ("tty"		, y_binning),
   ##    ("ttpt"   , ttpt_binning),
      ]
   discriminant = 'fullDiscr'
   
   for var, binning in to_fit:
      plotter.set_subdir(var)
      plotter.write_shapes(
         '', 'all_%s_%s' % (discriminant, var), full_discr_binning, binning.reco
         ) 
      plotter.card.add_systematic('lumi', 'lnN', '.*', '[^t]+.*', 1.05)
      #plotter.card.add_bbb_systematics(['*'], ['*'])
      plotter.save_card(var)
   
   #Migration matrices
   plotter.set_subdir('')
   fname = os.path.join(plotter.outputdir, 'migration_matrices.root')

   to_fit = [("hadtop_pt" , pt_binning)]

   with io.root_open(fname, 'recreate') as mfile:
      tt_view = plotter.get_view('ttJets_pu30')
      for var, binning in to_fit:
         mfile.mkdir('ptthad').cd() ##FIXME var) 
         matrix_view = plotter.rebin_view(tt_view, [pt_binning.gen, pt_binning.reco])
         mig_matrix = matrix_view.Get('RECO/truth_%s_matrix_fiducialtight' % var)
         mig_matrix.SetName('migration_matrix') ##FIXME var)
         mig_matrix.Write()
         thruth_distro = mig_matrix.ProjectionX() 
         #plotter.rebin_view(tt_view, pt_binning.gen).Get('TRUTH/truth_response_%s_truth' % var)
         thruth_distro.SetName('true_distribution')
         thruth_distro.Write()

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


   logging.info('file: %s written' % fname)
      
