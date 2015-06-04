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
from TTXSecPlotter import TTXSecPlotter

##################
#  DEFINITIONS
##################

pt_binning = Struct(
   gen = [0., 40., 75., 105., 135., 170., 220., 300., 1000.],
   reco = [0., 40., 60., 75., 90., 105., 120., 135., 150., 170., 195., 220., 260., 300., 500., 1000.],
   #reco = [0., 120., 1000.],
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


##################
#     PLOTS
##################

plotter = TTXSecPlotter()

plotter.views['ttJets_rightAssign'] = {
   'view' : plotter.create_tt_subsample(
      ['semilep_visible_right'], 
      'tt, right cmb',
      '#5555ab'
      )
   }
plotter.views['ttJets_rightThad'] = {
   'view' : plotter.create_tt_subsample(
      ['semilep_right_thad'], 
      'tt, right t_{h}',
      '#aaaad5'
      ),
   }
plotter.views['ttJets_rightTlep'] = {
   'view' : plotter.create_tt_subsample(
      ['semilep_right_tlep'], 
      'tt, right t_{l}',
      '#ab5555'
      )
   }
plotter.views['ttJets_wrongAssign'] = {
   'view' : plotter.create_tt_subsample(
      ['semilep_wrong'], 
      'tt, wrong cmb',
      '#d5aaaa'
      )
   }
plotter.views['ttJets_other'] = {
   'view' : plotter.create_tt_subsample(
      ['other_tt_decay'], 
      'Other tt decay',
      '#668db3',
      )
   }
plotter.views['ttJets_wrong'] = {
   'view' : urviews.MultifileView(
      **{'' : plotter.create_tt_subsample(
            ['semilep_wrong', 
             'other_tt_decay',
             'semilep_right_tlep',],
            'tt, wrong cmb',
            '#ab5555'
            ),
         'otherTT_ratio_up' : plotter.create_tt_subsample(
            ['semilep_wrong', 
             'other_tt_decay',
             'semilep_right_tlep',],
            'tt, wrong cmb',
            '#ab5555',
            relative_scale=[1., 1.5, 1.]
            ),
         'otherTT_ratio_down' : plotter.create_tt_subsample(
            ['semilep_wrong', 
             'other_tt_decay',
             'semilep_right_tlep',],
            'tt, wrong cmb',
            '#ab5555',
            relative_scale=[1., 0.5, 1.]
            )
         }
      )
   }
            
plotter.mc_samples = [
   '[WZ]Jets',
   'single*',
   'ttJets_rightAssign',
   
   'ttJets_rightThad',
   'ttJets_rightTlep',
   'ttJets_wrongAssign',
   
   #'ttJets_wrong',
   'ttJets_other'
   ]

plotter.card_names = {
   'vjets' : ['[WZ]Jets'],
   'single_top' : ['single*'],
   'tt_right' : ['ttJets_rightAssign'],
   'tt_wrong' : ['ttJets_wrong'],
   #'tt_other' : ['ttJets_other'],
   'only_thad_right' : ['ttJets_rightThad']
   }

plotter.systematics = {
   'lumi' : {
      'type' : 'lnN',
      'samples' : ['(?!tt_).*'],
      'categories' : ['.*'],
      'value' : 1.05,
      },
   'otherTT_ratio' : {
      'type' : 'shape',
      'samples' : ['tt_wrong'],
      'categories' : ['.*'],
      '+' : lambda x: 'otherTT_ratio_up/%s' % x,
      '-' : lambda x: 'otherTT_ratio_down/%s' % x,
      'value' : 1.00,
      'shape_only' : True,
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

plotter.initviews()
   

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
         '', var, 'all_%s_%s' % (discriminant, var), full_discr_binning, binning.reco
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
      
