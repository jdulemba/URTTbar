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
ROOT.gROOT.SetBatch(True)
from argparse import ArgumentParser
from URAnalysis.Utilities.struct import Struct
import re
from TTXSecPlotter import TTXSecPlotter

def run_module(**kwargs):
   ##################
   #  DEFINITIONS
   ##################
   opts = Struct(**kwargs)
   #when running optimization we do not want to
   #plot anything!
   if opts.optimize_binning and len(opts.optimize_binning):
      opts.noplots = True
   if opts.binning and len(opts.binning):
      opts.noplots = True

   discriminant =  'massDiscr'
   phase_space = 'fiducialtight'
   full_discr_binning = [4]#range(-15, 16)
      
   ## mass_binning = [250., 350., 370., 390., 410., 430., 450., 470., 490., 510., 530., 550., 575., 600., 630., 670., 720., 800., 900, 5000.]
   ## y_binning = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 3.]
   ## ttpt_binning = [0., 20., 30., 40., 50., 60., 70., 90., 110., 140., 180., 250., 1000.]

   vars_to_unfold = [
      Struct(
         var     = 'ptthad',
         binning = Struct(
            gen = [0., 60., 120., 180., 240., 300., 1000.],
            #gen = [0., 40., 75., 105., 135., 170., 220., 300., 1000.],
            reco = [30.0*i for i in range(11)]+[1000.],
            #reco = [0., 120., 1000.],
            ),
         xtitle  = 'p_{T}(t_{had})'
         ),
      Struct(
         var = 'pttlep',
         binning = Struct(
            gen = [0., 60., 120., 180., 240., 300., 1000.],
            #gen = [0., 40., 75., 105., 135., 170., 220., 300., 1000.],
            reco = [30.0*i for i in range(11)]+[1000.],
            ),
         xtitle = 'p_{T}(t_{lep})'
         ),
      Struct( 
         var = 'etatlep',
         binning = Struct(
            gen  = [0., 0.4, 0.8, 1.2, 1.6, 2.0, 8.0],
            #gen  = [0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.4,1.8, 2.3, 2.8, 8.0],
            reco = [0.2*i for i in range(11)]+[8.0],
            ),
         xtitle = '#eta(t_{lep})'
         ),
      Struct( 
         var = 'etathad',
         binning = Struct(
            gen  = [0., 0.4, 0.8, 1.2, 1.6, 2.0, 8.0],
            #gen  = [0., 0.2, 0.4, 0.6, 0.8, 1.0, 1.4,1.8, 2.3, 2.8, 8.0],
            reco = [0.2*i for i in range(11)]+[8.0],
            ),
         xtitle = '#eta(t_{had})'
         ),
   ]

   dir_postfix = ''
   if opts.optimize_binning and len(opts.optimize_binning):
      varname, prev, startbin, binning = tuple(opts.optimize_binning.split(':'))
      if not any(i.var == varname for i in vars_to_unfold):
         raise ValueError('Sample %s is not among the ones I have to unfold!' % sample)
      binning = eval(binning)
      prev = eval(prev)
      startbin = float(startbin)
      info = [i for i in vars_to_unfold if i.var == varname][0]
      #use ONLY the variable we want
      vars_to_unfold = []
      for stop in binning:
         clone = info.clone()
         clone.binning.reco = prev + [startbin, stop]
         clone.dir_postfix = '_'.join(['%.1f' % i for i in clone.binning.reco])
         vars_to_unfold.append( clone )

   if opts.binning and len(opts.binning):
      varname, binning = tuple(opts.binning.split(':'))
      if not any(i.var == varname for i in vars_to_unfold):
         raise ValueError('Sample %s is not among the ones I have to unfold!' % sample)
      binning = eval(binning)
      info = [i for i in vars_to_unfold if i.var == varname][0]
      info.binning.reco = binning
      info.dir_postfix = '_'.join(['%.1f' % i for i in info.binning.reco])
      #use ONLY the variable we want
      vars_to_unfold = [
         info
         ]

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
      

   ##################
   #     PLOTS
   ##################
   if not opts.noplots:
      to_plot = [
         ('all_lep_pt' , 4, 'p_{T}(l)'),
         ("all_ttm"    , 4, 'm(t#bar{t})'),
         ("all_tty"    , 4, 'y(t#bar{t})'),
         ("all_ttpt"   , 4, 'p_{T}(t#bar{t})'),
         ("all_costhetastar", 4, ''),
         ("all_njet", 4, '# of jets'),
         ('all_%s' % discriminant, 2, 'discriminant'),
         ]

      plotter.plot_mc_shapes(
         '', 'all_%s' % discriminant, rebin=2, xaxis=discriminant,
         leftside=False, normalize=True, show_err=True, xrange=(-8,1))
      plotter.save('%s_full_shape' % (discriminant), pdf=False)

      plotter.plot_mc_shapes(
         '', 'all_%s' % discriminant, rebin=2, xaxis=discriminant,
         leftside=False, normalize=True, show_err=True, xrange=(-8,1),
         use_only=set(['ttJets_rightTlep', 'ttJets_wrongAssign', 'ttJets_other']),
         ratio_range=1)
      plotter.save('%s_bkg_shape' % (discriminant), pdf=False)

      for var, rebin, xaxis in to_plot:
         plotter.plot_mc_vs_data('', var, leftside=False, rebin=rebin, xaxis=xaxis)
         plotter.save(var, pdf=False)

      for info in vars_to_unfold:
         var = info.var
         plotter.plot_mc_vs_data(
            '', 'all_%s' % var, 
            leftside=False, rebin=info.binning.reco, xaxis=info.xtitle)
         plotter.save(var, pdf=False)
      
         previous = info.binning.reco[0]
         plotter.set_subdir('%s/slices' % var)
         for idx, vbin in enumerate(info.binning.reco[1:]):
            plotter.plot_mc_vs_data(
               '', 'all_%s_%s' % (discriminant, var), leftside=False, 
               rebin = full_discr_binning[0],
               xaxis=discriminant,
               preprocess=lambda x: urviews.ProjectionView(x, 'X', [previous, vbin])
               )
            plotter.save('%s_slice_%i' % (discriminant, idx), pdf=False)
      
            plotter.plot_mc_shapes(
               '', 'all_%s_%s' % (discriminant, var), leftside=False, 
               rebin = full_discr_binning[0],
               xaxis=discriminant, normalize=True, show_err=True, xrange=(-8,1),
               preprocess=lambda x: urviews.ProjectionView(x, 'X', [previous, vbin]))
            plotter.save('%s_slice_%i_shape' % (discriminant, idx), pdf=False)
         
            plotter.plot_mc_shapes(
               '', 'all_%s_%s' % (discriminant, var), leftside=False, 
               rebin = full_discr_binning[0],
               xaxis=discriminant, normalize=True, show_err=True, xrange=(-8,1),
               preprocess=lambda x: urviews.ProjectionView(x, 'X', [previous, vbin]),
               use_only=set(['ttJets_rightTlep', 'ttJets_wrongAssign', 'ttJets_other']),
               ratio_range=1)
            plotter.save('%s_bkgslice_%i_shape' % (discriminant, idx), pdf=False)
            previous = vbin



   ##################
   #     CARDS
   ##################
   if not opts.noshapes:

      for info in vars_to_unfold:
         var = info.var
         plotter.set_subdir(
            os.path.join(
               var,
               opts.subdir,
               info.dir_postfix if hasattr(info, 'dir_postfix') else ''
               )
            )
         plotter.write_shapes(
            '', var, 'all_%s_%s' % (discriminant, var), 
            var_binning=info.binning.reco,
            disc_binning = lambda x, *args: full_discr_binning[0]
            ) 
         plotter.add_systematics()
         #plotter.card.add_systematic('lumi', 'lnN', '.*', '[^t]+.*', 1.05)
         plotter.save_card(var)
      
      #Migration matrices
      plotter.set_subdir('')
      fname = os.path.join(plotter.outputdir, 'migration_matrices.root')

      with io.root_open(fname, 'recreate') as mfile:
         for info in vars_to_unfold:
            var = info.var
            dirname = var
            if hasattr(info, 'dir_postfix'):
               dirname += '_%s' % info.dir_postfix
            mfile.mkdir(dirname).cd()
            matrix_path = 'RECO/truth_%s_matrix_%s' % (var, phase_space)
            tt_view = plotter.get_view(plotter.ttbar_to_use, 'unweighted_view')
            matrix_view_unscaled = plotter.rebin_view(
               tt_view, 
               [info.binning.gen, info.binning.reco]
               )
            mig_matrix_unscaled = matrix_view_unscaled.Get(matrix_path)
            mig_matrix_unscaled.SetName('migration_matrix')
            mig_matrix_unscaled.Write()

            tt_view = plotter.get_view(plotter.ttbar_to_use)
            matrix_view = plotter.rebin_view(
               tt_view,
               [info.binning.gen, info.binning.reco]
               )
            mig_matrix = matrix_view.Get(matrix_path)
            mig_matrix.SetName('migration_matrix_scaled')
            mig_matrix.Write()

            reco_distro = mig_matrix.ProjectionY() 
            reco_distro.SetName('reco_distribution')
            reco_distro.Write()
            
            prefit_view = plotter.rebin_view(
               plotter.get_view('ttJets_rightAssign'), 
               info.binning.reco
               )
            plotting.views.SumView.debug = True
            prefit_plot = prefit_view.Get('all_%s' % var)
            prefit_plot.name = 'prefit_distribution'
            prefit_plot.Write()

            thruth_distro = plotter.rebin_view(
               tt_view,
               info.binning.gen
               ).Get('RECO/truth_%s_gen_%s' % (var, phase_space))
            thruth_distro.SetName('true_distribution')
            thruth_distro.Write()         

      logging.info('file: %s written' % fname)
      


if __name__ == '__main__':
   parser = ArgumentParser()
   parser.add_argument('--noplots', dest='noplots', action='store_true',
                       help='skip plot making')
   parser.add_argument('--noshapes', dest='noshapes', action='store_true',
                       help='skip shape making')
   parser.add_argument('--optimize_binning', type=str,
                       help='map observable:[] first bin range for optimization'
                       )
   parser.add_argument('--binning', type=str,
                       help='map observable:[binning] binning for optimization'
                       )
   parser.add_argument('--subdir', type=str, default='',
                       help='sub directory to store shapes'
                       )
   opts = parser.parse_args()
   run_module(**dict(opts._get_kwargs()))
   
