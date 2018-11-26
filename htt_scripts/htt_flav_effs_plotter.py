#! /bin/env python

from URAnalysis.PlotTools.Plotter import Plotter, BasePlotter
from URAnalysis.PlotTools.BasePlotter import LegendDefinition
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
#from fnmatch import fnmatch
import URAnalysis.Utilities.prettyjson as prettyjson
#import URAnalysis.Utilities.quad as quad
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
from argparse import ArgumentParser
import math
#from uncertainties import ufloat
#from URAnalysis.Utilities.datacard import DataCard
#from URAnalysis.Utilities.tables import latex_table
#from URAnalysis.Utilities.latex import t2latex
#from URAnalysis.Utilities.roottools import Envelope
from URAnalysis.PlotTools.views.RebinView import RebinView
import re
import itertools
#import rootpy.stats as stats
import numpy as np

parser = ArgumentParser()
parser.add_argument('mode', choices=['electrons', 'muons'], help='choose leptonic decay type')
parser.add_argument('--preselection', action='store_true', help='')
parser.add_argument('--plots', action='store_true', help='')
parser.add_argument('--all', action='store_true', help='')
parser.add_argument('--njets', default='', help='choose 3, 4+')
parser.add_argument('--qcd_est', action='store_true', help='QCD estimation will be made using ABCD method')
parser.add_argument('--qcd_yields', action='store_true', help='Makes table comparing QCD yields from different methods, file with values used in ML Fit')
args = parser.parse_args()
arg_options = vars(args)


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

        self.tt_to_use = 'ttJets'
        self.tt_shifted = {
        #    'mtop_up'   : 'ttJets_mtopup',
        #    'mtop_down' : 'ttJets_mtopdown',
        #    'isr_up'   : 'ttJets_isrup',
        #    'isr_down' : 'ttJets_isrdown',
        #    'fsr_up'   : 'ttJets_fsrup',
        #    'fsr_down' : 'ttJets_fsrdown',
            }
        jobid = os.environ['jobid']

        files = glob.glob('results/%s/htt_qcd_est/*.root' % jobid)
        MC_files = [ fname for fname in files if not ('data' in fname ) ]
        data_files = [ fname for fname in files if (os.path.basename(fname).startswith('data') and mode[:-1] in os.path.basename(fname).lower() ) ]

        files = data_files+MC_files
        #logging.debug('files found %s' % files.__repr__())

        lumis = glob.glob('inputs/%s/*.lumi' % jobid)
        data_lumis = [ fname for fname in lumis if (os.path.basename(fname).startswith('data') and mode[:-1] in os.path.basename(fname).lower() ) ]
        MC_lumis = [ fname for fname in lumis if not ('data' in fname) ]
        lumis = data_lumis+MC_lumis
        #logging.debug('lumi files found %s' % lumis.__repr__())

        #self.tt_lhe_weights = prettyjson.loads(
        #    open('inputs/%s/ttJets.weights.json' % jobid).read()
        #    )
        outdir= 'plots/%s/htt_qcd_est/%s' % (jobid, mode)

        super(HTTPlotter, self).__init__(
            files, lumis, outdir, styles, None, lumi
            #defaults = {'save' : {'png' : True, 'pdf' : False}}
            )

        #set_trace()
        #select only mode subdir
        for info in self.views.itervalues():
            if args.njets == '3':
                info['view'] = views.SubdirectoryView(info['view'], mode+'/3Jets')
                info['unweighted_view'] = views.SubdirectoryView(info['unweighted_view'], mode+'/3Jets')
            elif args.njets == '4+':
                info['view'] = views.SubdirectoryView(info['view'], mode+'/4PJets')
                info['unweighted_view'] = views.SubdirectoryView(info['unweighted_view'], mode+'/4PJets')
            else:
                if args.flow and args.combinedflow:
                    info['view'] = views.SubdirectoryView(info['view'], '')
                    info['unweighted_view'] = views.SubdirectoryView(info['unweighted_view'], '')
                else:
                    info['view'] = views.SubdirectoryView(info['view'], mode)
                    info['unweighted_view'] = views.SubdirectoryView(info['unweighted_view'], mode)

        #set_trace()
        self.views['data']['view'] = urviews.BlindView(
            self.views['data']['view'], 
            '\w+/tight/MTHigh/csvPass/(:?(:?m_tt)|(:?.+_ctstar)|(:?.+_ctstar_abs)|(:?cdelta_ld)|(:?hframe_ctheta_d))'
            )

        self.defaults = {
            'blurb' : [13, self.views['data']['intlumi']],
            'save' : {'png' : True, 'pdf' : False}
            }
        self.jobid = jobid

        self.views['ttJets_preselection'] = self.views['ttJets']

        self.views['ttJets_right'] = {
            'view' : self.create_tt_subsample(
                ['right'],
                't#bar{t} matched',
                ROOT.kOrange + 9
                )
            }
        self.views['ttJets_matchable'] = {
            'view' : self.create_tt_subsample(
                ['matchable'], #,  'right_th', 'wrong'],
                't#bar{t} matchable',
                ROOT.kRed - 1
                )
            }
        self.views['ttJets_unmatchable'] = {
            'view' : self.create_tt_subsample(
                ['unmatchable'],#  'right_th', 'wrong'],
                't#bar{t} unmatchable',
                ROOT.kMagenta - 2
                )
            }

        self.views['ttJets_other'] = {
            'view' : self.create_tt_subsample(
                ['noslep'], 
                'Other t#bar{t}',
                ROOT.kCyan - 1
                )
            }

        ##
        ## General plotting views
        ##
        self.views['ttJets_generic'] = {
            'view' : views.TitleView(
                views.StyleView(
                    views.SumView(
                        self.views['ttJets_other'      ]['view'],
                        self.views['ttJets_unmatchable']['view'],
                        self.views['ttJets_matchable'  ]['view'],
                        self.views['ttJets_right'      ]['view'],                       
                        ),
                    fillcolor = 'r'
                    #fillcolor = ROOT.kOrange + 1
                    ),
                't#bar{t}'
                )
            }

        self.views['EWK'] = {
            'view' : views.TitleView(
                views.StyleView(
                    #views.SumView(*[self.get_view(i) for i in ['[WZ][WZ]', '[WZ]Jets', 'tt[WZ]*']]),
                    views.SumView(*[self.get_view(i) for i in ['[WZ][WZ]', 'ZJets', 'W[1-4]Jets', 'tt[WZ]*']]),
                    fillcolor = '#FFD700'
                    #fillcolor = ROOT.kGreen + 1
                    ),
                    'EW'
                )
            }

        self.views['AllButTT'] = {
            'view' : views.TitleView(
                views.StyleView(
                    views.SumView(
                        self.get_view('EWK'),
                        self.get_view('QCD*'),
                        self.get_view('single*'),
                        ),
                    fillcolor = ROOT.kGray
                    ),
                'Other'
                )
            }


        self.generic_mcs = [
            'EWK',
            'single*',          
            'ttJets_generic',
            'QCD*',
            ]

        self.mc_samples = self.generic_mcs
        #self.split_mcs = [
        #    'AllButTT',
        #    'ttJets_other',
        #    'ttJets_unmatchable',
        #    'ttJets_matchable',
        #    'ttJets_right',
        #    ]

        #non_qcd = ['TT', 'VV', 'TTV',   'WJets', 'ZJets',   'tChannel', 'tWChannel', 'sChannel']




    def create_subsample(self, baseview, subdirs, title, color='#9999CC'):
        return views.StyleView(
            views.TitleView(
                views.SumView(
                    *[views.SubdirectoryView(self.views[baseview]['view'], i) for i in subdirs]
                     ),
                title
                ),
            fillcolor = color
            )

    def create_tt_subsample(self, subdirs, title, color='#9999CC'):
        shifts = {shift : self.create_subsample(view, ['%s/nosys' % i for i in subdirs], title, color=color) for shift, view in self.tt_shifted.iteritems()}
        shifts[''] = self.create_subsample('ttJets', subdirs, title, color=color)
        return urviews.MultifileView(**shifts)


    def make_preselection_plot(self, *args, **kwargs):
        systematics = None
        if 'sys_effs' in kwargs:
            systematics = kwargs['sys_effs']
            del kwargs['sys_effs']
        mc_default = self.mc_samples
        self.mc_samples = [
            'EWK',
            'single*',          
            'ttJets_preselection',
            'QCD*',
            ]
        self.plot_mc_vs_data(*args, **kwargs)
        if kwargs['plot_unc']:
            mc_stack = [i for i in self.keep if isinstance(i, ROOT.THStack)][0]
            unc_hist = sum(mc_stack.hists)
            unc_hist.fillcolor = 'black'
            unc_hist.fillstyle = 3013
            unc_hist.title = 'Uncertainty'
            unc_hist.drawstyle = 'pe2'
            unc_hist.markerstyle = 0
            #set_trace()

            if kwargs['nodata']:
                plotter.overlay([mc_stack, unc_hist], legend_def=LegendDefinition(position='NE'),
                    logy=kwargs['logy'], x_range=kwargs['xrange']
                )
            else:
                #set_trace()
                data = [i for i in self.keep if isinstance(i, ROOT.TH1)][0]
                plotter.overlay_and_compare(
                    [mc_stack, unc_hist], data,
                    method='datamc', legend_def=LegendDefinition(position='NE'),
                    logy=kwargs['logy'], x_range=kwargs['xrange']
                )
            #set_trace()

        #if systematics is not None:
        #    data = [i for i in self.keep if isinstance(i, ROOT.TH1)][0]
        #    mc_stack = [i for i in self.keep if isinstance(i, ROOT.THStack)][0]
        #    stack_sum = sum(mc_stack.hists)
        #    set_trace()
        #    self.reset()
        #    dirname = args[0].split('/')[0]
        #    path = args[0]
        #    args = list(args)
        #    args[0] = path.replace(dirname, '%s_up' % systematics)
        #    kwargs['nodata'] = True
        #    self.plot_mc_vs_data(*args, **kwargs)
        #    stack_up = [i for i in self.keep if isinstance(i, ROOT.THStack)][0]
        #    self.reset()
        #    s_up = sum(stack_up.hists)
        #    for ibin, jbin in zip(stack_sum, s_up):
        #        ibin.error = quad.quad(ibin.error, abs(ibin.value - jbin.value))
        #    stack_sum.fillcolor = 'black'
        #    stack_sum.fillstyle = 3013
        #    stack_sum.title = 'uncertainty'
        #    stack_sum.drawstyle = 'pe2'
        #    stack_sum.markerstyle = 0
        #    plotter.overlay_and_compare(
        #        [mc_stack, stack_sum], data,
        #        xtitle = kwargs.get('xaxis',''),
        #        ytitle='Events', ignore_style=True,             
        #        method='ratio'
        #        )
        #    # Add legend
        #    self.pad.cd()
        #    self.add_legend(
        #        [mc_stack, stack_sum, data], kwargs.get('leftside', True), 
        #        entries=len(mc_stack.hists)+2
        #        )

        #self.mc_samples = mc_default




    def flav_fracs_and_effs(self, var, xtitle, name, xlims, nbinsx): ## find fractions and efficiencies of MC for b and prompt quarks

        #set_trace()
        ## find fractions and efficiencies of QCD MC for leps matched to b and prompt quarks
        tight_CSV_b = self.get_view('QCD*').Get('nosys/tight/csvPass/%s_B' % var )
        tight_CSV_b = RebinView.rebin(tight_CSV_b, nbinsx)
        loose_CSV_b = self.get_view('QCD*').Get('nosys/looseNOTTight/csvPass/%s_B' % var )
        loose_CSV_b = RebinView.rebin(loose_CSV_b, nbinsx)
        tight_CSV_p = self.get_view('QCD*').Get('nosys/tight/csvPass/%s_Prompt' % var )
        tight_CSV_p = RebinView.rebin(tight_CSV_p, nbinsx)
        loose_CSV_p = self.get_view('QCD*').Get('nosys/looseNOTTight/csvPass/%s_Prompt' % var )
        loose_CSV_p = RebinView.rebin(loose_CSV_p, nbinsx)
        tight_nonCSV_b = self.get_view('QCD*').Get('nosys/tight/csvFail/%s_B' % var )
        tight_nonCSV_b = RebinView.rebin(tight_nonCSV_b, nbinsx)
        loose_nonCSV_b = self.get_view('QCD*').Get('nosys/looseNOTTight/csvFail/%s_B' % var )
        loose_nonCSV_b = RebinView.rebin(loose_nonCSV_b, nbinsx)
        tight_nonCSV_p = self.get_view('QCD*').Get('nosys/tight/csvFail/%s_Prompt' % var )
        tight_nonCSV_p = RebinView.rebin(tight_nonCSV_p, nbinsx)
        loose_nonCSV_p = self.get_view('QCD*').Get('nosys/looseNOTTight/csvFail/%s_Prompt' % var )
        loose_nonCSV_p = RebinView.rebin(loose_nonCSV_p, nbinsx)

        N_b_CSV = tight_CSV_b+loose_CSV_b
        N_p_CSV = tight_CSV_p+loose_CSV_p
        N_b_nonCSV = tight_nonCSV_b+loose_nonCSV_b
        N_p_nonCSV = tight_nonCSV_p+loose_nonCSV_p

            ## format and plot lep pT hists for number of leps matched to b and prompt
        plotter.plot( N_b_CSV, x_range=xlims, ytitle='N_{b}^{CSV}', legend_def=LegendDefinition(position='NE'), drawstyle='E0 X0', legendstyle='p' )
        plotter.save( 'N_b_CSV_%s' % name )
        plotter.plot( N_p_CSV, x_range=xlims, ytitle='N_{prompt}^{CSV}', legend_def=LegendDefinition(position='NE'), drawstyle='E0 X0', legendstyle='p' )
        plotter.save( 'N_p_CSV_%s' % name )
        plotter.plot( N_b_nonCSV, x_range=xlims, ytitle='N_{b}^{nonCSV}', legend_def=LegendDefinition(position='NE'), drawstyle='E0 X0', legendstyle='p' )
        plotter.save( 'N_b_nonCSV_%s' % name )
        plotter.plot( N_p_nonCSV, x_range=xlims, ytitle='N_{prompt}^{nonCSV}', legend_def=LegendDefinition(position='NE'), drawstyle='E0 X0', legendstyle='p' )
        plotter.save( 'N_p_nonCSV_%s' % name )

        #set_trace()
            ## format and plot lep pT hists for fraction of leps matched to b and prompt
        f_b_CSV = N_b_CSV/(N_b_CSV+N_p_CSV)
        plotter.plot( f_b_CSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='f_{b}^{CSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')
        plotter.save( 'Fractions_b_CSV_%s' % name )
        f_p_CSV = N_p_CSV/(N_b_CSV+N_p_CSV)
        plotter.plot( f_p_CSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='f_{prompt}^{CSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')
        plotter.save( 'Fractions_p_CSV_%s' % name )
        #plotter.set_histo_style( f_p_CSV, name='Prompt', color='r', fillstyle=0, drawstyle='E0 X0')

        f_b_nonCSV = N_b_nonCSV/(N_b_nonCSV+N_p_nonCSV)
        plotter.plot( f_b_nonCSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='f_{b}^{nonCSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')
        plotter.save( 'Fractions_b_nonCSV_%s' % name )
        f_p_nonCSV = N_p_nonCSV/(N_b_nonCSV+N_p_nonCSV)
        plotter.plot( f_p_nonCSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='f_{prompt}^{nonCSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')
        plotter.save( 'Fractions_p_nonCSV_%s' % name )

        #set_trace()
            ## format and plot lep pT hists for efficiency of leps matched to b and prompt
        e_b_CSV = tight_CSV_b/N_b_CSV
        plotter.plot( e_b_CSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='#epsilon_{b}^{CSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')
        plotter.save( 'Efficiencies_b_CSV_%s' % name )
        e_p_CSV = tight_CSV_p/N_p_CSV
        plotter.plot( e_p_CSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='#epsilon_{prompt}^{CSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')
        plotter.save( 'Efficiencies_p_CSV_%s' % name )

        e_b_nonCSV = tight_nonCSV_b/N_b_nonCSV
        plotter.plot( e_b_nonCSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='#epsilon_{b}^{nonCSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')
        plotter.save( 'Efficiencies_b_nonCSV_%s' % name )
        e_p_nonCSV = tight_nonCSV_p/N_p_nonCSV
        plotter.plot( e_p_nonCSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='#epsilon_{prompt}^{nonCSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')
        plotter.save( 'Efficiencies_p_nonCSV_%s' % name )


        #set_trace()


    #def scale_factor(self, var, xtitle, name, xlims, nbinsx): ## find fractions and efficiencies of MC for b and prompt quarks

        data_Iso_nonCSV_hist = self.get_view('data').Get('nosys/tight/csvFail/%s_B' % var )+self.get_view('data').Get('nosys/tight/csvFail/%s_Prompt' % var )
        data_Iso_nonCSV_hist = RebinView.rebin(data_Iso_nonCSV_hist, nbinsx)
        data_nonIso_nonCSV_hist = self.get_view('data').Get('nosys/looseNOTTight/csvFail/%s_B' % var )+self.get_view('data').Get('nosys/looseNOTTight/csvFail/%s_Prompt' % var )
        data_nonIso_nonCSV_hist = RebinView.rebin(data_nonIso_nonCSV_hist, nbinsx)
        data_nonCSV_hist = data_Iso_nonCSV_hist+data_nonIso_nonCSV_hist

        r_data = data_Iso_nonCSV_hist/data_nonCSV_hist
        plotter.plot( r_data, x_range=xlims, xtitle=xtitle, ytitle='r=N_{data}^{Iso-nonCSV}/N_{data}^{nonCSV}', legend_def=LegendDefinition(position='NE'), legendstyle='p', drawstyle='E0 X0' )
        plotter.save( 'r_data_nonCSV_%s' % name )
        #set_trace()

        scale_factor = r_data/(f_b_nonCSV*e_b_nonCSV+f_p_nonCSV*e_p_nonCSV)
        plotter.plot( scale_factor, x_range=xlims, xtitle=xtitle, ytitle='SF^{nonCSV}', drawstyle='E0 X0' )
        plotter.save( 'SF_nonCSV_%s' % name )



plotter = HTTPlotter(args.mode)

preselection = [
    (True, "min_lep_jet_dr" , "min #DeltaR(lep, jets)",  1, (0, 0.5), False),
    (True, "lep_pt_B" , "p_{T}(l matched to b-quark) (GeV)",  20, (0,300), False),
    (True, "lep_pt_Prompt" , "p_{T}(l matched to prompt quark) (GeV)",  20, (0,300), False),
    (True, "lep_eta_B" , "#eta(l matched to b-quark)",  20, (-2.4, 2.4), False),
    (True, "lep_eta_Prompt" , "#eta(l matched to prompt quark)",  20, (-2.4, 2.4), False),
    (False, "MT" , "M_{T}",  10, (0, 300), False)
]

if args.preselection or args.all:
    if args.njets == '3':
        plotter.set_subdir('3Jets/preselection')
    elif args.njets == '4+':
        plotter.set_subdir('4PJets/preselection')
    else:
        raise RuntimeError('Your choice for --njets is invalid!')
        #plotter.set_subdir('Incl/preselection')
    #set_trace()
    for logy, var, axis, rebin, x_range, leftside in preselection:
    #for logy, var, axis, rebin, x_range, leftside in preselection + permutations:
        NoData = False
        if 'B' in var or 'Prompt' in var: NoData = True

        plotter.make_preselection_plot(
            'nosys/preselection', var, sort=True,
            xaxis=axis, leftside=leftside, rebin=rebin, 
            show_ratio=True, ratio_range=0.5, xrange=x_range, logy=True, nodata=NoData, plot_unc=True)        
        plotter.save(var+'_logy')

        plotter.make_preselection_plot(
            'nosys/preselection', var, sort=True,
            xaxis=axis, leftside=leftside, rebin=rebin, 
            show_ratio=True, ratio_range=0.5, xrange=x_range, logy=False, nodata=NoData, plot_unc=True)        
        plotter.save(var)

if args.qcd_yields or args.all:

    plotter.set_subdir( '%s/ISO_BTag_Est' % '3Jets' if args.njets == '3' else '4PJets' )
        
    lepIso_bTag = [
        ('lep_pt', 'p_{T}(l) (GeV)', 'Pt', (0, 300), 20),
        ('lep_eta', '#eta(l)', 'Eta', (-2.4, 2.4), 20)
    ]
    for var, overlay_xtitle, name_var, xlims, nbinsx in lepIso_bTag:
        #plotter.scale_factor(var, overlay_xtitle, name_var, xlims, nbinsx)

        plotter.flav_fracs_and_effs(var, overlay_xtitle, name_var, xlims, nbinsx)
    #plotter.set_subdir(plotter.base_out_dir)
    #yields_dict = {}
    #qcd_MC_scale, qcd_MC_error = plotter.QCD_est_from_MC('njets')
    #qcd_abcd_scale, qcd_abcd_error = plotter.QCD_est_from_abcd('njets')
    #qcd_mlFit_scale, qcd_mlFit_error = plotter.QCD_est_from_mlFit()

    ##set_trace()
    #rows = [
    #    ("Method", "Est. Value", "+", "-"),
    #    ("MC simulation", format(qcd_MC_scale, '.0f'), format(qcd_MC_error, '.0f'), format(qcd_MC_error, '.0f')),
    #    ("ABCD", format(qcd_abcd_scale, '.0f'), format(qcd_abcd_error, '.0f'), format(qcd_abcd_error, '.0f')),
    #    ("ML Fit", format(qcd_mlFit_scale, '.0f'), format(qcd_mlFit_error[0], '.0f'), format(qcd_mlFit_error[1], '.0f'))
    #]

    #plotter.print_table(rows, filename='%s/%s/QCD_Est_Results.raw_txt' % (plotter.outputdir, '3Jets' if args.njets == '3' else '4PJets'), print_output=True )
    #print '\n-----   Table comparing estimated QCD yields written to %s/%s/QCD_Est_Results.raw_txt   -----\n' % (plotter.outputdir, '3Jets' if args.njets == '3' else '4PJets')

    #    ## creates file to be used in estimating background with ml Fit
    #with open('%s/%s/%s_%s_yields.json' % (plotter.outputdir, '3Jets' if args.njets == '3' else '4PJets', plotter.outputdir.split('/')[-1], '3Jets' if args.njets == '3' else '4PJets'), 'w') as f: # write to same dir as plots
    #    f.write(prettyjson.dumps(yields_dict))
    #print '\n-----   File with values to be used in ML Fit for QCD yields written to %s/%s/%s_%s_yields.json   -----\n' %\
    #             (plotter.outputdir, '3Jets' if args.njets == '3' else '4PJets',\
    #              plotter.outputdir.split('/')[-1], '3Jets' if args.njets == '3' else '4PJets')



if args.plots or args.all:

    #qcd_est_scale, qcd_est_error = plotter.QCD_est_from_mlFit() if args.qcd_est else 1., (1.,1.)
    qcd_est_scale, qcd_est_error = 1., (1.,1.)

    #plotter.compare_QCD_estimations()
    #set_trace()
    for dirid in itertools.product(['looseNOTTight', 'tight'], ['csvPass', 'csvFail']):
        tdir = '%s/%s' % dirid
        if args.njets == '3':
            plotter.set_subdir('3Jets/'+tdir)
        elif args.njets == '4+':
            plotter.set_subdir('4PJets/'+tdir)
        else:
            raise RuntimeError('Your choice for --njets is invalid!')
        #plotter.set_subdir(tdir)
        #set_trace()
        #first = True

        #for var, xaxis, yaxis, xbins, ybins in vars2D:
        #    for i in plotter.mc_views(1, None, 'nosys/%s/csvPass' % tdir, False):
        #        hist = i.Get(var)
        #        hist = RebinView.newRebin2D(hist, xbins, ybins)
        #        #set_trace()
        #        plotter.set_histo_style(hist, xtitle=xaxis, ytitle=yaxis, drawstyle='colz')
        #        plotter.plot(hist)
        #        plotter.save('%s_%s' % (i.Get(var).title, var) )

            #set_trace()
        for logy, var, axis, rebin, x_range, leftside in preselection:
            NoData = False
            if 'B' in var or 'Prompt' in var: NoData = True

            plot_unc = True
            #set_trace()
            plotter.plot_mc_vs_data(
                'nosys/%s' % tdir, var, sort=False,
                #'nosys/%s' % tdir, var, sort=True,
                xaxis=axis, leftside=leftside, rebin=rebin,
                show_ratio=True, ratio_range=0.2, xrange=x_range,
                logy=False, qcd_renorm=False, qcd_scale=qcd_est_scale, nodata=NoData)

            #set_trace()
            if plot_unc:
                mc_stack = [i for i in plotter.keep if isinstance(i, ROOT.THStack)][0]
                unc_hist = sum(mc_stack.hists)
                unc_hist.fillcolor = 'black'
                unc_hist.fillstyle = 3013
                unc_hist.title = 'Uncertainty'
                unc_hist.drawstyle = 'pe2'
                unc_hist.markerstyle = 0
                #set_trace()

                if NoData:
                    plotter.overlay([mc_stack, unc_hist], legend_def=LegendDefinition(position='NE'),
                        logy=False, x_range=x_range
                    )
                else:
                    #set_trace()
                    data = [i for i in plotter.keep if isinstance(i, ROOT.TH1)][0]
                    plotter.overlay_and_compare(
                        [mc_stack, unc_hist], data,
                        method='datamc', legend_def=LegendDefinition(position='NE'),
                        logy=False, x_range=x_range
                    )

            #set_trace()
            print var
            plotter.save(var)

            plotter.plot_mc_vs_data(
                'nosys/%s' % tdir, var, sort=False,
                #'nosys/%s' % tdir, var, sort=True,
                xaxis=axis, leftside=leftside, rebin=rebin,
                show_ratio=True, ratio_range=0.2, xrange=x_range,
                logy=True, qcd_renorm=False, qcd_scale=qcd_est_scale, nodata=NoData)

            if plot_unc:
                mc_stack = [i for i in plotter.keep if isinstance(i, ROOT.THStack)][0]
                unc_hist = sum(mc_stack.hists)
                unc_hist.fillcolor = 'black'
                unc_hist.fillstyle = 3013
                unc_hist.title = 'Uncertainty'
                unc_hist.drawstyle = 'pe2'
                unc_hist.markerstyle = 0
                #set_trace()

                if NoData:
                    plotter.overlay([mc_stack, unc_hist], legend_def=LegendDefinition(position='NE'),
                        logy=True, x_range=x_range
                    )
                else:
                    #set_trace()
                    data = [i for i in plotter.keep if isinstance(i, ROOT.TH1)][0]
                    plotter.overlay_and_compare(
                        [mc_stack, unc_hist], data,
                        method='datamc', legend_def=LegendDefinition(position='NE'),
                        logy=True, x_range=x_range
                    )
            print var
            #set_trace()
            plotter.save(var+'_logy')
