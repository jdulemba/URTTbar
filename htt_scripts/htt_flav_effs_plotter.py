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
from rootpy.plotting import Hist
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
from rootpy import asrootpy
#from uncertainties import ufloat
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

        outdir= 'plots/%s/htt_qcd_est/%s' % (jobid, mode)

        super(HTTPlotter, self).__init__(
            files, lumis, outdir, styles, None, lumi
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


        self.fitopts_dict = {
            'sigmoid' : {
                            'fcn_string' : "p0/(1. + TMath::Exp( -p1*(x-p2) ) )",
                            'npars'      : 3
            },
            'erf'     : {
                            'fcn_string' : "p0*(TMath::Erf( p1*(x-p2)) ) + p3",
                            'npars'      : 4
            },
            'constant': {
                            'fcn_string' : "p0",
                            'npars'      : 1
            },
            'line'    : {
                            'fcn_string' : "p0 + p1*x",
                            'npars'      : 2
            },
            'parabola': {
                            'fcn_string' : "p0 + p1*x + p2*x*x",
                            'npars'      : 3
            },
        }


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



        ### checks overflow bins for hpass and htotal so errors with TEff can be avoided
    def check_overflow_bins(self, hpass, htotal):

        hpass_xoverflow = hpass.GetBinContent(hpass.GetNbinsX()+1)
        htotal_xoverflow = htotal.GetBinContent(hpass.GetNbinsX()+1)

        if hpass_xoverflow > htotal_xoverflow:
            hpass.SetBinContent(hpass.GetNbinsX()+1, 0.0)
            htotal.SetBinContent(htotal.GetNbinsX()+1, 0.0)

        return hpass, htotal


        ### extract pars from fit
    def print_fitpars(self, fitresult, numpars):

        fitpars = '#chi^{2}/ndof=%.3f/%i' % (fitresult.GetChisquare(), fitresult.GetNDF())
        for i in range(numpars):
            fitpars += '\np_{%i}=%.2f#pm%.2f' % (i, fitresult.GetParameter(i), fitresult.GetParError(i))

        return fitpars


    def teff2hist(self, teff, **kwargs):
        
        hpass = teff.GetCopyPassedHisto()

        mean_val = [teff.GetEfficiency(i) for i in range(1, hpass.GetNbinsX()+1) ]
        err_up = [teff.GetEfficiencyErrorUp(i) for i in range(1, hpass.GetNbinsX()+1) ]
        err_dw = [teff.GetEfficiencyErrorLow(i) for i in range(1, hpass.GetNbinsX()+1) ]

        ## make errors symmetric for Hist object (take max val between up/down)
        err = [max(err_up[i], err_dw[i]) for i in range(len(err_up))]

        xbins = sorted(set([hpass.GetXaxis().GetBinLowEdge(i) for i in range(1, hpass.GetNbinsX()+1)]+[hpass.GetXaxis().GetBinUpEdge(i) for i in range(1, hpass.GetNbinsX()+1)]))


        fig_name = kwargs['name']+'_fit' if 'name' in kwargs else 'test_fit'
        hist = Hist(xbins, name=fig_name, title=fig_name)
        for binx in range(1, len(xbins)):
            hist[binx].value = mean_val[binx-1]
            hist[binx].error = err[binx-1]

        hist = asrootpy(hist).Clone()

        xtitle = kwargs['x_title'] if 'x_title' in kwargs else ''
        ytitle = kwargs['y_title'] if 'y_title' in kwargs else ''

        plotter.set_histo_style(hist, xtitle=xtitle, ytitle=ytitle, drawstyle='E0 X0')
        plotter.plot(hist)

        if 'fit' in kwargs:

            with io.root_open('muons_3jets_fits.root', 'update') as out:
                hist.Write()

            if kwargs['fit'] == 'compare':
                fit = ROOT.TF1("const", "pol0", hpass.GetXaxis().GetXmin(), hpass.GetXaxis().GetXmax())
                fit.SetLineColor(2)
                fit.SetLineStyle(4)
                hist.Fit(fit)

                fit2 = ROOT.TF1("line", "pol1", hpass.GetXaxis().GetXmin(), hpass.GetXaxis().GetXmax())
                fit2.SetLineColor(3)
                fit2.SetLineStyle(4)
                hist.Fit(fit2, "+")

                fit3 = ROOT.TF1("line", "pol2", hpass.GetXaxis().GetXmin(), hpass.GetXaxis().GetXmax())
                fit3.SetLineColor(4)
                fit3.SetLineStyle(4)
                hist.Fit(fit3, "+")

            else:
                #set_trace()
                fit_model = plotter.fitopts_dict[kwargs['fit']]['fcn_string']
                num_pars = plotter.fitopts_dict[kwargs['fit']]['npars']
                par_string = ""
                fit_opts = "RMENS"
                fixpars = (False,)*num_pars # no parameters are fixed by default
                if 'Pt' in fig_name:
                    if 'name' in kwargs and 'QCD' in kwargs['name']:
                        par_string = "p0[0.5, 0, 1], p1[1.0], p2[1.0]"
                    #elif 'name' in kwargs and 'data' in kwargs['name']:
                    #    par_string = "p0[1.0, 0, 1], p1[1.0], p2[1.0]"
                    #elif 'name' in kwargs and 'Efficiencies_b' in kwargs['name']:
                    #    par_string = "p0[1.0, 0, 1], p1[1.0], p2[1.0]"
                    #elif 'name' in kwargs and 'Efficiencies_p' in kwargs['name']:
                    #    par_string = "p0[0.5, 0, 1], p1[0.0], p2[10.0]"
                    else:
                        par_string = "p0[1.0, 0, 1], p1[1.0], p2[1.0]"
                else:
                    #if 'name' in kwargs and 'QCD' in kwargs['name']:
                    #    par_string = "p0[0.5], p1[1.0], p2[1.0]"
                    ##elif 'name' in kwargs and 'data' in kwargs['name']:
                    ##    par_string = "p0[1.0, 0, 1], p1[1.0], p2[1.0]"
                    ##elif 'name' in kwargs and 'Efficiencies_b' in kwargs['name']:
                    ##    par_string = "p0[1.0, 0, 1], p1[1.0], p2[1.0]"
                    ##elif 'name' in kwargs and 'Efficiencies_p' in kwargs['name']:
                    ##    par_string = "p0[0.5, 0, 1], p1[0.0], p2[10.0]"
                    #else:
                    #    par_string = "p0[1.0], p1[1.0], p2[1.0]"
                    #par_string = "p0[1.0]"
                    #par_string = "p0[1.0], p1[1.0]"
                    par_string = "p0[1.0], p1[1.0], p2[1.0]"

                fit = plotter.fit_shape( hist, model=(fit_model, par_string),\
                    x_range=(hpass.GetXaxis().GetXmin(), hpass.GetXaxis().GetXmax()),\
                    fix_pars=fixpars, fit_opt=fit_opts )
                #set_trace()

                    ## extract pars and put them in text box
                fit_pars = plotter.print_fitpars(fit, num_pars)
                
                box1 = plotter.make_text_box( fit_pars, position='SE')
                box1.Draw()


        plotter.save(fig_name)


    def teff_comparisons(self, hpass, htotal, xtitle=None, ytitle=None, fig_name=None, **kwargs):

        lwidth = 1

        #set_trace()
        hpass, htotal = plotter.check_overflow_bins(hpass, htotal)
        ratio = ROOT.TEfficiency(hpass, htotal)

        ratio_norm = ratio.Clone()
        ratio_norm.SetStatisticOption(1)
        ratio_norm.SetLineColor(1)
        ratio_norm.SetLineWidth(lwidth)
        ratio_norm.SetMarkerColor(1)
        ratio_norm.SetTitle('Normal;%s;%s' % (xtitle, ytitle) )

        ratio_jeff = ratio.Clone()
        ratio_jeff.SetStatisticOption(5)
        ratio_jeff.SetLineColor(2)
        ratio_jeff.SetLineWidth(lwidth)
        ratio_jeff.SetMarkerColor(2)
        ratio_jeff.SetTitle('Jeffreys;%s;%s' % (xtitle, ytitle) )

        ratio_uni = ratio.Clone()
        ratio_uni.SetStatisticOption(6)
        ratio_uni.SetLineColor(3)
        ratio_uni.SetLineWidth(lwidth)
        ratio_uni.SetMarkerColor(3)
        ratio_uni.SetTitle('Uniform;%s;%s' % (xtitle, ytitle) )

        leg = ROOT.TLegend(0.895,0.14,0.965,0.25)
        leg.AddEntry(ratio_norm, "Normal")
        leg.AddEntry(ratio_jeff, "Jeffreys")
        leg.AddEntry(ratio_uni, "Uniform")

        ratio_uni.Draw()
        ratio_norm.Draw('same')
        ratio_jeff.Draw('same')

        leg.Draw('same')
        plotter.save( fig_name )

        #set_trace()
            ## different types of fitting
        if 'fit' in kwargs:
            plotter.teff2hist(ratio_jeff, x_title=xtitle, y_title=ytitle, name=fig_name, fit=kwargs['fit'])
            ##

        #set_trace()
        if 'intfrac' in kwargs:
            eff, err_up, err_low = ratio_jeff.GetEfficiency(1), ratio_jeff.GetEfficiencyErrorUp(1), ratio_jeff.GetEfficiencyErrorLow(1)
            return eff, err_up, err_low



    def sf_comparisons(self, hpass_data, htotal_data, hpass_QCD, htotal_QCD, xtitle=None, ytitle=None, fig_name=None, fit=False):

        lwidth = 1

        #set_trace()
        hpass_data, htotal_data = plotter.check_overflow_bins(hpass_data, htotal_data)
        r_data = ROOT.TEfficiency(hpass_data, htotal_data)

            ## use Jeffreys stats option (5)
        r_data_jeff = r_data.Clone()
        r_data_jeff.SetStatisticOption(5)
        r_data_jeff.SetLineColor(2)
        r_data_jeff.SetLineWidth(lwidth)
        r_data_jeff.SetMarkerColor(2)
        r_data_jeff.SetTitle('Jeffreys;%s;%s' % (xtitle, ytitle) )

        hpass_QCD, htotal_QCD = plotter.check_overflow_bins(hpass_QCD, htotal_QCD)
        r_QCD = ROOT.TEfficiency(hpass_QCD, htotal_QCD)

            ## use Jeffreys stats option (5)
        r_QCD_jeff = r_QCD.Clone()
        r_QCD_jeff.SetStatisticOption(5)
        r_QCD_jeff.SetLineColor(2)
        r_QCD_jeff.SetLineWidth(lwidth)
        r_QCD_jeff.SetMarkerColor(2)
        r_QCD_jeff.SetTitle('Jeffreys;%s;%s' % (xtitle, ytitle) )

        s_hat = [r_data_jeff.GetEfficiency(i)/r_QCD_jeff.GetEfficiency(i) for i in range(1, hpass_data.GetNbinsX()+1) ]
        s_err_up = [(r_data_jeff.GetEfficiencyErrorUp(i)**2+r_QCD_jeff.GetEfficiencyErrorUp(i)**2)**0.5 for i in range(1, hpass_data.GetNbinsX()+1) ]
        s_err_down = [(r_data_jeff.GetEfficiencyErrorLow(i)**2+r_QCD_jeff.GetEfficiencyErrorLow(i)**2)**0.5 for i in range(1, hpass_data.GetNbinsX()+1) ]

        xbins = sorted(set([hpass_data.GetXaxis().GetBinLowEdge(i) for i in range(1, hpass_data.GetNbinsX()+1)]+[hpass_data.GetXaxis().GetBinUpEdge(i) for i in range(1, hpass_data.GetNbinsX()+1)]))
        s_hist = Hist(xbins, name='SF', title='SF')
        for binx in range(len(xbins)):
            s_hist[binix].value = s_hat[binx]

 
        set_trace()

        r_data_jeff.Draw()
        r_QCD_jeff.Draw('same')        

        leg = ROOT.TLegend(0.895,0.14,0.965,0.25)
        leg.AddEntry(r_data_jeff, "Data")
        leg.AddEntry(r_QCD_jeff, "QCD")


    def sf_simplifications(self, var, xtitle, name, xlims, nbinsx): ## find fractions and efficiencies of MC for b and prompt quarks

        ### combine CSV, nonCSV regions

        plotter.set_subdir( '%s/ISO_BTag_Est_Simps/%s' % ('3Jets' if args.njets == '3' else '4PJets', name) )


        N_Iso_CSV_b = self.get_view('QCD*').Get('nosys/tight/csvPass/%s_B' % var )
        N_Iso_CSV_b = RebinView.rebin(N_Iso_CSV_b, nbinsx)
        N_nonIso_CSV_b = self.get_view('QCD*').Get('nosys/looseNOTTight/csvPass/%s_B' % var )
        N_nonIso_CSV_b = RebinView.rebin(N_nonIso_CSV_b, nbinsx)

        N_Iso_nonCSV_b = self.get_view('QCD*').Get('nosys/tight/csvFail/%s_B' % var )
        N_Iso_nonCSV_b = RebinView.rebin(N_Iso_nonCSV_b, nbinsx)
        N_nonIso_nonCSV_b = self.get_view('QCD*').Get('nosys/looseNOTTight/csvFail/%s_B' % var )
        N_nonIso_nonCSV_b = RebinView.rebin(N_nonIso_nonCSV_b, nbinsx)

        N_Iso_CSV_p = self.get_view('QCD*').Get('nosys/tight/csvPass/%s_Prompt' % var )
        N_Iso_CSV_p = RebinView.rebin(N_Iso_CSV_p, nbinsx)
        N_nonIso_CSV_p = self.get_view('QCD*').Get('nosys/looseNOTTight/csvPass/%s_Prompt' % var )
        N_nonIso_CSV_p = RebinView.rebin(N_nonIso_CSV_p, nbinsx)

        N_Iso_nonCSV_p = self.get_view('QCD*').Get('nosys/tight/csvFail/%s_Prompt' % var )
        N_Iso_nonCSV_p = RebinView.rebin(N_Iso_nonCSV_p, nbinsx)
        N_nonIso_nonCSV_p = self.get_view('QCD*').Get('nosys/looseNOTTight/csvFail/%s_Prompt' % var )
        N_nonIso_nonCSV_p = RebinView.rebin(N_nonIso_nonCSV_p, nbinsx)

        N_b_CSV = N_Iso_CSV_b+N_nonIso_CSV_b
        N_p_CSV = N_Iso_CSV_p+N_nonIso_CSV_p
        N_b_nonCSV = N_Iso_nonCSV_b+N_nonIso_nonCSV_b
        N_p_nonCSV = N_Iso_nonCSV_p+N_nonIso_nonCSV_p

        #set_trace()
            ## format and plot lep pT hists for number of leps matched to b and prompt
        plotter.plot( N_b_CSV, x_range=xlims, ytitle='N_{b}^{CSV}', legend_def=LegendDefinition(position='NE'), drawstyle='E0 X0', legendstyle='p' )
        plotter.save( 'N_b_CSV_%s' % name )
        plotter.plot( N_p_CSV, x_range=xlims, ytitle='N_{prompt}^{CSV}', legend_def=LegendDefinition(position='NE'), drawstyle='E0 X0', legendstyle='p' )
        plotter.save( 'N_p_CSV_%s' % name )
        plotter.plot( N_b_nonCSV, x_range=xlims, ytitle='N_{b}^{nonCSV}', legend_def=LegendDefinition(position='NE'), drawstyle='E0 X0', legendstyle='p' )
        plotter.save( 'N_b_nonCSV_%s' % name )
        plotter.plot( N_p_nonCSV, x_range=xlims, ytitle='N_{prompt}^{nonCSV}', legend_def=LegendDefinition(position='NE'), drawstyle='E0 X0', legendstyle='p' )
        plotter.save( 'N_p_nonCSV_%s' % name )

        #if 'Eta' in name: set_trace()
        #set_trace()

            ## format and plot hists for fractions of leps matched to b and prompt
        single_bin = [nbinsx[0], nbinsx[-1]] if 'Pt' in name else [-2.4, 2.4]

        N_Iso_CSV_b_fracs = RebinView.rebin(N_Iso_CSV_b, single_bin)
        N_nonIso_CSV_b_fracs = RebinView.rebin(N_nonIso_CSV_b, single_bin)
        N_Iso_nonCSV_b_fracs = RebinView.rebin(N_Iso_nonCSV_b, single_bin)
        N_nonIso_nonCSV_b_fracs = RebinView.rebin(N_nonIso_nonCSV_b, single_bin)

        N_Iso_CSV_p_fracs = RebinView.rebin(N_Iso_CSV_p, single_bin)
        N_nonIso_CSV_p_fracs = RebinView.rebin(N_nonIso_CSV_p, single_bin)
        N_Iso_nonCSV_p_fracs = RebinView.rebin(N_Iso_nonCSV_p, single_bin)
        N_nonIso_nonCSV_p_fracs = RebinView.rebin(N_nonIso_nonCSV_p, single_bin)

        N_CSV_b_fracs = N_Iso_CSV_b_fracs + N_nonIso_CSV_b_fracs
        N_nonCSV_b_fracs = N_Iso_nonCSV_b_fracs + N_nonIso_nonCSV_b_fracs
        N_CSV_p_fracs = N_Iso_CSV_p_fracs + N_nonIso_CSV_p_fracs
        N_nonCSV_p_fracs = N_Iso_nonCSV_p_fracs + N_nonIso_nonCSV_p_fracs

        f_b_CSV_val, f_b_CSV_errup, f_b_CSV_errlow = plotter.teff_comparisons(N_CSV_b_fracs, N_CSV_b_fracs+N_CSV_p_fracs,\
                    xtitle=xtitle, ytitle='f_{b}^{CSV}', fig_name='Fractions_b_CSV_%s' % name, intfrac=True )

        f_b_nonCSV_val, f_b_nonCSV_errup, f_b_nonCSV_errlow = plotter.teff_comparisons(N_nonCSV_b_fracs, N_nonCSV_b_fracs+N_nonCSV_p_fracs,\
                    xtitle=xtitle, ytitle='f_{b}^{nonCSV}', fig_name='Fractions_b_nonCSV_%s' % name, intfrac=True )

        plotter.teff_comparisons(N_CSV_p_fracs, N_CSV_b_fracs+N_CSV_p_fracs,\
                    xtitle=xtitle, ytitle='f_{prompt}^{CSV}', fig_name='Fractions_p_CSV_%s' % name, intfrac=True )

        plotter.teff_comparisons(N_nonCSV_p_fracs, N_nonCSV_b_fracs+N_nonCSV_p_fracs,\
                    xtitle=xtitle, ytitle='f_{prompt}^{nonCSV}', fig_name='Fractions_p_nonCSV_%s' % name, intfrac=True )

            ## format and plot hists for efficiencies of leps matched to b and prompt
        N_Iso_b = N_Iso_CSV_b + N_Iso_nonCSV_b
        N_nonIso_b = N_nonIso_CSV_b + N_nonIso_nonCSV_b
        N_Iso_p = N_Iso_CSV_p + N_Iso_nonCSV_p
        N_nonIso_p = N_nonIso_CSV_p + N_nonIso_nonCSV_p

        plotter.teff_comparisons(N_Iso_b, N_Iso_b + N_nonIso_b, xtitle=xtitle,\
            ytitle='#epsilon_{b}=(N_{b}^{Iso-CSV}+N_{b}^{Iso-nonCSV})/(N_{b}^{CSV}+N_{b}^{nonCSV})', fig_name='Efficiencies_b_%s' % name, fit='sigmoid' if 'Pt' in name else 'parabola' )

        plotter.teff_comparisons(N_Iso_p, N_Iso_p + N_nonIso_p, xtitle=xtitle,\
            ytitle='#epsilon_{p}=(N_{p}^{Iso-CSV}+N_{p}^{Iso-nonCSV})/(N_{p}^{CSV}+N_{p}^{nonCSV})', fig_name='Efficiencies_p_%s' % name, fit='sigmoid' if 'Pt' in name else 'parabola' )

        plotter.teff_comparisons(f_b_nonCSV_val*N_Iso_b*(N_Iso_p+N_nonIso_p)+(1-f_b_nonCSV_val)*N_Iso_p*(N_Iso_b+N_nonIso_b),\
            (N_Iso_b+N_nonIso_b)*(N_Iso_p+N_nonIso_p), xtitle=xtitle,\
            ytitle='r_{QCD}^{nonCSV}=f_{b}^{nonCSV}#epsilon_{b}+(1-f_{b}^{nonCSV})#epsilon_{prompt}', fig_name='r_QCD_nonCSV_%s' % name, fit='sigmoid' if 'Pt' in name else 'parabola' )


        #set_trace()

            ### data hists
        data_Iso_nonCSV_hist = self.get_view('data').Get('nosys/tight/csvFail/%s_B' % var )+self.get_view('data').Get('nosys/tight/csvFail/%s_Prompt' % var )
        data_Iso_nonCSV_hist = RebinView.rebin(data_Iso_nonCSV_hist, nbinsx)
        data_nonIso_nonCSV_hist = self.get_view('data').Get('nosys/looseNOTTight/csvFail/%s_B' % var )+self.get_view('data').Get('nosys/looseNOTTight/csvFail/%s_Prompt' % var )
        data_nonIso_nonCSV_hist = RebinView.rebin(data_nonIso_nonCSV_hist, nbinsx)
        data_nonCSV_hist = data_Iso_nonCSV_hist+data_nonIso_nonCSV_hist

        data_Iso_CSV_hist = self.get_view('data').Get('nosys/tight/csvPass/%s_B' % var )+self.get_view('data').Get('nosys/tight/csvPass/%s_Prompt' % var )
        data_Iso_CSV_hist = RebinView.rebin(data_Iso_CSV_hist, nbinsx)
        data_nonIso_CSV_hist = self.get_view('data').Get('nosys/looseNOTTight/csvPass/%s_B' % var )+self.get_view('data').Get('nosys/looseNOTTight/csvPass/%s_Prompt' % var )
        data_nonIso_CSV_hist = RebinView.rebin(data_nonIso_CSV_hist, nbinsx)
        data_CSV_hist = data_Iso_CSV_hist+data_nonIso_CSV_hist

            ## get non-QCD MC hists
        MC_Iso_nonCSV = sum([ i.Get('%s_B' % var) for i in plotter.mc_views(1, None, 'nosys/tight/csvFail', False) if i.Get('%s_B' % var).title != 'QCD' ])+\
                        sum([ i.Get('%s_Prompt' % var) for i in plotter.mc_views(1, None, 'nosys/tight/csvFail', False) if i.Get('%s_Prompt' % var).title != 'QCD' ])
        MC_Iso_nonCSV = RebinView.rebin(MC_Iso_nonCSV, nbinsx)
        MC_nonIso_nonCSV = sum([ i.Get('%s_B' % var) for i in plotter.mc_views(1, None, 'nosys/looseNOTTight/csvFail', False) if i.Get('%s_B' % var).title != 'QCD' ])+\
                        sum([ i.Get('%s_Prompt' % var) for i in plotter.mc_views(1, None, 'nosys/looseNOTTight/csvFail', False) if i.Get('%s_Prompt' % var).title != 'QCD' ])
        MC_nonIso_nonCSV = RebinView.rebin(MC_nonIso_nonCSV, nbinsx)
        MC_nonCSV_hist = MC_Iso_nonCSV+MC_nonIso_nonCSV

        MC_Iso_CSV = sum([ i.Get('%s_B' % var) for i in plotter.mc_views(1, None, 'nosys/tight/csvPass', False) if i.Get('%s_B' % var).title != 'QCD' ])+\
                     sum([ i.Get('%s_Prompt' % var) for i in plotter.mc_views(1, None, 'nosys/tight/csvPass', False) if i.Get('%s_Prompt' % var).title != 'QCD' ])
        MC_Iso_CSV = RebinView.rebin(MC_Iso_CSV, nbinsx)
        MC_nonIso_CSV = sum([ i.Get('%s_B' % var) for i in plotter.mc_views(1, None, 'nosys/looseNOTTight/csvPass', False) if i.Get('%s_B' % var).title != 'QCD' ])+\
                        sum([ i.Get('%s_Prompt' % var) for i in plotter.mc_views(1, None, 'nosys/looseNOTTight/csvPass', False) if i.Get('%s_Prompt' % var).title != 'QCD' ])
        MC_nonIso_CSV = RebinView.rebin(MC_nonIso_CSV, nbinsx)
        MC_CSV_hist = MC_Iso_CSV+MC_nonIso_CSV


        Iso_nonCSV_hist = data_Iso_nonCSV_hist - MC_Iso_nonCSV
        plotter.plot( Iso_nonCSV_hist, xtitle=xtitle, ytitle='N_{data-nonQCD MC}^{Iso-nonCSV}', drawstyle='E0 X0')#, legend_def=LegendDefinition(position='NE'), legendstyle='p', drawstyle='E0 X0' )
        plotter.save( 'N_data_minus_MC_Iso_nonCSV_%s' % name )
        
        Iso_CSV_hist = data_Iso_CSV_hist - MC_Iso_CSV
        plotter.plot( Iso_CSV_hist, xtitle=xtitle, ytitle='N_{data-nonQCD MC}^{Iso-CSV}', drawstyle='E0 X0')#, legend_def=LegendDefinition(position='NE'), legendstyle='p', drawstyle='E0 X0' )
        plotter.save( 'N_data_minus_MC_Iso_CSV_%s' % name )

        nonCSV_hist = data_nonCSV_hist - MC_nonCSV_hist
        plotter.plot( nonCSV_hist, xtitle=xtitle, ytitle='N_{data-nonQCD MC}^{nonCSV}', drawstyle='E0 X0')#, legend_def=LegendDefinition(position='NE'), legendstyle='p', drawstyle='E0 X0' )
        plotter.save( 'N_data_minus_MC_nonCSV_%s' % name )

        CSV_hist = data_CSV_hist - MC_CSV_hist
        plotter.plot( CSV_hist, xtitle=xtitle, ytitle='N_{data-nonQCD MC}^{CSV}', drawstyle='E0 X0')#, legend_def=LegendDefinition(position='NE'), legendstyle='p', drawstyle='E0 X0' )
        plotter.save( 'N_data_minus_MC_CSV_%s' % name )


        plotter.teff_comparisons(Iso_nonCSV_hist, nonCSV_hist, xtitle=xtitle,\
            ytitle='r_{data}^{nonCSV}=N_{data}^{Iso-nonCSV}/N_{data}^{Iso-nonCSV}', fig_name='r_data_nonCSV_%s' % name, fit='sigmoid' if 'Pt' in name else 'parabola' )


        ##plotter.sf_comparisons(Iso_nonCSV_hist, nonCSV_hist, N_Iso_nonCSV_p+N_Iso_nonCSV_b, N_p_nonCSV+N_b_nonCSV, xtitle=xtitle,\
        ##    ytitle='s^{nonCSV}=r_{data-nonQCD MC}^{nonCSV}/r_{QCD}^{nonCSV}', fig_name='SF_nonCSV_%s' % name)
        #set_trace()

        ##scale_factor = r_data/(f_b_nonCSV*e_b_nonCSV+f_p_nonCSV*e_p_nonCSV)
        ##plotter.plot( scale_factor, x_range=xlims, xtitle=xtitle, ytitle='SF^{nonCSV}', drawstyle='E0 X0' )
        ##plotter.save( 'SF_nonCSV_%s' % name )


    def flav_fracs_and_effs(self, var, xtitle, name, xlims, nbinsx): ## find fractions and efficiencies of MC for b and prompt quarks

        plotter.set_subdir( '%s/ISO_BTag_Est/%s' % ('3Jets' if args.njets == '3' else '4PJets', name) )
        
        #set_trace()
        ## find fractions and efficiencies of QCD MC for leps matched to b and prompt quarks
        N_Iso_CSV_b = self.get_view('QCD*').Get('nosys/tight/csvPass/%s_B' % var )
        N_Iso_CSV_b = RebinView.rebin(N_Iso_CSV_b, nbinsx)
        N_nonIso_CSV_b = self.get_view('QCD*').Get('nosys/looseNOTTight/csvPass/%s_B' % var )
        N_nonIso_CSV_b = RebinView.rebin(N_nonIso_CSV_b, nbinsx)
        N_Iso_CSV_p = self.get_view('QCD*').Get('nosys/tight/csvPass/%s_Prompt' % var )
        N_Iso_CSV_p = RebinView.rebin(N_Iso_CSV_p, nbinsx)
        N_nonIso_CSV_p = self.get_view('QCD*').Get('nosys/looseNOTTight/csvPass/%s_Prompt' % var )
        N_nonIso_CSV_p = RebinView.rebin(N_nonIso_CSV_p, nbinsx)
        N_Iso_nonCSV_b = self.get_view('QCD*').Get('nosys/tight/csvFail/%s_B' % var )
        N_Iso_nonCSV_b = RebinView.rebin(N_Iso_nonCSV_b, nbinsx)
        N_nonIso_nonCSV_b = self.get_view('QCD*').Get('nosys/looseNOTTight/csvFail/%s_B' % var )
        N_nonIso_nonCSV_b = RebinView.rebin(N_nonIso_nonCSV_b, nbinsx)
        N_Iso_nonCSV_p = self.get_view('QCD*').Get('nosys/tight/csvFail/%s_Prompt' % var )
        N_Iso_nonCSV_p = RebinView.rebin(N_Iso_nonCSV_p, nbinsx)
        N_nonIso_nonCSV_p = self.get_view('QCD*').Get('nosys/looseNOTTight/csvFail/%s_Prompt' % var )
        N_nonIso_nonCSV_p = RebinView.rebin(N_nonIso_nonCSV_p, nbinsx)

        N_b_CSV = N_Iso_CSV_b+N_nonIso_CSV_b
        N_p_CSV = N_Iso_CSV_p+N_nonIso_CSV_p
        N_b_nonCSV = N_Iso_nonCSV_b+N_nonIso_nonCSV_b
        N_p_nonCSV = N_Iso_nonCSV_p+N_nonIso_nonCSV_p

        #set_trace()
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
            ## format and plot hists for fraction of leps matched to b and prompt
        plotter.teff_comparisons(N_b_CSV, N_b_CSV+N_p_CSV, xtitle=xtitle, ytitle='f_{b}^{CSV}', fig_name='Fractions_b_CSV_%s' % name )
        #plotter.plot( f_b_CSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='f_{b}^{CSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')

        plotter.teff_comparisons(N_p_CSV, N_b_CSV+N_p_CSV, xtitle=xtitle, ytitle='f_{prompt}^{CSV}', fig_name='Fractions_p_CSV_%s' % name )
        #plotter.plot( f_p_CSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='f_{prompt}^{CSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')

        plotter.teff_comparisons(N_b_nonCSV, N_b_nonCSV+N_p_nonCSV, xtitle=xtitle, ytitle='f_{b}^{nonCSV}', fig_name='Fractions_b_nonCSV_%s' % name )
        #plotter.plot( f_b_nonCSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='f_{b}^{nonCSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')

        plotter.teff_comparisons(N_p_nonCSV, N_b_nonCSV+N_p_nonCSV, xtitle=xtitle, ytitle='f_{prompt}^{nonCSV}', fig_name='Fractions_p_nonCSV_%s' % name )
        #plotter.plot( f_p_nonCSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='f_{prompt}^{nonCSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')

            ## format and plot hists for efficiencies of leps matched to b and prompt
        ##set_trace()
        plotter.teff_comparisons(N_Iso_CSV_b, N_b_CSV, xtitle=xtitle, ytitle='#epsilon_{b}^{CSV}', fig_name='Efficiencies_b_CSV_%s' % name )
        #plotter.plot( e_b_CSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='#epsilon_{b}^{CSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')

        plotter.teff_comparisons(N_Iso_CSV_p, N_p_CSV, xtitle=xtitle, ytitle='#epsilon_{prompt}^{CSV}', fig_name='Efficiencies_p_CSV_%s' % name )
        #plotter.plot( e_p_CSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='#epsilon_{prompt}^{CSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')

        plotter.teff_comparisons(N_Iso_nonCSV_b, N_b_nonCSV, xtitle=xtitle, ytitle='#epsilon_{b}^{nonCSV}', fig_name='Efficiencies_b_nonCSV_%s' % name )
        #plotter.plot( e_b_nonCSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='#epsilon_{b}^{nonCSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')

        plotter.teff_comparisons(N_Iso_nonCSV_p, N_p_nonCSV, xtitle=xtitle, ytitle='#epsilon_{prompt}^{nonCSV}', fig_name='Efficiencies_p_nonCSV_%s' % name )
        #plotter.plot( e_p_nonCSV, legend_def=LegendDefinition(position='NE'), legendstyle='p', x_range=xlims, ytitle='#epsilon_{prompt}^{nonCSV}', xtitle=xtitle, fillstyle=0, drawstyle='E0 X0')

        #set_trace()
        plotter.teff_comparisons(N_Iso_nonCSV_p+N_Iso_nonCSV_b, N_p_nonCSV+N_b_nonCSV, xtitle=xtitle,\
            ytitle='r_{QCD}=f_{b}^{nonCSV}#epsilon_{b}^{nonCSV}+f_{prompt}^{nonCSV}#epsilon_{prompt}^{nonCSV}', fig_name='r_QCD_nonCSV_%s' % name, fit=True )




        data_Iso_nonCSV_hist = self.get_view('data').Get('nosys/tight/csvFail/%s_B' % var )+self.get_view('data').Get('nosys/tight/csvFail/%s_Prompt' % var )
        data_Iso_nonCSV_hist = RebinView.rebin(data_Iso_nonCSV_hist, nbinsx)
        data_nonIso_nonCSV_hist = self.get_view('data').Get('nosys/looseNOTTight/csvFail/%s_B' % var )+self.get_view('data').Get('nosys/looseNOTTight/csvFail/%s_Prompt' % var )
        data_nonIso_nonCSV_hist = RebinView.rebin(data_nonIso_nonCSV_hist, nbinsx)
        data_nonCSV_hist = data_Iso_nonCSV_hist+data_nonIso_nonCSV_hist

            ## get non-QCD MC hists
        MC_Iso_nonCSV = sum([ i.Get('%s_B' % var) for i in plotter.mc_views(1, None, 'nosys/tight/csvFail', False) if i.Get('%s_B' % var).title != 'QCD' ])+\
                        sum([ i.Get('%s_Prompt' % var) for i in plotter.mc_views(1, None, 'nosys/tight/csvFail', False) if i.Get('%s_Prompt' % var).title != 'QCD' ])
        MC_Iso_nonCSV = RebinView.rebin(MC_Iso_nonCSV, nbinsx)
        MC_nonIso_nonCSV = sum([ i.Get('%s_B' % var) for i in plotter.mc_views(1, None, 'nosys/looseNOTTight/csvFail', False) if i.Get('%s_B' % var).title != 'QCD' ])+\
                        sum([ i.Get('%s_Prompt' % var) for i in plotter.mc_views(1, None, 'nosys/looseNOTTight/csvFail', False) if i.Get('%s_Prompt' % var).title != 'QCD' ])
        MC_nonIso_nonCSV = RebinView.rebin(MC_nonIso_nonCSV, nbinsx)
        MC_nonCSV_hist = MC_Iso_nonCSV+MC_nonIso_nonCSV


        Iso_nonCSV_hist = data_Iso_nonCSV_hist - MC_Iso_nonCSV
        plotter.plot( Iso_nonCSV_hist, xtitle=xtitle, ytitle='N_{data-nonQCD MC}^{Iso-nonCSV}', drawstyle='E0 X0')#, legend_def=LegendDefinition(position='NE'), legendstyle='p', drawstyle='E0 X0' )
        plotter.save( 'N_data_minus_MC_Iso_nonCSV_%s' % name )

        nonCSV_hist = data_nonCSV_hist - MC_nonCSV_hist
        plotter.plot( nonCSV_hist, xtitle=xtitle, ytitle='N_{data-nonQCD MC}^{nonCSV}', drawstyle='E0 X0')#, legend_def=LegendDefinition(position='NE'), legendstyle='p', drawstyle='E0 X0' )
        plotter.save( 'N_data_minus_MC_nonCSV_%s' % name )



        plotter.teff_comparisons(Iso_nonCSV_hist, nonCSV_hist, xtitle=xtitle,\
            ytitle='r_{data-nonQCD MC}^{nonCSV}=N_{data-nonQCD MC}^{Iso-nonCSV}/N_{data-nonQCD MC}^{nonCSV}', fig_name='r_data_nonCSV_%s' % name, fit=True )


        ##plotter.sf_comparisons(Iso_nonCSV_hist, nonCSV_hist, N_Iso_nonCSV_p+N_Iso_nonCSV_b, N_p_nonCSV+N_b_nonCSV, xtitle=xtitle,\
        ##    ytitle='s^{nonCSV}=r_{data-nonQCD MC}^{nonCSV}/r_{QCD}^{nonCSV}', fig_name='SF_nonCSV_%s' % name)
        #set_trace()

        #scale_factor = r_data/(f_b_nonCSV*e_b_nonCSV+f_p_nonCSV*e_p_nonCSV)
        #plotter.plot( scale_factor, x_range=xlims, xtitle=xtitle, ytitle='SF^{nonCSV}', drawstyle='E0 X0' )
        #plotter.save( 'SF_nonCSV_%s' % name )



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

    with io.root_open('muons_3jets_fits.root', 'w') as out:
        out.cd()

    lepIso_bTag = [
        ('lep_pt', 'p_{T}(l) (GeV)', 'Pt', (0, 300), [20, 40, 60, 80, 100, 500]),
        #('lep_eta', '#eta(l)', 'Eta', (-2.4, 2.4), 30)
        ('lep_eta', '#eta(l)', 'Eta', (-2.4, 2.4), [-2.4, -1.8, -1.2, -0.6, 0.0, 0.6, 1.2, 1.8, 2.4])
    ]
    for var, overlay_xtitle, name_var, xlims, nbinsx in lepIso_bTag:
        plotter.sf_simplifications(var, overlay_xtitle, name_var, xlims, nbinsx)
        #plotter.flav_fracs_and_effs(var, overlay_xtitle, name_var, xlims, nbinsx)



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
