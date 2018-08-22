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
from URAnalysis.Utilities.roottools import Envelope
import re
import itertools
import rootpy.stats as stats
import functions as fncts

parser = ArgumentParser()
parser.add_argument('mode', choices=['electrons', 'muons'], help='choose leptonic decay type')
parser.add_argument('--preselection', action='store_true', help='')
parser.add_argument('--plots', action='store_true', help='')
parser.add_argument('--flow', action='store_true', help='')
parser.add_argument('--combinedflow', action='store_true', help='Full cut_flow with mu/e info')
parser.add_argument('--shapes', action='store_true', help='')
parser.add_argument('--all', action='store_true', help='')
parser.add_argument('--btag', action='store_true', help='')
parser.add_argument('--card', action='store_true', help='')
parser.add_argument('--sysplots', action='store_true', help='dumps systematics plots, valid only if --card')
parser.add_argument('--smoothsys', action='store_true', help='')
parser.add_argument('--njets', default='', help='choose 3, 4+')
#parser.add_argument('--pdfs', action='store_true', help='make plots for the PDF uncertainties')
args = parser.parse_args()

if args.sysplots and not args.card:
    raise RuntimeError('I can only draw systematics plots if I am creating the card!')

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
        #filtering = lambda x: 'HtoTT' not in x or not os.path.basename(x).startswith('data') or \
        #     (os.path.basename(x).startswith('data') and mode[:-1] in os.path.basename(x).lower())# or \
        #     #'HtoTT' not in x

        self.tt_to_use = 'ttJets'
        self.tt_shifted = {
            'mtop_up'   : 'ttJets_mtopup',
            'mtop_down' : 'ttJets_mtopdown',
            'isr_up'   : 'ttJets_isrup',
            'isr_down' : 'ttJets_isrdown',
            'fsr_up'   : 'ttJets_fsrup',
            'fsr_down' : 'ttJets_fsrdown',
            }
        jobid = os.environ['jobid']

        files = glob.glob('results/%s/htt_simple/*.root' % jobid)
        data_files = [ fname for fname in files if (os.path.basename(fname).startswith('data') and mode[:-1] in os.path.basename(fname).lower() ) ]
        MC_files = [ fname for fname in files if not ('data' in fname ) ]
        #MC_files = [ fname for fname in files if not ('data' in fname or 'HtoTT' in fname) ]
        #files = filter(filtering, files)
        files = data_files+MC_files
        #logging.debug('files found %s' % files.__repr__())

        lumis = glob.glob('inputs/%s/*.lumi' % jobid)
        data_lumis = [ fname for fname in lumis if (os.path.basename(fname).startswith('data') and mode[:-1] in os.path.basename(fname).lower() ) ]
        MC_lumis = [ fname for fname in lumis if not ('data' in fname) ]
        #MC_lumis = [ fname for fname in lumis if not ('data' in fname or 'HtoTT' in fname) ]
        #lumis = filter(filtering, lumis)
        lumis = data_lumis+MC_lumis
        #logging.debug('lumi files found %s' % lumis.__repr__())

        self.tt_lhe_weights = prettyjson.loads(
            open('inputs/%s/ttJets.weights.json' % jobid).read()
            )
        outdir= 'plots/%s/htt/%s' % (jobid, mode)

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

            #info['view'] = views.SubdirectoryView(info['view'], mode)
            #info['unweighted_view'] = views.SubdirectoryView(info['unweighted_view'], mode)
        self.views['data']['view'] = urviews.BlindView(
            self.views['data']['view'], 
            '\w+/tight/MTHigh/(:?(:?m_tt)|(:?.+_ctstar)|(:?.+_ctstar_abs)|(:?cdelta_ld)|(:?hframe_ctheta_d))'
            )

        self.defaults = {
            'blurb' : [13, self.views['data']['intlumi']],
            'save' : {'png' : True, 'pdf' : False}
            }
        self.jobid = jobid

        #set_trace()
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
                    views.SumView(*[self.get_view(i) for i in ['[WZ][WZ]', '[WZ]Jets', 'tt[WZ]*']]),
                    #views.SumView(*[self.get_view(i) for i in ['[WZ][WZ]', 'ZJets', 'W[1-4]Jets', 'tt[WZ]*']]),
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

        ##
        ## signal sub-samples
        ## 
        added_samples = []
        #set_trace()
        for sample in self.views.keys():
            #if sample.startswith('HtoTT_'): 
            #    raise ValueError("I did not implement it yet remember to swap the H and A")
            if sample.startswith('AtoTT_') or sample.startswith('HtoTT_'):
                _, mass, width, pI = tuple(sample.split('_'))               
                samtype = 'int' if pI == 'Int' else 'sgn'
                bostype = 'ggA' if _ == 'AtoTT' else 'ggH'
                #set_trace()
                if pI == 'Int':
                    sub_name = '%s_neg-%s-%s-%s' % (bostype, samtype, width, mass)
                    self.views[sub_name] = {
                        'view' : views.ScaleView(
                            self.create_subsample(sample, ['negative'], '%s negative' % sample, color='#9999CC'),
                            -1
                            )
                        }
                    added_samples.append(sub_name)
                else:
                    sub_name = '%s_pos-%s-%s-%s' % (bostype, samtype, width, mass)
                    self.views[sub_name] = {
                        'view' : self.create_subsample(sample, ['positive'], '%s positive' % sample, color='#9999CC')
                        }
                    added_samples.append(sub_name)

        self.generic_mcs = [
            'QCD*',
            'EWK',
            'single*',          
            'ttJets_generic',
            ]

        self.mc_samples = self.generic_mcs
        self.split_mcs = [
            'AllButTT',
            'ttJets_other',
            'ttJets_unmatchable',
            'ttJets_matchable',
            'ttJets_right',
            ]

        self.card_names = {
            #synced
            'QCDmujets' if mode == 'muons' else 'QCDejets' : ['QCD*'],
            'TT' : ['ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right'],
            'VV' : ['[WZ][WZ]'],
            'TTV': ['tt[WZ]*'],
            'WJets' : ['W[0-9]Jets'],
            'ZJets' : ['ZJets'],
            'tChannel' : ['singlet_[st]channel'],
            'tWChannel' : ['single*_tW'],
            'data_obs'  : ['data']
            }
        non_qcd = ['TT', 'VV', 'TTV',   'WJets', 'ZJets',   'tChannel', 'tWChannel']

        self.card_names.update({i : [i] for i in added_samples})

        self.systematics = {
            #JES JER
            'CMS_scale_j_13TeV' : {
                'samples' : ['.*'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'jes_up'),
                '-' : lambda x: x.replace('nosys', 'jes_down'),
                'value' : 1.00,
                },
            'CMS_res_j_13TeV' : {
                'samples' : ['.*'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'jer_up'),
                '-' : lambda x: x.replace('nosys', 'jer_down'),
                'constants' : ('jer_down', 'jer_up'),
                'value' : 1.00,             
                },
            #BTAG
            'CMS_eff_b_13TeV' : {
                'samples' : ['.*'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'beff_up'),
                '-' : lambda x: x.replace('nosys', 'beff_down'),
                'value' : 1.00,
                },
            'CMS_fake_b_13TeV' : {
                'samples' : ['.*'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'bfake_up'),
                '-' : lambda x: x.replace('nosys', 'bfake_down'),
                'value' : 1.00,
                },
            #LEP Effs
            'CMS_eff_m' if mode == 'muons' else 'CMS_eff_e' : {
                'samples' : ['.*'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'lepeff_up'),
                '-' : lambda x: x.replace('nosys', 'lepeff_down'),
                'value' : 1.00,
                },

            #SCALE/PS
            'QCDscaleMERenorm_TT' : {
                'samples' : ['TT$'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'renorm_up'),
                '-' : lambda x: x.replace('nosys', 'renorm_down'),
                'value' : 1.00,
                'scales' : (1./self.tt_lhe_weights['3'], 1./self.tt_lhe_weights['6']),
                },
            'QCDscaleMEFactor_TT' : {
                'samples' : ['TT$'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'factor_up'),
                '-' : lambda x: x.replace('nosys', 'factor_down'),
                'value' : 1.00,
                'scales' : (1./self.tt_lhe_weights['1'], 1./self.tt_lhe_weights['2']),
                },
            'QCDscaleMERenormFactor_TT' : {
                'samples' : ['TT$'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'renfactor_up'),
                '-' : lambda x: x.replace('nosys', 'renfactor_down'),
                'value' : 1.00,
                'scales' : (1./self.tt_lhe_weights['4'], 1./self.tt_lhe_weights['8']),
                },
            #'QCDscalePS_TT' : {
            #    'samples' : ['TT$'],
            #    'categories' : ['.*'],
            #    'type' : 'shape',
            #    '+' : lambda x: x.replace('nosys', 'scale_up'),
            #    '-' : lambda x: x.replace('nosys', 'scale_down'),
            #    'value' : 1.00,
            #    },
            'Hdamp_TT' : {
                'samples' : ['TT$'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'hdamp_up'),
                '-' : lambda x: x.replace('nosys', 'hdamp_down'),
                'scales' : (1./self.tt_lhe_weights['240'], 1./self.tt_lhe_weights['231']),
                'value' : 1.00,
                },
            #TMass
            'TMass' : {
                'samples' : ['TT$'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'mtop_up'),
                '-' : lambda x: x.replace('nosys', 'mtop_down'),
                'value' : 1.00,
                },
            #OTHER
            'CMS_pileup' : {
                'samples' : ['.*'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'pu_up'),
                '-' : lambda x: x.replace('nosys', 'pu_down'),
                },
            'CMS_METunclustered_13TeV' : {
                'samples' : ['.*'],
                'categories' : ['.*'],
                'type' : 'shape',
                '+' : lambda x: x.replace('nosys', 'met_up'),
                '-' : lambda x: x.replace('nosys', 'met_down'),
                },
            #'pdf' : {
            #    'samples' : ['TT$'],
            #    'categories' : ['.*'],
            #    'type' : 'pdf',
            #    'pdftype' : 'nnpdf'
            #    },

            }
        self.card = None
        self.binning = {
            }

    def save_card(self, name, bbb=True):
        if not self.card:
            raise RuntimeError('There is no card to save!')
        self.card.clamp_negative_bins('.*','ZJets')
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

        if bbb:
            mcsamples = [i for i in self.card_names if i != 'data_obs']
            for catname in self.card.categories.keys():
                category = self.card.categories[catname]
                #inflate bin errors of tt for each sample
                hists = [category[i] for i in mcsamples]
                for bins in zip(*hists):
                    if bins[0].overflow: continue
                    unc = quad.quad(*[i.error for i in bins[1:]])
                    idx = bins[0].idx
                    hup = category['TT'].Clone()
                    hup[idx].value += unc
                    hdw = category['TT'].Clone()
                    hdw[idx].value -= unc
                    category['TT_CMS_httbar_%s_MCstatBin%dUp'   % (catname, idx)] = hup
                    category['TT_CMS_httbar_%s_MCstatBin%dDown' % (catname, idx)] = hdw
        self.card.save(name, self.outputdir)
    
    def pdf_unc_histo(self, view, var, pdf, central):
        nnpdf = range(9, 109)
        ## ct10 = range(112, 165)
        ## mmht = range(167, 218)
        ## allpdf = [0]+range(10, 218)
        vartemplate = '%s_mcws_pdf_%d'
        
        def to_up_down(central):
            up, down = central.Clone(), central.Clone()
            up.Reset()
            down.Reset()
            for ubin, dbin, cbin in zip(up, down, central):
                ubin.value = cbin.value + cbin.error
                dbin.value = cbin.value - cbin.error
            return up, down

        def get_and_scale(view, var, idx):
            scale = 1./self.tt_lhe_weights['%d' % idx]
            value = view.Get(vartemplate % (var, idx))
            value.Scale(scale)          
            return value

        def full_envelope(view, var, idrange):
            envelope = Envelope()
            for idx in idrange:
                envelope.add(
                    get_and_scale(view, var, idx)
                    )
            return to_up_down(envelope.one_sigma)

        def rms_envelope(view, var, idrange):
            central = view.Get(vartemplate % (var, idrange[0]))
            sqsum = central.Clone()
            sqsum.Reset()
            for idx in idrange[1:]:
                histo = get_and_scale(view, var, idx)
                for sbin, cbin, hbin in zip(sqsum, central, histo):
                    sbin.value += (cbin.value-hbin.value)**2

            nreplicas = len(idrange)-1
            for sbin, cbin in zip(sqsum, central):
                err = math.sqrt(sbin.value/nreplicas)
                cbin.error = err
            return to_up_down(central)
                    
        def quad_envelope(view, var, idrange):
            central = view.Get(vartemplate % (var, idrange[0]))
            sqsum = central.Clone()
            sqsum.Reset()
            for pos, idx in enumerate(idrange[1:]):
                if pos % 2: continue
                h1 = get_and_scale(view, var, idx)
                h2 = get_and_scale(view, var, idrange[pos+1])
                for sbin, cbin, h1bin, h2bin in zip(sqsum, central, h1, h2):
                    d1 = abs(cbin.value-h1bin.value)
                    d2 = abs(cbin.value-h2bin.value)                    
                    sbin.value = quad.quad(sbin.value, max(d1,d2))

            for sbin, cbin in zip(sqsum, central):              
                cbin.error = sbin.value
            return to_up_down(central)

        if pdf == 'nnpdf':
            up, down = full_envelope(view, var, nnpdf)
            alpha_up = get_and_scale(view, var, 110)
            alpha_dw = get_and_scale(view, var, 109)
            ret_up = central.Clone()
            ret_dw = central.Clone()
            for idx in range(len(central)):
                pdfu_delta = up[idx].value - central[idx].value
                lphu_delta = alpha_up[idx].value - central[idx].value
                ret_up[idx].value += quad.quad(pdfu_delta, lphu_delta)
                pdfd_delta = down[idx].value - central[idx].value
                lphd_delta = alpha_dw[idx].value - central[idx].value
                ret_dw[idx].value -= quad.quad(pdfd_delta, lphd_delta)
            return ret_up, ret_dw
        else:
            raise RuntimeError('invalid pdf type')

    def write_shapes(self, category_name, folder, variable,
                                     rebin=1, preprocess=None):
        '''
        should use get_shape above in some memorized form in the future
        '''
        if not self.card: self.card = DataCard('TT')
        self.card.add_category(category_name)
        category = self.card[category_name]

        path = os.path.join(folder, variable)
        card_views = {}
        for name, samples in self.card_names.iteritems():
            card_views[name] = self.rebin_view(
                views.SumView(
                    *[self.get_view(i) for i in samples]
                    ),
                rebin
                )
            if 'ggA_neg' in name or 'ggH_neg' in name:
                card_views[name] = views.ScaleView(
                    card_views[name], 1)

        if preprocess:
            for name in card_views:
                card_views[name] = preprocess(card_views[name])

        for name, view in card_views.iteritems():
            histo = view.Get(path)
            integral = histo.Integral()
            
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
                    hu, hd = self.pdf_unc_histo(
                        views.SubdirectoryView(view, folder),
                        variable, info['pdftype'], category[name]
                        )
                    category['%s_%sUp'   % (name, sys_name)] = hu
                    category['%s_%sDown' % (name, sys_name)] = hd                   
                    plotter.card.add_systematic(sys_name, 'shape', category_name, name, 1.00)
                elif info['type'] == 'shape' or 'value' not in info:
                    if 'mass_sf' in info:
                        #use SF instead of the shape for smoother templates
                        #get central values
                        hup = category[name].Clone()
                        hdw = category[name].Clone()
                        idx=0;
                        nsfs = len(info['mass_sf']['+'])
                        for ubin, dbin in zip(hup, hdw):
                            if ubin.overflow: continue
                            ubin.value *= info['mass_sf']['+'][idx % nsfs]
                            dbin.value *= info['mass_sf']['-'][idx % nsfs]
                            idx += 1
                    else:
                        path_up = info['+'](path) if '+' in info else None
                        path_dw = info['-'](path) if '-' in info else None
                        hup = view.Get(path_up) if path_up else None
                        hdw = view.Get(path_dw) if path_dw else None
                    if 'scales' in info:
                        sup, sdw = info['scales']
                        hup.Scale(sup)
                        hdw.Scale(sdw)

                    if hup is None and hdw is None:
                        raise RuntimeError('%s systematic does not define neither "+" nor "-" values' % sys_name)
                    elif hup is None or hdw is None: 
                        mirrored = histo.Clone()
                        multiplier = integral 
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

                    #shapes only: store shape in root file                  
                    if info['type'] == 'shape':
                        if hup.Integral() == 0.:
                            rootpy.log["/"].warning('%s Up for %s/%s has normalization == 0, forcing it to 10**-6' %(sys_name, category_name, name))
                            mbin = category[name].GetMaximumBin()
                            hup[mbin].value = 10**-6
                        if hdw.Integral() == 0.:
                            rootpy.log["/"].warning('%s Down for %s/%s has normalization == 0, forcing it to 10**-6' %(sys_name, category_name, name))
                            mbin = category[name].GetMaximumBin()
                            hdw[mbin].value = 10**-6
                        category['%s_%sUp'  % (name, sys_name)] = hup
                        category['%s_%sDown' % (name, sys_name)] = hdw                  


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

        if args.combinedflow:
            data = self.get_view('data').Get('cut_flow')
            for idx in range(1,len(data)):
                cflow[idx][data.title] = data[idx].value
            smin = min(stack.min(), data.min(), 1.2)
            smax = max(stack.max(), data.max())
            histo.yaxis.range_user = smin*0.8, smax*1.2
            xbin_min = [ binx for binx in range(1, histo.GetNbinsX()+1) if 'trigger' in histo.GetXaxis().GetBinLabel(binx)][0]
            xbin_max = [ binx for binx in range(1, histo.GetNbinsX()+1) if '%s/object selection' %  args.mode in histo.GetXaxis().GetBinLabel(binx)][0]
            histo.xaxis.range_user = xbin_min-1, xbin_max
            stack.Draw('same')
            data.Draw('same')
            self.keep.append(data)
            self.add_legend([stack, data], False, entries=len(views_to_flow)+1)
            self.pad.SetLogy()
            self.add_ratio_plot(data, stack, x_range=(xbin_min-1, xbin_max), ratio_range=0.4)
            self.lower_pad.SetLogy(False)

        if not args.combinedflow:
            data = self.get_view('data').Get('cut_flow')
            for idx in range(1,len(samples[0])):
                cflow[idx][data.title] = data[idx].value
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

        set_trace()

        return cflow
            #cut_flow.GetYaxis().SetRangeUser(1, 10**7)

    def make_preselection_plot(self, *args, **kwargs):
        systematics = None
        if 'sys_effs' in kwargs:
            systematics = kwargs['sys_effs']
            del kwargs['sys_effs']
        mc_default = self.mc_samples
        self.mc_samples = [
            'QCD*',
            'EWK',
            'single*',          
            'ttJets_preselection'
            ]
        #set_trace()
        self.plot_mc_vs_data(*args, **kwargs)
        if systematics is not None:
            data = [i for i in self.keep if isinstance(i, ROOT.TH1)][0]
            mc_stack = [i for i in self.keep if isinstance(i, ROOT.THStack)][0]
            stack_sum = sum(mc_stack.hists)
            set_trace()
            self.reset()
            dirname = args[0].split('/')[0]
            path = args[0]
            args = list(args)
            args[0] = path.replace(dirname, '%s_up' % systematics)
            kwargs['nodata'] = True
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
            # Add legend
            self.pad.cd()
            self.add_legend(
                [mc_stack, stack_sum, data], kwargs.get('leftside', True), 
                entries=len(mc_stack.hists)+2
                )

        self.mc_samples = mc_default

    def print_table(self, lines, separate_head=True):
        #"""Prints a formatted table given a 2 dimensional array"""
        #Count the column width
        widths = []
        for line in lines:
            for i, size in enumerate([len(x) for x in line]):
                while i <= len(widths):
                    widths.append(0)
                if size < widths[i]:
                    widths[i] = size
        
        #Generate the format string to pad the columns
        print_string = ""
        for i, width in enumerate(widths):
            print_string += "{" + str(i) + ":" + str(width) + "} | "
        if (len(print_string) == 0):
            return
        print_string = print_string[:-3]
        
        #Print the actual data
        for i, line in enumerate(lines):
            print(print_string.format(*line))
            if (i==0 and separate_head):
                print("-"*(sum(widths)+3*(len(widths)-1)))
        

    def make_sample_table(self, threshold=None, absolute=False, fname='yields.raw_txt'):
        stack = [i for i in self.keep if isinstance(i, plotting.HistStack)][0]
        names = [i.title for i in stack.hists]
        yields = [i.Integral() for i in stack.hists]
        #set_trace()
        data = [i for i in self.keep if hasattr(i, 'title') and i.title == 'Observed'][0]
        def ratio(lst):
            tot = sum(lst)/100. if not absolute else 1.
            return [i/tot for i in lst] if tot else [0. for _ in lst]
        yields = ratio(yields)
        if absolute:
            yields.append(data.Integral())
            names.append('Observed')


#       rows = []
#       rows.append(("Sample", "Fraction"))
#       for i in range(len(names)):
#           print "Name: ", names[i]
#           print "Yield: ", yields[i]
#           #rows.append((names[i], format(yields[i], '.1f')))
#       #plotter.print_table(rows)

        mlen = max(len(i) for i in names)
        format = '%'+str(mlen)+('s     %.1f%%\n' if not absolute else 's     %.0f\n')
        header = ('%'+str(mlen)+'s    frac') % 'sample'
        if threshold is not None:
            binid = stack.hists[0].xaxis.FindBin(threshold)
            fbin = stack.hists[0].nbins()+1
            less = ratio([i.Integral(1, binid-1) for i in stack.hists])
            above = ratio([i.Integral(binid, fbin) for i in stack.hists])
            if absolute:
                less.append(data.Integral(1, binid-1))
                above.append(data.Integral(binid, fbin))
                format = format.replace('\n', '    %.0f    %.0f\n')
            else:
                format = format.replace('\n', '    %.1f%%    %.1f%%\n')
            yields = [i for i in zip(yields, less, above)]
            header += '    <%.0f %s     >%.0f %s' % (threshold, stack.hists[0].xaxis.title, threshold, stack.hists[0].xaxis.title)
        header += '\n'
        fullname = os.path.join(self.outputdir, fname)
        with open(fullname, 'w') as f:
            f.write('Lumi :%.1f' % self.views['data']['intlumi'])
            f.write(header)
            for i in zip(names, yields):
                info = i
                if threshold is not None:
                    #repack info
                    a, b = i
                    c, d, e = b
                    info = (a, c, d, e)
                try:
                    f.write(format % info)
                except:
                    set_trace()

    def get_yields(self, thr):
        stack = [i for i in self.keep if isinstance(i, plotting.HistStack)][0]
        obs = self.keep[1]
        assert(obs.title == 'Observed')
        hists = stack.hists
        hists.append(obs)
        ret = []
        for hist in hists:
            name = hist.title
            binid = hist.xaxis.FindBin(thr)
            fbin = hist.nbins()+1
            loe = ROOT.Double()
            lov = hist.IntegralAndError(1, binid-1, loe)
            hie = ROOT.Double()
            hiv = hist.IntegralAndError(binid, fbin, hie)
            ret.append(
                (name, (lov, loe), (hiv, hie))
                )
        return ret

    #def create_mtt_root(self, fname):
    #    # create and then move to correct directory to get full tdir path
    #    with io.root_open(fname, 'w') as out:
    #        for lepiso in ['looseNOTTight', 'tight']:
    #            out.mkdir(lepiso).cd()
    #            for mt in ['MTHigh', 'MTLow']:
    #                out.mkdir(lepiso+'/'+mt).cd()


    #def make_mtt_plots(self, tdir, fname):
    #    hists = [ i.Get('njets') for i in plotter.mc_views(1, None, 'nosys/%s' % tdir) ]
    #    if tdir != 'tight/MTHigh':
    #        data = self.get_view('data').Get('nosys/%s/njets' % tdir)
    #        hists.append(data)
    #    names = [i.title for i in hists]
    #    yields = [i.Integral() for i in hists]
    #    [ hists[i].set_name(hists[i].title) for i in range(len(hists))] # set hist name to QCD, EW, Single...
    #    mtt_dict[tdir] = hists
    #    #set_trace()
    #    #open root file
    #    with io.root_open(fname, 'update') as out:
    #        #set_trace()
    #        #move to correct directory to get full tdir path
    #        out.cd(tdir)
    #        for hist in hists:
    #            hist.Write()



    def get_abcd_scale(self, var):

        N_B = 0
        N_C = 0
        N_D = 0
        for dirid in itertools.product(['looseNOTTight', 'tight'], ['MTHigh', 'MTLow']):
            tdir = '%s/%s' % dirid
            if tdir == 'tight/MTHigh':
                continue
            #set_trace()
            mc_hists = [ i.Get(var) for i in plotter.mc_views(1, None, 'nosys/%s' % tdir, False) if i.Get(var).title != 'QCD' ]
            mc_names = [i.title for i in mc_hists]
            [ mc_hists[i].set_name(mc_hists[i].title) for i in range(len(mc_hists))] # set hist name to QCD, EW, Single...
            data_hist = self.get_view('data').Get('nosys/%s/%s' % (tdir, var))
            data_hist.set_name(data_hist.title)

            ### subtract all mc (except QCD) from data in all regions
            N = data_hist.Integral() - sum(mc_hists).Integral()
            if tdir == 'tight/MTLow':
                N_B = N
            if tdir == 'looseNOTTight/MTLow':
                N_D = N
            if tdir == 'looseNOTTight/MTHigh':
                N_C = N
            #set_trace()

        #set_trace()
        #### get scale Na = Nb*(Nc/Nd)
        N_A = N_B*(N_C/N_D)
        return N_A
        #set_trace()


    #def qcd_norm(self, fname, scale, var):
    #    shape_region = 'tight/MTLow'
    #    signal_region = 'tight/MTHigh'
    #    qcd_shape_hist = [ i.Get(var) for i in plotter.mc_views(1, None, 'nosys/%s' % shape_region) if i.Get(var).title == 'QCD' ][0]
    #    #N_C = 0
    #    #N_D = 0
    #    set_trace()
    #    #for key in mtt_dict.keys():

    #    #    if key == 'tight/MTHigh':
    #    #        continue
    #    #    ## find indices for observed and ttbar hists in dictionary
    #    #    obs_idx = [ n for n,i in enumerate(mtt_dict[key]) if i.name == 'Observed' ][0]
    #    #    tt_idx = [ n for n,i in enumerate(mtt_dict[key]) if i.name == 't#bar{t}' ][0]
    #    #    singlet_idx = [ n for n,i in enumerate(mtt_dict[key]) if i.name == 'single top' ][0]
    #    #    ew_idx = [ n for n,i in enumerate(mtt_dict[key]) if i.name == 'EW' ][0]

    #    #    ### subtract all mc (except QCD) from data in all regions
    #    #    N = (mtt_dict[key][obs_idx]-(mtt_dict[key][tt_idx]+mtt_dict[key][singlet_idx]+mtt_dict[key][ew_idx])).Integral()

    #    #    if key == 'tight/MTLow':
    #    #        N_B = N
    #    #    if key == 'looseNOTTight/MTLow':
    #    #        N_D = N
    #    #    if key == 'looseNOTTight/MTHigh':
    #    #        N_C = N
    #    #        #set_trace()

    #    #### get scale Na = Nb*(Nc/Nd)
    #    #N_A = N_B*(N_C/N_D)
    #    ##set_trace()
    #    #qcd_shape = [ i for i in mtt_dict['looseNOTTight/MTHigh'] if i.name == 'QCD' ][0]  # get shape of QCD from region C
    #    #qcd_hist = qcd_shape*(N_A/qcd_shape.Integral())
    #    #qcd_hist.set_name('QCD_est')
    #    #qcd_hist.set_title('QCD Estimation')
    #    #with io.root_open(fname, 'update') as out:
    #    #    qcd_hist.Write()
    #    #set_trace()


plotter = HTTPlotter(args.mode)

variables = [
  (False, "m_tt"    , "m(t#bar{t}) (GeV)", 1, None, False),       
  (False, "tlep_ctstar" , "cos #theta^{*}(t_{lep})", 2, None, False),      
  (False, "thad_ctstar" , "cos #theta^{*}(t_{had})", 2, None, False),      
  (False, "tlep_ctstar_abs" , "|cos #theta^{*}(t_{lep})|", 2, None, False),      
  (False, "thad_ctstar_abs" , "|cos #theta^{*}(t_{had})|", 2, None, False),      
  (False, "cdelta_ld", "cos #delta(ld)", 2, None, False),       
  (False, "hframe_ctheta_d", "cos #theta(d-jet)", 2, None, False),      
  (False, "lead_jet_pt" , 'leading jet p_{T}', 10, (0, 300), False),
  (False, "lead_jet_eta", 'leading jet #eta' ,  5, None, False),
  (False, "pt_thad" , "p_{T}(t_{had}) (GeV)", 5, (0, 1000), False),     
  (False, "eta_thad"    , "#eta_{T}(t_{had}) (GeV)", 5, None, False),       
  (False, "pt_tlep" , "p_{T}(t_{lep}) (GeV)", 5, (0, 1000), False),     
  (False, "eta_tlep"    , "#eta_{T}(t_{lep}) (GeV)", 5, None, False),       
  (False, "pt_tt"   , "p_{T}(t#bar{t}) (GeV)", 5, (0, 1000), False),      
  (False, "eta_tt"  , "#eta_{T}(t#bar{t}) (GeV)", 2, None, False),        
  (False, "full_discriminant_3j"    , "#lambda_{comb} 3 jets",   2, (0, 25), False),        
  (False, "full_discriminant_4j"    , "#lambda_{comb} 4 jets",   2, (5, 30), False),        
  (False, "full_discriminant_5j"    , "#lambda_{comb} 5 jets",   2, (5, 30), False),        
  (False, "full_discriminant_6Pj"   , "#lambda_{comb} #geq 6 jets", 2, (5, 30), False),     
]

preselection = [
    (False, "MT" , "M_{T}",  10, (0, 300), False),
    (False, "njets"  , "# of selected jets", range(13), None, False),
    (False, "jets_eta", "#eta(jet)", 5, None, False),
    (False, "jets_pt", "p_{T}(jet)", 10, (0, 300), False),
    (False, "lep_eta", "#eta(l)", 5, None, False),
    (False, "lep_pt"    , "p_{T}(l) (GeV)", 20, (0, 300), False),       
    (False, "nvtx", "# of reconstructed vertices", range(41), None, False),
    (False, "rho", "#rho", range(40), None, False),
    (False, "lep_iso", 'l rel Iso', 1, [0,1], False),
    (False, "lep_wp" , "electron wp", 1, None, False),
    (True   , "csvv2"    , "csv",  1, None, False),
    (True   , "csvv2_p11", "csv^{11}", 1, None, False),
    (False, "METPhi", "MET #varphi", 4, None, False),
    (False, "MET"   , "MET E_{T}"  , 1, [0, 400], False),
]

permutations = [
    (False, "bjets_pt", "p_{T} (b-jets)", 10, (0, 300), False),
    (False, "wjets_pt", "p_{T} (W-jets)", 10, (0, 200), False),
    (False, "full_discriminant", "#lambda_{C}", 1, (0, 30), False),
    (False, "nu_discriminant", "#chi^{2}_{#nu}", 1, None, False),
    (False, "mass_discriminant", "#lambda_{M}", 1, None, False),
    (False, "Wmasshad", "M(W_{h})", 2, (0,400), False),
    (False, "tmasshad", "M(t_{h})", 2, None, False),
    (False, "lbratio" , "p_{T}(l)/p_{T}(b)", 1, (0, 6), False),
    (False, "j2bratio", "p_{T}(j)/p_{T}(b)", 1, (0, 6), False),
]

jet_categories = ["3jets", "4jets", "5Pjets"]

#cut flow

if args.flow or args.all:
   flow = plotter.cut_flow()
   if args.combinedflow:
      plotter.set_outdir('plots/%s/htt' % plotter.jobid)
   plotter.save('cut_flow')
   rootpy.log["/"].warning('Be careful interpreting these plots, bins may not be in the same order!!')


if args.preselection or args.all:
    if args.njets == '3':
        plotter.set_subdir('3Jets/preselection')
    elif args.njets == '4+':
        plotter.set_subdir('4PJets/preselection')
    else:
        raise RuntimeError('Your choice for --njets is invalid!')
        #plotter.set_subdir('Incl/preselection')
    for logy, var, axis, rebin, x_range, leftside in preselection + permutations:
        plotter.make_preselection_plot(
            'nosys/preselection', var, sort=True,
            xaxis=axis, leftside=leftside, rebin=rebin, 
            show_ratio=True, ratio_range=0.5, xrange=x_range, logy=logy)        
        plotter.save(var)

if args.plots or args.all:

    #qcd_renorm = True
    qcd_renorm = False
    qcd_scale = 1.
    if qcd_renorm:
        qcdd_scale = plotter.get_abcd_scale('njets')
    #set_trace()
    #    ## create file to save mtt hists for abcd method
    #if args.njets == '3':
    #    plotter.create_mtt_root(fname='3J_mtt_plots.root')
    #elif args.njets == '4+':
    #    plotter.create_mtt_root(fname='4PJ_mtt_plots.root')
    #else:
    #    plotter.create_mtt_root(fname='Incl_mtt_plots.root')
    #mtt_dict = {}
    vals = []
    #for dirid in itertools.product(['looseNOTTight', 'tight'], ['MTHigh']):
    for dirid in itertools.product(['looseNOTTight', 'tight'], ['MTHigh', 'MTLow']):
        tdir = '%s/%s' % dirid
        if args.njets == '3':
            plotter.set_subdir('3Jets/'+tdir)
        elif args.njets == '4+':
            plotter.set_subdir('4PJets/'+tdir)
        else:
            raise RuntimeError('Your choice for --njets is invalid!')
            #plotter.set_subdir('Incl/'+tdir)
        #plotter.set_subdir(tdir)
        #set_trace()
        first = True
        #for logy, var, axis, rebin, x_range, leftside in variables:
        for logy, var, axis, rebin, x_range, leftside in preselection+variables+permutations:
            if 'discriminant' in var:
                plotter.mc_samples = plotter.split_mcs
                if 'mass' in var: rebin = [0, 5, 10, 15, 20]
            plotter.plot_mc_vs_data(
                'nosys/%s' % tdir, var, sort=True,
                xaxis=axis, leftside=leftside, rebin=rebin,
                show_ratio=True, ratio_range=0.2, xrange=x_range,
                logy=logy, qcd_renorm=qcd_renorm, qcd_scale=qcd_scale)
            #set_trace()
            if first:
                first = False
                plotter.make_sample_table(threshold=50)
                plotter.make_sample_table(threshold=50, absolute=True, fname='yields_abs.raw_txt')
                vals.append(
                    (tdir, plotter.get_yields(50))
                    )   
            #    #set_trace()
            if '3j' in var:
                plotter.make_sample_table(threshold=50, fname='Discriminant_3J_yields.raw_txt')
                vals.append(
                    (tdir, plotter.get_yields(50))
                    )   
            elif '4j' in var:
                plotter.make_sample_table(threshold=50, fname='Discriminant_4J_yields.raw_txt')
                vals.append(
                    (tdir, plotter.get_yields(50))
                    )   
            elif '5j' in var:
                plotter.make_sample_table(threshold=50, fname='Discriminant_5J_yields.raw_txt')
                vals.append(
                    (tdir, plotter.get_yields(50))
                    )   
            elif '6Pj' in var:
                plotter.make_sample_table(threshold=50, fname='Discriminant_6PJ_yields.raw_txt')
                vals.append(
                    (tdir, plotter.get_yields(50))
                    )   
            print var
            #set_trace()
            plotter.save(var)
            if 'discriminant'in var:
                plotter.mc_samples = plotter.generic_mcs

    jmap = {}
    for tdir, sams in vals:
        for sam, lo, hi in sams:
            if sam not in jmap:
                jmap[sam] = {}
            jmap[sam]['%s/%s' % (tdir, 'hi')] = hi
            jmap[sam]['%s/%s' % (tdir, 'lo')] = lo
    with open('yields.json', 'w') as out:
        out.write(prettyjson.dumps(jmap))


if args.shapes or args.all:
    for peak, dname in [('*', 'shapes'), ('Peak', 'shapes_peak'), ('Int', 'shapes_interference')]:
        if args.njets == '3':
            plotter.set_subdir('3Jets/'+dname)
        elif args.njets == '4+':
            plotter.set_subdir('4PJets/'+dname)
        else:
            raise RuntimeError('Your choice for --njets is invalid!')
        #plotter.set_subdir(dname)
        for boson in ['AtoTT', 'HtoTT']:
            for mass in [400, 500, 600, 750]:
                histos = []
                for width, color in zip([5, 10, 25, 50], ['#f9a505', '#2aa198', '#0055ff', '#6666b3']):
                    htt_view = plotter.get_view('%s_M%d_%dpc_%s' % (boson, mass, width, peak))
                    histos.append(
                        sum(
                            htt_view.Get('%s/nosys/tight/MTHigh/m_tt' % i) \
                                for i in ['positive', 'negative']
                                #for i in ['right', 'matchable', 'unmatchable', 'noslep']
                            )
                        )
                    histos[-1].Rebin(2)
                    histos[-1].linecolor = color
                    if histos[-1].Integral():
                        histos[-1].Scale(1/abs(histos[-1].Integral()))
                plotter.overlay(histos, xtitle='m(tt) (GeV)', ytitle='a.u.', x_range=(0, 1300))
                plotter.add_legend(histos)
                plotter.save('%s_M%s' % (boson, mass) )

if args.btag:
    mc_default = plotter.mc_samples
    plotter.mc_samples = [
        '[WZ][WZ]', 'QCD*', '[WZ]Jets', 
        'single*', 'ttJets_preselection'
        ]
    hists = [i.Get('jets_csvv2_WP') for i in plotter.mc_views(1, None, 'nosys/preselection')]
    plotter.mc_samples = mc_default
    ttb = hists[-1]
    bkg = sum(hists[:-1])
    qcd = [i for i in hists if i.title == 'QCD'][0]
    rest = sum(i for i in hists if i.title != 'QCD')
    sig_sqrt_bkg = bkg.Clone()
    sig_sqrt_sb  = bkg.Clone()
    qcd_part = bkg.Clone()
    ttfrac = bkg.Clone()
    labels = [None, 'None', 'Loose', 'Medium', 'Tight']
    rows = []
    rows.append(("WP", "WP", "QCD", "Other MC", "QCD Frac"))
    for xstrt, ystrt in itertools.product(range(1,5), range(1,5)):
        s, b = 0, 0
        q, e = 0, 0
        sig_sqrt_sb.xaxis.SetBinLabel(xstrt, labels[xstrt])
        sig_sqrt_sb.yaxis.SetBinLabel(ystrt, labels[ystrt])
        sig_sqrt_bkg.xaxis.SetBinLabel(xstrt, labels[xstrt])
        sig_sqrt_bkg.yaxis.SetBinLabel(ystrt, labels[ystrt])
        qcd_part.xaxis.SetBinLabel(xstrt, labels[xstrt])
        qcd_part.yaxis.SetBinLabel(ystrt, labels[ystrt])
        if not ttb[xstrt, ystrt].value and not bkg[xstrt, ystrt].value: continue
        for x, y in itertools.product(range(xstrt,5), range(ystrt,5)):
            s += ttb[x,y].value
            b += bkg[x,y].value
            q += qcd[x,y].value
            e += rest[x,y].value
        print labels[xstrt], labels[ystrt], q, e, q / (q+e)
        rows.append( (labels[xstrt], labels[ystrt], format(q), format(e), format(q / (q+e)) ) )
        sig_sqrt_bkg[xstrt, ystrt].value = s/math.sqrt(b)
        sig_sqrt_sb[ xstrt, ystrt].value = s/math.sqrt(s+b)
        ttfrac[xstrt, ystrt].value = s
        qcd_part[xstrt, ystrt].value = q / (q+e)
    ttfrac.Scale(1/ttb.Integral())
    if args.njets == '3':
        plotter.set_subdir('3Jets/preselection/btag')
    elif args.njets == '4+':
        plotter.set_subdir('4PJets/preselection/btag')
    else:
        raise RuntimeError('Your choice for --njets is invalid!')
    plotter.plot(sig_sqrt_bkg, drawstyle='colz')
    plotter.save('sig_sqrt_bkg')
    plotter.plot(sig_sqrt_sb, drawstyle='colz')
    plotter.save('sig_sqrt_sb')
    plotter.plot(qcd_part, drawstyle='colz')
    plotter.save('qcd_contamination')
    plotter.plot(ttfrac, drawstyle='colz')
    plotter.save('fraction_tt')
    fncts.print_table(rows, filename='%s/WP_Fracs.raw_txt' % plotter.outputdir )

binnind2D = (
    [250.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 610.0, 640.0, 680.0, 730.0, 800.0, 920.0, 1200.0], #~3k events each mtt bin
    [1] #[0, 0.2, 0.4, 0.6, 0.8, 1.0]
)
if args.card:
    correction_factors = {}
    category = 'mujets' if args.mode == 'muons' else 'ejets'
    njets = '3Jets' if args.njets == '3' else '4PJets'
    if args.smoothsys:
        #systematics = args.smoothsys.split(',')
        systematics = plotter.systematics.keys()
        #set_trace()
        for shift in systematics:
            if not shift in plotter.systematics:
                raise KeyError(
                    'there is no systematic called:'
                    ' %s, all I have is: %s' % (shift, ', '.join(plotter.systematics.keys())),
                    )
            
            path = 'nosys/tight/MTHigh/m_tt'
            view = plotter.rebin_view(
                views.SumView(
                    *[plotter.get_view(i) for i in plotter.card_names['TT']]
                    ),
                binnind2D[0]
                )
            central = view.Get(path)
            up   = view.Get(plotter.systematics[shift]['+'](path))
            up /= central
            ufcn = plotting.F1('pol4')
            ufcn.decorate(linecolor='blue', linewidth=2)
            up.fit(ufcn)
            down = view.Get(plotter.systematics[shift]['-'](path))
            down /= central
            dfcn = plotting.F1('pol4')
            dfcn.decorate(linecolor='red', linewidth=2)
            down.fit(dfcn)
            newd = down*-1
            plotter.set_subdir('%s/shapes/%s' % (njets, category) )
            with io.root_open('smoothsys_test.root', 'w') as out:
                newd.name = 'down'
                up.name = 'up'
                out.WriteTObject(newd, '%s_down' % shift)
                out.WriteTObject(up, '%s_up' % shift)
            
            #set_trace()
            plotter.overlay(
                [up, down], linecolor=['blue', 'red'], y_range=(0.8,1.2),
                markercolor=['blue', 'red'], title=['up', 'down'],
                legend_def=LegendDefinition(position='NW'), legendstyle='p',
                markerstyle=20, fillstyle='hollow',
                drawstyle='E0', markersize=0.5,
                xtitle='m(t#bar{t}) (GeV)')
            plotter.save('%s_unsmoothed' % shift)
            print "Symmetry test:"
            for idx in range(ufcn.GetNpar()):
                delta = abs(ufcn[idx].value+dfcn[idx].value)/quad.quad(ufcn[idx].error,dfcn[idx].error)
                print '  %d: %.2f' % (idx, delta)
            upsf, downsf = [], []
            for bin in central:
                if bin.overflow: continue
                usf = ufcn(bin.x.center)
                dsf = dfcn(bin.x.center)
                if usf > 1 and dsf > 1:
                    print 'warning'
                if usf < 1 and dsf < 1:
                    print 'warning'                 
                upsf.append(ufcn(bin.x.center))
                downsf.append(dfcn(bin.x.center))
            plotter.systematics[shift]['mass_sfs'] = {
                '+' : upsf, '-' : downsf
                }
    #raise ValueError()
    plotter.write_shapes(       
        category,
        'nosys/tight/MTHigh', 'mtt_tlep_ctstar_abs',
        rebin=binnind2D, preprocess=urviews.LinearizeView)  
    ## plotter.write_shapes(        
    ##  'mujets' if args.mode == 'muons' else 'ejets',
    ##  'nosys/tight/MTHigh', 'm_tt',
    ##  rebin=binnind2D[0])
    ## plotter.write_shapes(        
    ##  'mujets' if args.mode == 'muons' else 'ejets',
    ##  'nosys/tight/MTHigh', 'tlep_ctstar',
    ##  rebin=10)
    plotter.set_subdir(njets)
    #set_trace()
    plotter.save_card(args.mode)
    if args.sysplots:
        categories = plotter.card.categories.keys()
        for category in categories:
            plotter.set_subdir('%s/shapes/%s' % (njets, category) )
            category = plotter.card.categories[category]
            samples  = category.keys()
            shifted_tt = [i for i in samples if i.startswith('TT_') and i.endswith('Up')]
            ttsample = category.TT
            ttsample.linecolor = 'black'
            ttsample.linewidth = 2
            ttsample.drawstyle = 'hist'
            ttsample.fillstyle = 'hollow'
            for i in ttsample:
                i.error = 0
            for shifted in shifted_tt:
                if 'MCstatBin' in shifted: continue
                name = shifted[3:-2]
                up = category[shifted]
                down = category[shifted.replace('Up','Down')]
                center = (up+down)/2.
                delta = (up-down)/2.
                for cbin, dbin in zip(center, delta):
                    cbin.error = dbin.value
                center.fillcolor = '#b3ecfd'
                center.linewidth  = 0
                center.markersize = 0
                center.drawstyle = 'E2'
                plotter.overlay_and_compare(
                        [center], ttsample, method='ratio', lower_y_range='auto',
                        legend_def=None, xtitle='mass discriminant', ytitle='events',
                        ignore_style=False, ratio_style = {
                        'fillcolor' : '##b3ecfd',
                        'linewidth' : 0,
                        'markersize': 0,
                        'drawstyle' : 'E2',
                        'fillstyle' : 'solid'
                        }
                        )
                plotter.save(name)

