#! /bin/env python

import math
import ROOT
import rootpy.plotting as plotting
from rootpy.plotting import Hist
import rootpy
import rootpy.io as io
from URAnalysis.PlotTools.views import EnvelopeView
from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
from URAnalysis.Utilities.roottools import ArgSet, ArgList
from pdb import set_trace
import logging
from os.path import join
from labels import set_pretty_label
import os
import URAnalysis.Utilities.prettyjson as prettyjson
import re

asrootpy = rootpy.asrootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(11111)
ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)
log = rootpy.log["/toy_diagnostics"]
log.setLevel(rootpy.log.INFO)

tau_nbins = 20
pull_nbins = 20
pull_mean_nbins = 20
pull_sigma_nbins = 20

delta_nbins = 20
delta_mean_nbins = 20
delta_sigma_nbins = 20

unfolded_nbins = 20



def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]


def create_and_save_canvas(histo, xname):
    canvas = plotting.Canvas()
    histo.SetStats(True)
    histo.GetXaxis().SetName(xname)
    histo.Draw()
    canvas.Update()
    canvas.SaveAs('%s.png' % histo.GetName())
    canvas.SaveAs('%s.pdf' % histo.GetName())

def unfolding_toy_diagnostics(indir, variable):
    
    plotter = BasePlotter(
        defaults={
            'clone' : False,
            'name_canvas' : True,
            'show_title' : True,
            'save' : {'png' : True, 'pdf' : False}
            },
        )
    styles = {
        'dots' : {
            'linestyle' : 0, 
            'markerstyle' : 21, 
            'markercolor' : 1
            },
        'compare' : {
            'linesstyle' : [1,0],
            'markerstyle' : [0,21],
            'markercolor' : [2,1],
            'linecolor' : [2,1],
            'drawstyle' : ['hist', 'pe'],
            'legendstyle' : ['l', 'p']
            }
        }
    
    xaxislabel = set_pretty_label(variable)
    
    true_distribution = None
    

    curdir = os.getcwd()
    os.chdir(indir)
    toydirs = get_immediate_subdirectories(".")
    
    methods = []
    pulls_lists = {}
    pull_means_lists = {}
    pull_mean_errors_lists = {}
    pull_sums_lists = {}
    pull_sigmas_lists = {}
    pull_sigma_errors_lists = {}
    deltas_lists = {}
    delta_means_lists = {}
    delta_mean_errors_lists = {}
    delta_sigmas_lists = {}
    delta_sigma_errors_lists = {}
    ratio_sums_lists = {}
    nneg_bins_lists = {}
    unfoldeds_lists = {}
    unfolded_sigmas_lists = {}
    taus_lists = {}
    
    histos_created = False
    lists_created = False
    idir = 0
    true_distro = None
    #loop over toys
    for directory in toydirs:
        if not directory.startswith('toy_'): continue
        os.chdir(directory)
        log.debug('Inspecting toy %s' % directory)
        idir = idir + 1
        i = 0
        if not os.path.isfile("result_unfolding.root"):
            raise ValueError('root file not found in %s' % os.getcwd())
        with io.root_open("result_unfolding.root") as inputfile:
            log.debug('Iteration %s over the file' % i)
            i = i+1
            if not methods:
                keys = [i.name for i in inputfile.keys()]
                for key in keys:
                    if hasattr(getattr(inputfile, key), "hdata_unfolded"):
                        methods.append(key)
                    
            unfolded_hists = [inputfile.get('%s/hdata_unfolded' % i) for i in methods]
            unfolded_wps_hists = [inputfile.get('%s/hdata_unfolded_ps_corrected' % i) for i in methods]
            for unf, unfps, method in zip(unfolded_hists, unfolded_wps_hists, methods):
                unf.name = method
                unfps.name = method
            if true_distro is None:
                true_distribution = inputfile.true_distribution
                ROOT.TH1.AddDirectory(False)
                true_distro = true_distribution.Clone()
            taus = prettyjson.loads(inputfile.best_taus.GetTitle())
            if len(taus_lists) == 0:
                taus_lists = dict((i, []) for i in taus)
            for i, t in taus.iteritems():
                taus_lists[i].append(t)

            for histo in unfolded_hists:
                #create pull/delta containers during first iteration
                name = histo.name
                nbins = histo.nbins()
                log.debug("name = %s, n bins = %s" % (name,nbins))
                if not lists_created:
                    for ibin in range(1,nbins+1):
                        outname = "pull_" + name + "_bin" + str(ibin)
                        pulls_lists[outname] = []
                        outname = "delta_" + name + "_bin" + str(ibin)
                        deltas_lists[outname] = []
                        outname = "unfolded_" + name + "_bin" + str(ibin)
                        unfoldeds_lists[outname] = []
                        unfolded_sigmas_lists[outname] = []
                    outname = "pull_" + name
                    pull_means_lists[outname] = {}
                    pull_mean_errors_lists[outname] = {}
                    pull_sigmas_lists[outname] = {}
                    pull_sigma_errors_lists[outname] = {}
                                            
                    outname = "delta_" + name
                    delta_means_lists[outname] = {}
                    delta_mean_errors_lists[outname] = {}
                    delta_sigmas_lists[outname] = {}
                    delta_sigma_errors_lists[outname] = {}
                    
                    
                for ibin in range(1,nbins+1):
                    outname = "pull_" + name + "_bin" + str(ibin)
                    unfolded_bin_content = histo.GetBinContent(ibin)
                    unfolded_bin_error = histo.GetBinError(ibin)
                    true_bin_content = true_distro.GetBinContent(ibin)
                    true_bin_error = true_distro.GetBinError(ibin)
                    total_bin_error = math.sqrt(unfolded_bin_error**2) #???
                    if(total_bin_error != 0):
                        pull = (unfolded_bin_content-true_bin_content)/total_bin_error
                    else:
                        pull = 9999
                    log.debug('unfolded bin content %s +/- %s, true bin content %s, pull %s' % (unfolded_bin_content, unfolded_bin_error, true_bin_content, pull))
                    pulls_lists[outname].append(pull)
                    outname = "delta_" + name + "_bin" + str(ibin)
                    delta = unfolded_bin_content-true_bin_content
                    log.debug('unfolded bin content %s +/- %s, true bin content %s, delta %s' % (unfolded_bin_content, unfolded_bin_error, true_bin_content, delta))
                    deltas_lists[outname].append(delta)
                    outname = "unfolded_" + name + "_bin" + str(ibin)
                    unfoldeds_lists[outname].append(unfolded_bin_content)
                    unfolded_sigmas_lists[outname].append(unfolded_bin_error)
            
            nneg_bins_hists = [i for i in inputfile.keys() if i.GetName().startswith("nneg_bins")]
            nneg_bins_hists = [asrootpy(i.ReadObj()) for i in nneg_bins_hists]
            for histo in nneg_bins_hists:
                #create pull/delta containers during first iteration
                name = histo.name
                nbins = histo.nbins()
                log.debug("name = %s, n bins = %s" % (name,nbins))
                if not lists_created:
                    outname = name
                    nneg_bins_lists[outname] = []
                outname = name
                nneg_bins_lists[outname].append(histo.GetBinContent(1))
            
            pull_sums_hists = [i for i in inputfile.keys() if i.GetName().startswith("sum_of_pulls")]
            pull_sums_hists = [asrootpy(i.ReadObj()) for i in pull_sums_hists]
            for histo in pull_sums_hists:
                #create pull/delta containers during first iteration
                name = histo.name
                nbins = histo.nbins()
                log.debug("name = %s, n bins = %s" % (name,nbins))
                if not lists_created:
                    outname = name
                    pull_sums_lists[outname] = []
                outname = name
                pull_sums_lists[outname].append(histo.GetBinContent(1))
 
            ratio_sums_hists = [i for i in inputfile.keys() if i.GetName().startswith("sum_of_ratios")]
            ratio_sums_hists = [asrootpy(i.ReadObj()) for i in ratio_sums_hists]
            for histo in ratio_sums_hists:
                #create ratio/delta containers during first iteration
                name = histo.name
                nbins = histo.nbins()
                log.debug("name = %s, n bins = %s" % (name,nbins))
                if not lists_created:
                    outname = name
                    ratio_sums_lists[outname] = []
                outname = name
                ratio_sums_lists[outname].append(histo.GetBinContent(1))
 
            #after the first iteration on the file all the lists are created
            lists_created=True

        os.chdir("..")
       
    #create histograms
    #histo containers
    taus = {}
    for name, vals in taus_lists.iteritems():
        ROOT.TH1.AddDirectory(False) #repeat, you never know
        val_min = min(vals)
        val_min = 0.8*val_min if val_min > 0 else 1.2*val_min
        val_max = max(vals)
        val_max = 0.8*val_max if val_max < 0 else 1.2*val_max
        if val_min == val_max:
            if tau_nbins % 2: #if odd
                val_min, val_max = val_min-0.01, val_min+0.01
            else:
                brange = 0.02
                bwidth = brange/tau_nbins
                val_min, val_max = val_min-0.01+bwidth/2., val_min+0.01+bwidth/2.
        title = '#tau choice - %s ;#tau;N_{toys}' % (name)
        histo = Hist(tau_nbins, val_min, val_max, name=name, title=title)
        for val in vals:
            histo.Fill(val)
        taus[name] = histo        

    pulls = {}
    for name, vals in pulls_lists.iteritems():
        ROOT.TH1.AddDirectory(False) #repeat, you never know
        val_min = min(vals)
        val_min = 0.8*val_min if val_min > 0 else 1.2*val_min
        val_max = max(vals)
        val_max = 0.8*val_max if val_max < 0 else 1.2*val_max
        abs_max = max(abs(val_min),abs(val_max))
        if 'L_curve' in name:
            method = 'L_curve'
            binno = name.split('_')[-1]
        else:
            _, method, binno = tuple(name.split('_'))
        title = 'Pulls - %s - %s ;Pull;N_{toys}' % (binno, method)
        histo = Hist(pull_nbins, -abs_max, abs_max, name=name, title=title)
        for val in vals:
            histo.Fill(val)
        pulls[name] = histo

    deltas = {}
    for name, vals in deltas_lists.iteritems():
        ROOT.TH1.AddDirectory(False) #repeat, you never know
        val_min = min(vals)
        val_min = 0.8*val_min if val_min > 0 else 1.2*val_min
        val_max = max(vals)
        val_max = 0.8*val_max if val_max < 0 else 1.2*val_max
        if 'L_curve' in name:
            method = 'L_curve'
            binno = name.split('_')[-1]
        else:
            _, method, binno = tuple(name.split('_'))
        title = 'Deltas - %s - %s ;Delta;N_{toys}' % (binno, method)
        histo = Hist(delta_nbins, val_min, val_max, name=name, title=title)
        for val in vals:
            histo.Fill(val)
        deltas[name] = histo

    unfoldeds = {}
    for name, vals in unfoldeds_lists.iteritems():
        ROOT.TH1.AddDirectory(False) #repeat, you never know     
        val_min = min(vals)
        val_min = 0.8*val_min if val_min > 0 else 1.2*val_min
        val_max = max(vals)
        val_max = 0.8*val_max if val_max < 0 else 1.2*val_max
        if 'L_curve' in name:
            method = 'L_curve'
            binno = name.split('_')[-1]
        else:
            _, method, binno = tuple(name.split('_'))
        title = 'Unfoldeds - %s - %s ;Unfolded;N_{toys}' % (binno, method)
        histo = Hist(unfolded_nbins, val_min, val_max, name=name, title=title)
        for val in vals:
            histo.Fill(val)
        unfoldeds[name] = histo
        
    nneg_bins = {}
    for name, vals, in nneg_bins_lists.iteritems():
        ROOT.TH1.AddDirectory(False) #repeat, you never know     
        val_min = min(vals)
        val_min = 0 if val_min > 0 else val_min-1
        val_max = max(vals)
        val_max = 0 if val_max < 0 else val_max+1
        if 'L_curve' in name:
            method = 'L_curve'
        else:
            set_trace()
            _, method, _ = tuple(name.split('_'))
        title = 'N of negative bins - %s ;N. neg bins;N_{toys}' % method
        histo = Hist(int(val_max-val_min+1), val_min, val_max, name=name, title=title)
        for val in vals:
            histo.Fill(val)
        nneg_bins[name] = histo
        
    pull_sums = {}
    for name, vals in pull_sums_lists.iteritems():
        ROOT.TH1.AddDirectory(False) #repeat, you never know     
        val_min = min(vals)
        val_min = 0.8*val_min if val_min > 0 else 1.2*val_min
        val_max = max(vals)
        val_max = 0.8*val_max if val_max < 0 else 1.2*val_max
        if 'L_curve' in name:
            method = 'L_curve'
        else:
            set_trace()
            _, _, _, _, _, method = tuple(name.split('_'))
        title = 'Pull sums - %s ;#Sigma(pull)/N_{bins};N_{toys}' % method
        histo = Hist(unfolded_nbins, val_min, val_max, name=name, title=title)
        for val in vals:
            histo.Fill(val)
        pull_sums[name] = histo

    ratio_sums = {}
    for name, vals in ratio_sums_lists.iteritems():
        ROOT.TH1.AddDirectory(False) #repeat, you never know     
        val_min = min(vals)
        val_min = 0.8*val_min if val_min > 0 else 1.2*val_min
        val_max = max(vals)
        val_max = 0.8*val_max if val_max < 0 else 1.2*val_max
        if 'L_curve' in name:
            method = 'L_curve'
            binno = name.split('_')[-1]
        else:
            set_trace()
            _, _, _, _, _, method = tuple(name.split('_'))
        title = 'Ratio sums - %s;#Sigma(ratio)/N_{bins};N_{toys}' % method
        histo = Hist(unfolded_nbins, val_min, val_max, name=name, title=title)
        for val in vals:
            histo.Fill(val)
        ratio_sums[name] = histo

    unfolded_sigmas = {}
    for name, vals in unfolded_sigmas_lists.iteritems():
        ROOT.TH1.AddDirectory(False) #repeat, you never know     
        val_min = min(vals)
        val_min = 0.8*val_min if val_min > 0 else 1.2*val_min
        val_max = max(vals)
        val_max = 0.8*val_max if val_max < 0 else 1.2*val_max
        if 'L_curve' in name:
            method = 'L_curve'
            binno = name.split('_')[-1]
        else:
            _, method, binno = tuple(name.split('_'))
        title = 'Unfolded uncertainties - %s - %s ;Uncertainty;N_{toys}' % (binno, method)
        histo = Hist(unfolded_nbins, val_min, val_max, name=name, title=title)
        for val in vals:
            histo.Fill(val)
        unfolded_sigmas[name] = histo

    for name, histo in pulls.iteritems():
        log.debug("name is %s and object type is %s" % (name, type(histo)))
        histo.Fit("gaus",'Q')
        if not histo.GetFunction("gaus"):
            log.warning("Function not found for histogram %s" % name)
            continue
        mean = histo.GetFunction("gaus").GetParameter(1)
        meanError = histo.GetFunction("gaus").GetParError(1)
        sigma = histo.GetFunction("gaus").GetParameter(2)
        sigmaError = histo.GetFunction("gaus").GetParError(2)

        general_name, idx = tuple(name.split('_bin'))
        idx = int(idx)

        pull_means_lists[general_name][idx] = mean
        pull_mean_errors_lists[general_name][idx] = meanError
        pull_sigmas_lists[general_name][idx] = sigma
        pull_sigma_errors_lists[general_name][idx] = sigmaError
                
    for name, histo in deltas.iteritems():
        log.debug("name is %s and object type is %s" % (name, type(histo)))
        histo.Fit("gaus",'Q')
        if not histo.GetFunction("gaus"):
            log.warning("Function not found for histogram %s" % name)
            continue
        mean = histo.GetFunction("gaus").GetParameter(1)
        meanError = histo.GetFunction("gaus").GetParError(1)
        sigma = histo.GetFunction("gaus").GetParameter(2)
        sigmaError = histo.GetFunction("gaus").GetParError(2)

        general_name, idx = tuple(name.split('_bin'))
        idx = int(idx)

        delta_means_lists[general_name][idx] = mean
        delta_mean_errors_lists[general_name][idx] = meanError
        delta_sigmas_lists[general_name][idx] = sigma
        delta_sigma_errors_lists[general_name][idx] = sigmaError
        
    
    outfile = rootpy.io.File("unfolding_diagnostics.root", "RECREATE")
    outfile.cd()

    pull_means = {}
    pull_sigmas = {}
    pull_means_summary = {}
    pull_sigmas_summary = {}
    delta_means = {}
    delta_sigmas = {}
    delta_means_summary = {}
    delta_sigmas_summary = {}
    
    for outname, pmeans in pull_means_lists.iteritems():
        outname_mean = outname + "_mean"
        outtitle = "Pull means - " + outname + ";Pull mean; N_{toys}"
        pull_mean_min = min(pmeans.values())
        pull_mean_max = max(pmeans.values())
        pull_mean_newmin = pull_mean_min - (pull_mean_max-pull_mean_min)*0.5
        pull_mean_newmax = pull_mean_max + (pull_mean_max-pull_mean_min)*0.5
        pull_means[outname] = plotting.Hist(pull_mean_nbins, pull_mean_newmin, pull_mean_newmax, name = outname_mean, title=outtitle)

        outname_mean_summary = outname + "_mean_summary"
        outtitle_mean_summary = "Pull mean summary - " + outname
        histocloned = true_distro.Clone(outname_mean_summary)
        histocloned.Reset()
        histocloned.xaxis.title = xaxislabel
        histocloned.yaxis.title = 'Pull mean'
        histocloned.title = outtitle_mean_summary
        pull_means_summary[outname] = histocloned

        for idx, pmean in pmeans.iteritems():
            pull_means[outname].Fill(pmean)
            histocloned[idx].value = pmean
            histocloned[idx].error = pull_mean_errors_lists[outname][idx]
        histocloned.yaxis.SetRangeUser(min(pmeans.values()), max(pmeans.values()))

    
    for outname, psigmas in pull_sigmas_lists.iteritems():
        outname_sigma = outname + "_sigma"
        outtitle_sigma = "Pull #sigma's - " + outname + ";Pull #sigma; N_{toys}"
        pull_sigma_min = min(psigmas.values())
        pull_sigma_max = max(psigmas.values())
        pull_sigma_newmin = pull_sigma_min - (pull_sigma_max-pull_sigma_min)*0.5
        pull_sigma_newmax = pull_sigma_max + (pull_sigma_max-pull_sigma_min)*0.5
        pull_sigmas[outname] = plotting.Hist(pull_sigma_nbins, pull_sigma_newmin, pull_sigma_newmax, name = outname_sigma, title=outtitle_sigma)

        outname_sigma_summary = outname + "_sigma_summary"
        outtitle_sigma_summary = "Pull #sigma summary - " + outname
        histocloned = true_distro.Clone(outname_sigma_summary)
        histocloned.Reset()
        histocloned.xaxis.title = xaxislabel
        histocloned.yaxis.title = 'Pull #sigma'
        histocloned.title = outtitle_sigma_summary
        pull_sigmas_summary[outname] = histocloned

        for idx, psigma in psigmas.iteritems():
            pull_sigmas[outname].Fill(psigma)
            histocloned[idx].value = psigma
            histocloned[idx].error = pull_sigma_errors_lists[outname][idx]
        histocloned.yaxis.SetRangeUser(min(psigmas.values()), max(psigmas.values()))

    for outname, dmeans in delta_means_lists.iteritems():
        outname_mean = outname + "_mean"
        outtitle = "Delta means - " + outname + ";Delta mean; N_{toys}"
        delta_mean_min = min(dmeans.values())
        delta_mean_max = max(dmeans.values())
        delta_mean_newmin = delta_mean_min - (delta_mean_max-delta_mean_min)*0.5
        delta_mean_newmax = delta_mean_max + (delta_mean_max-delta_mean_min)*0.5
        delta_means[outname] = plotting.Hist(delta_mean_nbins, delta_mean_newmin, delta_mean_newmax, name = outname_mean, title=outtitle)

        outname_mean_summary = outname + "_mean_summary"
        outtitle_mean_summary = "Delta mean summary - " + outname
        histocloned = true_distro.Clone(outname_mean_summary)
        histocloned.Reset()
        histocloned.xaxis.title = xaxislabel
        histocloned.yaxis.title = 'Delta mean'
        histocloned.title = outtitle_mean_summary
        delta_means_summary[outname] = histocloned

        for idx, dmean in dmeans.iteritems():
            delta_means[outname].Fill(dmean)
            histocloned[idx].value = dmean
            histocloned[idx].error = delta_mean_errors_lists[outname][idx]
        histocloned.yaxis.SetRangeUser(min(dmeans.values()), max(dmeans.values()))
        
    for outname, dsigmas in delta_sigmas_lists.iteritems():
        outname_sigma = outname + "_sigma"
        outtitle_sigma = "Delta #sigma's - " + outname + ";Delta #sigma; N_{toys}"
        delta_sigma_min = min(dsigmas.values())
        delta_sigma_max = max(dsigmas.values())
        delta_sigma_newmin = delta_sigma_min - (delta_sigma_max-delta_sigma_min)*0.5
        delta_sigma_newmax = delta_sigma_max + (delta_sigma_max-delta_sigma_min)*0.5
        delta_sigmas[outname] = plotting.Hist(delta_sigma_nbins, delta_sigma_newmin, delta_sigma_newmax, name = outname_sigma, title=outtitle_sigma)

        outname_sigma_summary = outname + "_sigma_summary"
        outtitle_sigma_summary = "Delta #sigma summary - " + outname
        histocloned = true_distro.Clone(outname_sigma_summary)
        histocloned.Reset()
        histocloned.xaxis.title = xaxislabel
        histocloned.yaxis.title = 'Delta #sigma'
        histocloned.title = outtitle_sigma_summary
        delta_sigmas_summary[outname] = histocloned

        for idx, dsigma in dsigmas.iteritems():
            delta_sigmas[outname].Fill(dsigma)
            histocloned[idx].value = dsigma
            histocloned[idx].error = delta_sigma_errors_lists[outname][idx]
        histocloned.yaxis.SetRangeUser(min(dsigmas.values()), max(dsigmas.values()))

    unfolded_summary = {}
    unfolded_average = {}
    unfolded_envelope = {}
    for name, histo in unfoldeds.iteritems():
        log.debug("name is %s and object type is %s" % (name, type(histo)))
        histo.Fit("gaus",'Q')
        if not histo.GetFunction("gaus"):
            log.warning("Function not found for histogram %s" % name)
            continue
        mean = histo.GetFunction("gaus").GetParameter(1)
        meanError = histo.GetFunction("gaus").GetParError(1)
        sigma = histo.GetFunction("gaus").GetParameter(2)
        sigmaError = histo.GetFunction("gaus").GetParError(2)

        general_name, idx = tuple(name.split('_bin'))
        idx = int(idx)

        if general_name not in unfolded_summary:
            histo = true_distro.Clone("%s_unfolded_summary" % general_name)
            outtitle_unfolded_summary = "Unfolded summary - " + general_name
            histo.Reset()
            histo.xaxis.title = xaxislabel
            histo.yaxis.title = 'N_{events}'
            histo.title = outtitle_unfolded_summary
            unfolded_summary[general_name] = histo

            unfolded_envelope[general_name] = histo.Clone("%s_unfolded_envelope" % general_name)
            unfolded_average[general_name] = histo.Clone("%s_unfolded_average" % general_name)

        unfolded_summary[general_name][idx].value = mean
        unfolded_summary[general_name][idx].error = meanError

        unfolded_envelope[general_name][idx].value = mean
        unfolded_envelope[general_name][idx].error = sigma
        
        unfolded_average[general_name][idx].value = mean
        unfolded_average[general_name][idx].error = \
           unfolded_sigmas['%s_bin%i' % (general_name, idx)].GetMean()
    
    plotter.set_subdir('taus')
    for name, histo in taus.iteritems():
        #canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        plotter.canvas.cd()
        histo = plotter.plot(
            histo, **styles['dots']
            )
        histo.SetStats(True)

        info = plotter.make_text_box(
            'mode #tau = %.5f' % histo[histo.GetMaximumBin()].x.center,
            position=(plotter.pad.GetLeftMargin(), plotter.pad.GetTopMargin(), 0.3, 0.025)
            )
        info.Draw()
        
        plotter.save()
        histo.Write()
        plotter.canvas.Write()


    plotter.set_subdir('pulls')
    for name, histo in pulls.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.SetStats(True)
        plotter.save()
        histo.Write()
        plotter.canvas.Write()
    for name, histo in pull_means.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.Write()
        plotter.save()
    for name, histo in pull_sigmas.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.Write()
        plotter.save()

    plotter.set_subdir('pull_summaries')
    for name, histo in pull_means_summary.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        #histo.SetStats(True)
        line = ROOT.TLine(histo.GetBinLowEdge(1),0,histo.GetBinLowEdge(histo.GetNbinsX()+1),0)
        line.Draw("same")
        plotter.save()
        histo.Write()
        plotter.canvas.Write()
    for name, histo in pull_sigmas_summary.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        #histo.SetStats(True)
        line = ROOT.TLine(histo.GetBinLowEdge(1),1,histo.GetBinLowEdge(histo.GetNbinsX()+1),1)
        line.Draw("same")
        plotter.save()
        histo.Write()
        plotter.canvas.Write()
        
    plotter.set_subdir('deltas')
    for name, histo in deltas.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.SetStats(True)
        plotter.save()
        histo.Write()
        plotter.canvas.Write()
    for name, histo in delta_means.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.Write()
        plotter.save()
    for name, histo in delta_sigmas.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.Write()
        plotter.save()

    plotter.set_subdir('delta_summaries')
    for name, histo in delta_means_summary.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        #histo.SetStats(True)
        plotter.save()
        histo.Write()
        plotter.canvas.Write()
    for name, histo in delta_sigmas_summary.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        #histo.SetStats(True)
        plotter.save()
        histo.Write()
        plotter.canvas.Write()

    plotter.set_subdir('unfolding_unc')
    for name, histo in unfolded_sigmas.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.SetStats(True)
        plotter.save()
        histo.Write()
        plotter.canvas.Write()
        
    plotter.set_subdir('unfolded')
    for name, histo in unfoldeds.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.SetStats(True)
        plotter.save()
        histo.Write()
        plotter.canvas.Write()

    plotter.set_subdir('unfolded_summaries')
    for name, histo in unfolded_summary.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.SetStats(True)
        plotter.save()
        histo.Write()
        plotter.canvas.Write()

    for name, histo in unfolded_summary.iteritems():
        leg = LegendDefinition("Unfolding comparison", 'NE', labels=['Truth', 'Unfolded'])
        plotter.overlay_and_compare([true_distro], histo, legend_def=leg, **styles['compare'])
        plotter.canvas.name = 'Pull_'+name
        plotter.save()
        plotter.canvas.Write()
        plotter.overlay_and_compare(
            [true_distro], histo, legend_def=leg, method='ratio', **styles['compare']
            )
        plotter.canvas.name = 'Ratio_'+name
        plotter.save()
        plotter.canvas.Write()

    plotter.set_subdir('unfolded_average')
    for name, histo in unfolded_average.iteritems():
        leg = LegendDefinition("Unfolding comparison", 'NE', labels=['Truth', 'Unfolded'])
        #set_trace()
        plotter.overlay_and_compare([true_distro], histo, legend_def=leg, **styles['compare'])
        plotter.canvas.name = 'Pull_'+name
        plotter.save()
        plotter.canvas.Write()
        plotter.overlay_and_compare([true_distro], histo, legend_def=leg, method='ratio', **styles['compare'])
        plotter.canvas.name = 'Ratio_'+name
        plotter.save()
        plotter.canvas.Write()

    plotter.set_subdir('unfolded_envelope')
    for name, histo in unfolded_envelope.iteritems():
        leg = LegendDefinition("Unfolding comparison", 'NE', labels=['Truth', 'Unfolded'])
        plotter.overlay_and_compare([true_distro], histo, legend_def=leg, **styles['compare'])
        plotter.canvas.name = 'Pull_'+name
        plotter.save()
        plotter.canvas.Write()
        plotter.overlay_and_compare([true_distro], histo, legend_def=leg, method='ratio', **styles['compare'])
        plotter.canvas.name = 'Ratio_'+name
        plotter.save()
        plotter.canvas.Write()
        
    plotter.set_subdir('figures_of_merit')
    for name, histo in nneg_bins.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.SetStats(True)
        plotter.save()
        histo.Write()
        plotter.canvas.Write()
    for name, histo in pull_sums.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.SetStats(True)
        plotter.save()
        histo.Write()
        plotter.canvas.Write()
    for name, histo in ratio_sums.iteritems():
        histo = plotter.plot(histo, **styles['dots'])
        histo.SetStats(True)
        plotter.save()
        histo.Write()
        plotter.canvas.Write()


    outfile.close()
    os.chdir(curdir)

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('input_dir', type=str, help='input directory')
    parser.add_argument('variable', type=str, help='unfolding variable')
    opts = parser.parse_args()
    unfolding_toy_diagnostics(opts.input_dir, opts.variable)
    
