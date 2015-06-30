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


def set_pretty_label(variable):
    if variable == 'ptthad':
        result = 'p_{T}(t_{had}) [GeV]'
    elif variable == 'pttlep':
        result = 'p_{T}(t_{lep}) [GeV]'
    elif variable == 'etathad':
        result == '|#eta(t_{had})|'
    elif variable == 'etatlep':
        result == '|#eta(t_{lep})|'
    else:
        result == variable
    return result


def unfolding_toy_diagnostics(indir, variable):
    
    plotter = BasePlotter()
    
    xaxislabel = set_pretty_label(variable)
    
    true_distribution = None
    

    curdir = os.getcwd()
    os.chdir(indir)
    toydirs = get_immediate_subdirectories(".")
    
    pulls_lists = {}
    pull_means_lists = {}
    pull_mean_errors_lists = {}
    pull_sigmas_lists = {}
    pull_sigma_errors_lists = {}
    deltas_lists = {}
    delta_means_lists = {}
    delta_mean_errors_lists = {}
    delta_sigmas_lists = {}
    delta_sigma_errors_lists = {}
    unfoldeds_lists = {}
    
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
        with io.root_open("result_unfolding.root") as inputfile:
            log.debug('Iteration %s over the file' % i)
            i = i+1
            keys = inputfile.keys()
            unfolded_hists = [i for i in inputfile.keys() if i.GetName().startswith("hdata_unfolded")]
            unfolded_hists = [asrootpy(i.ReadObj()) for i in unfolded_hists]
            true_distribution = inputfile.true_distribution
            ROOT.TH1.AddDirectory(False)
            true_distro = true_distribution.Clone()

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
                    true_bin_content = true_distribution.GetBinContent(ibin)
                    true_bin_error = true_distribution.GetBinError(ibin)
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
            #after the first iteration on the file all the lists are created
            lists_created=True
        os.chdir("..")
       
    #create histograms
    #histo containers
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
            _, _, _, method, binno = tuple(name.split('_'))
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
            _, _, _, method, binno = tuple(name.split('_'))
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
            _, _, _, method, binno = tuple(name.split('_'))
        title = 'Unfoldeds - %s - %s ;Unfolded;N_{toys}' % (binno, method)
        histo = Hist(unfolded_nbins, val_min, val_max, name=name, title=title)
        for val in vals:
            histo.Fill(val)
        unfoldeds[name] = histo

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
    unfolded_summary = {}
    
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

        unfolded_summary[general_name][idx].value = mean
        unfolded_summary[general_name][idx].error = meanError
    
    subdir = 'pulls'
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    for name, histo in pulls.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))
    for name, histo in pull_means.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo)
        histo.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))
    for name, histo in pull_sigmas.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo)
        histo.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))


    subdir = 'pull_summaries'
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    for name, histo in pull_means_summary.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        #histo.SetStats(True)
        line = ROOT.TLine(histo.GetBinLowEdge(1),0,histo.GetBinLowEdge(histo.GetNbinsX()+1),0)
        line.Draw("same")
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))
    for name, histo in pull_sigmas_summary.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        #histo.SetStats(True)
        line = ROOT.TLine(histo.GetBinLowEdge(1),1,histo.GetBinLowEdge(histo.GetNbinsX()+1),1)
        line.Draw("same")
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))
        
    subdir = 'deltas'
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    for name, histo in deltas.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))
    for name, histo in delta_means.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo)
        histo.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))
    for name, histo in delta_sigmas.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo)
        histo.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))

    subdir = 'delta_summaries'
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    for name, histo in delta_means_summary.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        #histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))
    for name, histo in delta_sigmas_summary.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        #histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))

    subdir = 'unfolded'
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    for name, histo in unfoldeds.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))

    subdir = 'unfolded_summaries'
    if not os.path.isdir(subdir):
        os.makedirs(subdir)
    for name, histo in unfolded_summary.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s/%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s/%s.pdf' % (subdir, canvas.GetName()))
    for name, histo in unfolded_summary.iteritems():
        leg = LegendDefinition("Unfolding comparison", ['Truth', 'Unfolded'], 'NE')
        canvas = plotter.create_and_write_canvas_with_comparison('Pull_'+name,[1,0],[0,21],[2,1], leg, False, False, [true_distribution, histo], write=False, comparison = 'pull')
        canvas.Update()
        canvas.Write()
        canvas.SaveAs('%s%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s%s.pdf' % (subdir, canvas.GetName()))
        canvas = plotter.create_and_write_canvas_with_comparison('Ratio_'+name,[1,0],[0,21],[2,1], leg, False, False, [true_distribution, histo], write=False, comparison = 'ratio')
        canvas.Update()
        canvas.Write()
        canvas.SaveAs('%s%s.png' % (subdir, canvas.GetName()))
        canvas.SaveAs('%s%s.pdf' % (subdir, canvas.GetName()))
        
    outfile.close()
    os.chdir(curdir)

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('input_dir', type=str, help='input directory')
    parser.add_argument('variable', type=str, help='unfolding variable')
    opts = parser.parse_args()
    unfolding_toy_diagnostics(opts.input_dir, opts.variable)
    
