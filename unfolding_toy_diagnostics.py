#! /bin/env python

import math
import ROOT
import rootpy.plotting as plotting
import rootpy
import rootpy.io as io
from URAnalysis.PlotTools.views import EnvelopeView
from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
from URAnalysis.Utilities.roottools import ArgSet, ArgList
from pdb import set_trace
import logging
from os.path import join
import os

asrootpy = rootpy.asrootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(11111)
ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)
log = rootpy.log["/toy_diagnostics"]

pull_nbins = 100

pull_mean_nbins = 100
pull_mean_min = -10
pull_mean_max = 10

pull_sigma_nbins = 100

delta_nbins = 100

delta_mean_nbins = 100

delta_sigma_nbins = 100

unfolded_nbins = 100



import URAnalysis.Utilities.prettyjson as prettyjson

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
    pulls = {}
    pull_means = {}
    pull_sigmas = {}
    pull_means_summary = {}
    pull_sigmas_summary = {}
    deltas = {}
    delta_means = {}
    delta_sigmas = {}
    delta_means_summary = {}
    delta_sigmas_summary = {}
    unfoldeds = {}
    unfolded_summary = {}
    
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
    for directory in toydirs:
        os.chdir(directory)
        print('Iteration %s over the directories' % idir)
        idir = idir + 1
        i = 0

        with io.root_open("result_unfolding.root") as inputfile:
            print ('Iteration %s over the file' % i)
            i = i+1
            keys = inputfile.keys()
            
            
            if not lists_created:
                for key in keys:
                    if key.GetName().startswith("hdata_unfolded"):
                        histo = asrootpy(key.ReadObj())
                        name = histo.GetName()
                        nbins = histo.GetXaxis().GetNbins()
                        print ("name = %s, n bins = %s" % (name,nbins))
                        for ibin in range(1,nbins+1):
                            outname = "pull_" + name + "_bin" + str(ibin)
                            pulls_lists[outname] = []
                            outname = "delta_" + name + "_bin" + str(ibin)
                            deltas_lists[outname] = []
                            outname = "unfolded_" + name + "_bin" + str(ibin)
                            unfoldeds_lists[outname] = []
                        outname = "pull_" + name
                        pull_means_lists[outname] = []
                        pull_mean_errors_lists[outname] = []
                        pull_sigmas_lists[outname] = []
                        pull_sigma_errors_lists[outname] = []
                                                
                        outname = "delta_" + name
                        delta_means_lists[outname] = []
                        delta_mean_errors_lists[outname] = []
                        delta_sigmas_lists[outname] = []
                        delta_sigma_errors_lists[outname] = []
                        
                        outname = 'unfolded_' + name
                                            
                lists_created = True
            
            
            
            #true_distribution = keys.FindObject("true_distribution").ReadObj()
            true_distribution = filter(lambda x: x.GetName() == 'true_distribution', keys)[0].ReadObj()
            j = 0
            for key in keys:
                print('Iteration %s over the keys' % j)
                j = j+1
                if key.GetName().startswith("hdata_unfolded"):
                    histo = key.ReadObj()
                    name = histo.GetName()
                    nbins = histo.GetXaxis().GetNbins()
                    print ("name = %s, n bins = %s" % (name,nbins))
                    for ibin in range(1,nbins+1):
                        outname = "pull_" + name + "_bin" + str(ibin)
                        unfolded_bin_content = histo.GetBinContent(ibin)
                        unfolded_bin_error = histo.GetBinError(ibin)
                        true_bin_content = true_distribution.GetBinContent(ibin)
                        true_bin_error = true_distribution.GetBinError(ibin)
                        #total_bin_error = math.sqrt(unfolded_bin_error**2 + true_bin_error**2)
                        total_bin_error = math.sqrt(unfolded_bin_error**2)
                        if(total_bin_error != 0):
                            pull = (unfolded_bin_content-true_bin_content)/total_bin_error
                        else:
                            pull = 9999
                        print ('unfolded bin content %s +/- %s, true bin content %s, pull %s' % (unfolded_bin_content, unfolded_bin_error, true_bin_content, pull))
                        #print outname
                        #print pulls
                        pulls_lists[outname].append(pull)
                        outname = "delta_" + name + "_bin" + str(ibin)
                        delta = unfolded_bin_content-true_bin_content
                        print ('unfolded bin content %s +/- %s, true bin content %s, delta %s' % (unfolded_bin_content, unfolded_bin_error, true_bin_content, delta))
                        #print outname
                        #print deltas
                        deltas_lists[outname].append(delta)
                        outname = "unfolded_" + name + "_bin" + str(ibin)
                        unfoldeds_lists[outname].append(unfolded_bin_content)
                        #print outname
                        #print unfoldeds
                        
            

        os.chdir("..")
       
    histos_created = False
    for directory in toydirs:
        os.chdir(directory)
        print('Iteration %s over the directories' % idir)
        idir = idir + 1
        i = 0

        with io.root_open("result_unfolding.root") as inputfile:
            print ('Iteration %s over the file' % i)
            i = i+1
            keys = inputfile.keys()
        
        if not histos_created:
            for key in keys:
                if key.GetName().startswith("hdata_unfolded"):
                    histo = asrootpy(key.ReadObj())
                    name = histo.GetName()
                    nbinpullss = histo.GetXaxis().GetNbins()
                    print ("name = %s, n bins = %s" % (name,nbins))
                    for ibin in range(1,nbins+1):
                        ROOT.TH1.AddDirectory(False)
                        outname = "pull_" + name + "_bin" + str(ibin)
                        outtitle = 'Pulls - bin ' + str(ibin) + ' - ' + name + ';Pull;N_{toys}'
                        pull_minmax = max(abs(min(pulls_lists[outname]))*1.2,abs(max(pulls_lists[outname]))*1.2)
                        pulls[outname] = plotting.Hist(pull_nbins, -pull_minmax, pull_minmax, name = outname, title=outtitle)
                        print ("outname is %s, we are in bin %s" % (outname,ibin))
                        outname = "delta_" + name + "_bin" + str(ibin)
                        outtitle = 'Deltas - bin ' + str(ibin) + ' - ' + name + ';Delta;N_{toys}'
                        delta_minmax = max(abs(min(deltas_lists[outname]))*1.2,abs(max(deltas_lists[outname]))*1.2)
                        deltas[outname] = plotting.Hist(delta_nbins, -delta_minmax, delta_minmax, name = outname, title=outtitle)
                        print ("outname is %s, we are in bin %s" % (outname,ibin))
                        outname = "unfolded_" + name + "_bin" + str(ibin)
                        outtitle = 'Unfoldeds - bin ' + str(ibin) + ' - ' + name + ';Unfolded;N_{toys}'
                        unfoldeds_minmax = max(abs(min(unfoldeds_lists[outname]))*1.2,abs(max(unfoldeds_lists[outname]))*1.2)
                        unfoldeds[outname] = plotting.Hist(unfolded_nbins, -unfolded_minmax, unfolded_minmax, name = outname, title=outtitle)
                        print ("outname is %s, we are in bin %s" % (outname,ibin))
    
                    
            histos_created = True
        
        if histos_created:
            break
        
    for name, pull in pulls_lists:
        pulls[name].Fill(pull)
    
    for name, delta in deltas_lists:
        deltas[name].Fill(delta)
    
    for name, unfolded in unfoldeds_lists:
        unfoldeds[name].Fill(unfolded)    
    
    outfile = rootpy.io.File("unfolding_diagnostics.root", "RECREATE")
    
    outfile.cd()

    for name, histo in pulls.iteritems():
        #set_trace()
        print("name is %s and object type is %s" % (name, type(histo)))
        histo.Fit("gaus")
        if not histo.GetFunction("gaus"):
            log.warning("Function not found for histogram %s" % name)
            continue
        mean = histo.GetFunction("gaus").GetParameter(1)
        meanError = histo.GetFunction("gaus").GetParError(1)
        sigma = histo.GetFunction("gaus").GetParameter(2)
        sigmaError = histo.GetFunction("gaus").GetParError(2)
        if filter(str.isdigit, name) != '':
            index = int(filter(str.isdigit, name))
        else:
            index = 1
        print index
        for name_mean, list_mean in pull_means_lists.iteritems():
            if name.startswith(name_mean):
                list_mean.append(mean)
                pull_mean_errors_lists[name_mean].append(meanError)
        for name_sigma, list_sigma in pull_sigmas_lists.iteritems():
            if name.startswith(name_sigma):
                list_sigma.append(sigma)
                pull_sigma_errors_lists[name_sigma].append(sigmaError)
                
    for name, histo in deltas.iteritems():
        #set_trace()
        print("name is %s and object type is %s" % (name, type(histo)))
        histo.Fit("gaus")
        if not histo.GetFunction("gaus"):
            log.warning("Function not found for histogram %s" % name)
            continue
        mean = histo.GetFunction("gaus").GetParameter(1)
        meanError = histo.GetFunction("gaus").GetParError(1)
        sigma = histo.GetFunction("gaus").GetParameter(2)
        sigmaError = histo.GetFunction("gaus").GetParError(2)
        if filter(str.isdigit, name) != '':
            index = int(filter(str.isdigit, name))
        else:
            index = 1
        print index
        for name_mean, list_mean in delta_means_lists.iteritems():
            if name.startswith(name_mean):
                list_mean.append(mean)
                delta_mean_errors_lists[name_mean].append(meanError)
        for name_sigma, list_sigma in delta_sigmas_lists.iteritems():
            if name.startswith(name_sigma):
                list_sigma.append(sigma)
                delta_sigma_errors_lists[name_sigma].append(sigmaError)
                
    
    histos_created = False
    for directory in toydirs:
        os.chdir(directory)
        print('Iteration %s over the directories' % idir)
        idir = idir + 1
        i = 0

        with io.root_open("result_unfolding.root") as inputfile:
            print ('Iteration %s over the file' % i)
            i = i+1
            keys = inputfile.keys()
        
        if not histos_created:
            for key in keys:
                if key.GetName().startswith("hdata_unfolded"):
                    histo = asrootpy(key.ReadObj())
                    name = histo.GetName()
                    nbinpullss = histo.GetXaxis().GetNbins()
                    print ("name = %s, n bins = %s" % (name,nbins))

                    outname = "pull_" + name
                    outname_mean = outname + "_mean"
                    outtitle = "Pull means - " + name + ";Pull mean; N_{toys}"
                    pull_mean_min = min(pulls_means_lists[outname])
                    pull_mean_max = max(pulls_means_lists[outname])
                    pull_mean_newmin = pull_mean_min - (pull_mean_max-pull_mean_min)*0.2
                    pull_mean_newmax = pull_mean_max + (pull_mean_max-pull_mean_min)*0.2
                    pull_means[outname] = plotting.Hist(pull_mean_nbins, pull_mean_newmin, pull_mean_newmax, name = outname_mean, title=outtitle)
                    outname_sigma = outname + "_sigma"
                    outtitle_sigma = "Pull #sigma's - " + name + ";Pull #sigma; N_{toys}"
                    pull_sigma_min = min(pulls_sigmas_lists[outname])
                    pull_sigma_max = max(pulls_sigmas_lists[outname])
                    pull_sigma_newmin = pull_sigma_min - (pull_sigma_max-pull_sigma_min)*0.2
                    pull_sigma_newmax = pull_sigma_max + (pull_sigma_max-pull_sigma_min)*0.2
                    pull_sigmas[outname] = plotting.Hist(pull_sigma_nbins, pull_sigma_newmin, pull_sigma_newmax, name = outname_sigma, title=outtitle_sigma)

                    outname_mean_summary = outname + "_mean_summary"
                    outtitle_mean_summary = "Pull mean summary - " + name
                    histocloned = histo.Clone(outname_mean_summary)
                    histocloned.Reset()
                    histocloned.xaxis.title = xaxislabel
                    histocloned.yaxis.title = 'Pull mean'
                    histocloned.title = outtitle_mean_summary
                    pull_means_summary[outname] = histocloned
                    outname_sigma_summary = outname + "_sigma_summary"
                    outtitle_sigma_summary = "Pull #sigma summary - " + name
                    histocloned = histo.Clone(outname_sigma_summary)
                    histocloned.Reset()
                    histocloned.xaxis.title = xaxislabel
                    histocloned.yaxis.title = 'Pull #sigma'
                    histocloned.title = outtitle_sigma_summary
                    pull_sigmas_summary[outname] = histocloned
                    
                    outname = "delta_" + name
                    outname_mean = outname + "_mean"
                    outtitle = "Delta means - " + name + ";Delta mean; N_{toys}"
                    delta_mean_min = min(deltas_means_lists[outname])
                    delta_mean_max = max(deltas_means_lists[outname])
                    delta_mean_newmin = delta_mean_min - (delta_mean_max-delta_mean_min)*0.2
                    delta_mean_newmax = delta_mean_max + (delta_mean_max-delta_mean_min)*0.2
                    delta_means[outname] = plotting.Hist(delta_mean_nbins, delta_mean_newmin, delta_mean_newmax, name = outname_mean, title=outtitle)
                    outname_sigma = outname + "_sigma"
                    outtitle_sigma = "Delta #sigma's - " + name + ";Delta #sigma; N_{toys}"
                    delta_sigma_min = min(deltas_sigmas_lists[outname])
                    delta_sigma_max = max(deltas_sigmas_lists[outname])
                    delta_sigma_newmin = delta_sigma_min - (delta_sigma_max-delta_sigma_min)*0.2
                    delta_sigma_newmax = delta_sigma_max + (delta_sigma_max-delta_sigma_min)*0.2
                    delta_sigmas[outname] = plotting.Hist(delta_sigma_nbins, delta_sigma_newmin, delta_sigma_newmax, name = outname_sigma, title=outtitle_sigma)
                    
                    outname_mean_summary = outname + "_mean_summary"
                    outtitle_mean_summary = "Delta mean summary - " + name
                    histocloned = histo.Clone(outname_mean_summary)
                    histocloned.Reset()
                    histocloned.xaxis.title = xaxislabel
                    histocloned.yaxis.title = 'Delta mean'
                    histocloned.title = outtitle_mean_summary
                    delta_means_summary[outname] = histocloned
                    outname_sigma_summary = outname + "_sigma_summary"
                    outtitle_sigma_summary = "Delta #sigma summary - " + name
                    histocloned = histo.Clone(outname_sigma_summary)
                    histocloned.Reset()
                    histocloned.xaxis.title = xaxislabel
                    histocloned.yaxis.title = 'Delta #sigma'
                    histocloned.title = outtitle_sigma_summary
                    delta_sigmas_summary[outname] = histocloned
                    
                    outname = 'unfolded_' + name
                    outname_unfolded_summary = name + "_unfolded_summary"
                    outtitle_unfolded_summary = "Unfolded summary - " + name
                    histocloned = histo.Clone(outname_unfolded_summary)
                    histocloned.Reset()
                    histocloned.xaxis.title = xaxislabel
                    histocloned.yaxis.title = 'N_{events}'
                    histocloned.title = outtitle_unfolded_summary
                    unfolded_summary[outname] = histocloned
            histos_created = True
        
        if histos_created:
            break  

    for name, value in pull_means_lists.iteritems():
        pull_means[name].Fill(value)
    for name, value in pull_sigmas_lists.iteritems():
        pull_sigmas[name].Fill(value)
    for name, value in delta_means_lists.iteritems():
        delta_means[name].Fill(value)
    for name, value in delta_sigmas_lists.iteritems():
        delta_sigmas[name].Fill(value)
    
    for name, histo in pulls.iteritems():
        print("name is %s and object type is %s" % (name, type(histo)))
        histo.Fit("gaus")
        if not histo.GetFunction("gaus"):
            log.warning("Function not found for histogram %s" % name)
            continue
        mean = histo.GetFunction("gaus").GetParameter(1)
        meanError = histo.GetFunction("gaus").GetParError(1)
        sigma = histo.GetFunction("gaus").GetParameter(2)
        sigmaError = histo.GetFunction("gaus").GetParError(2)
        if filter(str.isdigit, name) != '':
            index = int(filter(str.isdigit, name))
        else:
            index = 1
        print index
        for name_mean_summary, list_mean_summary in pull_means_summary.iteritems():
            if name.startswith(name_mean_summary):
                pull_means_summary[name_mean_summary].SetBinContent(index,mean)
                pull_means_summary[name_mean_summary].SetBinError(index,meanError)
        for name_sigma_summary, histo_sigma_summary in pull_sigmas_summary.iteritems():
            if name.startswith(name_sigma_summary):
                pull_sigmas_summary[name_sigma_summary].SetBinContent(index,sigma)
                pull_sigmas_summary[name_sigma_summary].SetBinError(index,sigmaError)
    
    for name, histo in deltas.iteritems():
        #set_trace()
        print("name is %s and object type is %s" % (name, type(histo)))
        histo.Fit("gaus")
        if not histo.GetFunction("gaus"):
            log.warning("Function not found for histogram %s" % name)
            continue
        mean = histo.GetFunction("gaus").GetParameter(1)
        meanError = histo.GetFunction("gaus").GetParError(1)
        sigma = histo.GetFunction("gaus").GetParameter(2)
        sigmaError = histo.GetFunction("gaus").GetParError(2)
        if filter(str.isdigit, name) != '':
            index = int(filter(str.isdigit, name))
        else:
            index = 1
        print index
        for name_mean_summary, histo_mean_summary in delta_means_summary.iteritems():
            if name.startswith(name_mean_summary):
                delta_means_summary[name_mean_summary].SetBinContent(index,mean)
                delta_means_summary[name_mean_summary].SetBinError(index,meanError)
        for name_sigma_summary, histo_sigma_summary in delta_sigmas_summary.iteritems():
            if name.startswith(name_sigma_summary):
                delta_sigmas_summary[name_sigma_summary].SetBinContent(index,sigma)
                delta_sigmas_summary[name_sigma_summary].SetBinError(index,sigmaError)
    
    for name, histo in unfoldeds.iteritems():
        #set_trace()
        print("name is %s and object type is %s" % (name, type(histo)))
        histo.Fit("gaus")
        if not histo.GetFunction("gaus"):
            log.warning("Function not found for histogram %s" % name)
            continue
        mean = histo.GetFunction("gaus").GetParameter(1)
        meanError = histo.GetFunction("gaus").GetParError(1)
        sigma = histo.GetFunction("gaus").GetParameter(2)
        sigmaError = histo.GetFunction("gaus").GetParError(2)
        if filter(str.isdigit, name) != '':
            index = int(filter(str.isdigit, name))
        else:
            index = 1
        print index
        for name_unfolded_summary, histo_unfolded_summary in unfolded_summary.iteritems():
            #set_trace()
            if name.startswith(name_unfolded_summary):
                unfolded_summary[name_unfolded_summary].SetBinContent(index, mean)
                unfolded_summary[name_unfolded_summary].SetBinError(index, meanError)   
        
    for name, histo in pulls.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
    for name, histo in pull_means.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo)
        histo.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
    for name, histo in pull_sigmas.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo)
        histo.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
    for name, histo in pull_means_summary.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        #histo.SetStats(True)
        line = ROOT.TLine(histo.GetBinLowEdge(1),0,histo.GetBinLowEdge(histo.GetNbinsX()+1),0)
        line.Draw("same")
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
    for name, histo in pull_sigmas_summary.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        #histo.SetStats(True)
        line = ROOT.TLine(histo.GetBinLowEdge(1),1,histo.GetBinLowEdge(histo.GetNbinsX()+1),1)
        line.Draw("same")
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
        

    for name, histo in deltas.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
    for name, histo in delta_means.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo)
        histo.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
    for name, histo in delta_sigmas.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo)
        histo.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
    for name, histo in delta_means_summary.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        #histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
    for name, histo in delta_sigmas_summary.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        #histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())


    for name, histo in unfoldeds.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
    for name, histo in unfolded_summary.iteritems():
        canvas = plotter.create_and_write_canvas_single(0, 21, 1, False, False, histo, write=False)
        histo.SetStats(True)
        canvas.Update()
        histo.Write()
        canvas.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
    for name, histo in unfolded_summary.iteritems():
        print name
        leg = LegendDefinition("Unfolding comparison", ['Truth', 'Unfolded'], 'NE')
        canvas = plotter.create_and_write_canvas_with_comparison('Pull_'+name,[1,0],[0,21],[2,1], leg, False, False, [true_distribution, histo], write=False, comparison = 'pull')
        canvas.Update()
        canvas.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
        canvas = plotter.create_and_write_canvas_with_comparison('Ratio_'+name,[1,0],[0,21],[2,1], leg, False, False, [true_distribution, histo], write=False, comparison = 'ratio')
        canvas.Update()
        canvas.Write()
        canvas.SaveAs('%s.png' % canvas.GetName())
        canvas.SaveAs('%s.pdf' % canvas.GetName())
        
        
    outfile.close()
    #os.chdir("..")
    
    
    
    os.chdir(curdir)
    #set_trace()

    #if not args.noshapes:
        ##Overlays the prefit values of the different shapes with the envelope of 
        ##what is fitted by the toys
        #out = os.path.join(args.out, 'shapes')
        #os.makedirs(out)
        #with io.root_open(args.harvested) as harvest:
            #has_prefit = hasattr(harvest, 'prefit')
            #prefit = harvest.prefit if has_prefit else None
            #toys = EnvelopeView(
                #*[harvest.get(i.GetName()).get(args.variable) 
                  #for i in harvest.keys() 
                  #if i.GetName().startswith('toy_')]
                #)
            #shapes = [i.GetName() for i in prefit.keys()]
      
        #for shape in shapes:
            #canvas = plotting.Canvas()
            #legend = plotting.Legend(4, rightmargin=0.07, topmargin=0.05, leftmargin=0.45)
            #legend.SetBorderSize(0)
            #if has_prefit:
                #pre_shape = prefit.Get(shape)
                #pre_shape.title = 'input shape'
                #pre_shape.legendstyle = 'p'
            #toy_shape = toys.Get(shape)
         
            #toy_shape.Draw()
            #if has_prefit:
                #pre_shape.Draw('P same')

            #legend.AddEntry(toy_shape.two_sigma)
            #legend.AddEntry(toy_shape.one_sigma)
            #legend.AddEntry(toy_shape.median)
            #if has_prefit:
                #legend.AddEntry(pre_shape)
            #legend.Draw()
        
            #canvas.Update()
            #canvas.SaveAs('%s/%s.png' % (out, shape))
            #canvas.SaveAs('%s/%s.pdf' % (out, shape))

        #def make_hist(key, vals, out, prefix=''):
            #canvas = plotting.Canvas()
            #low = min(vals)
        #low = low*0.8 if low > 0 else low*1.2
        #hi = max(vals)*1.2
        #hist = plotting.Hist(50, low, hi, title='')
        #for v in vals:
            #hist.Fill(v)
        #hist.Draw()
        #log.debug('%s has %i entries, %.0f integral' % (key, hist.GetEntries(), hist.Integral()))
        #canvas.SaveAs('%s/%s%s.png' % (out, prefix, key.replace('/','_')))
        ##canvas.SaveAs('%s/%s.pdf' % (out, key.replace('/','_')))


    #if not args.nopars:
        ##Plots the post-fit distribution of the POI and nuisances
        #out = os.path.join(args.out, 'floating_parameters')
        #os.makedirs(out)
        #with io.root_open(args.mlfit) as mlfit:
            #pars = {}
            #yields = {}
            #first = True
            #toys = [i.GetName() for i in mlfit.keys() if i.GetName().startswith('toy_')]
            #log.info('examining %i toys' % len(toys))
            #for toy in toys:
                #norms = ArgSet(
                    #mlfit.Get(
                        #join(
                            #toy,
                            #'norm_fit_s'
                            #)
                    #)
                #)
                #norms = [i for i in norms]

                #fit_result = mlfit.Get(
                    #join(
                        #toy,
                        #'fit_s'
                    #)
                #)
                #fit_pars = ArgList(fit_result.floatParsFinal())
         
                #if first:
                    #first = False
                    #pars.update( 
                        #dict(
                            #(i.GetName(), []) for i in fit_pars
                        #)
                    #)
                    #yields.update(
                        #dict(
                            #(i.GetName(), []) for i in norms
                        #)
                    #)

                #for i in norms:
                    #yields[i.GetName()].append(i.getVal())
                
                #for i in fit_pars:
                    #pars[i.GetName()].append(i.getVal())

            #for i, j in yields.iteritems():
                #make_hist(i, j, out, prefix='yield_')
    
            #for i, j in pars.iteritems():
                #make_hist(i, j, out, prefix='par_')

        #if not args.postpulls:
            #pass

