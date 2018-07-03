from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os, glob, sys
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
import rootpy.plotting as plotting
from rootpy.plotting import views, Graph, Hist, Hist2D
from argparse import ArgumentParser
import URAnalysis.Utilities.roottools as urtools
import ROOT
from URAnalysis.Utilities.tables import latex_table
from URAnalysis.PlotTools.views.RebinView import RebinView
import argparse
import matplotlib.pyplot as plt
import numpy as np
from styles import styles
import URAnalysis.Utilities.prettyjson as prettyjson


'''
This script opens for files from the most recent 2016 and 2017 jobids to make
plots comparing various btag and muon sf histograms.
'''

project = os.environ['URA_PROJECT']

plotter = BasePlotter(
    '%s/plots/2016_2017_Comparison' % project,
    #defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
    defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}},
)

comp_file_name = 'sf_comparison.root'

### check that files with MC hists exist
jobid_2016 = '2016Attempt'
jobid_2017 = 'deepctag_AllRuns'

if not (comp_file_name in os.listdir( '%s/plots/%s/ctageff/All_Runs/mass_discriminant/notag/notag' % (project, jobid_2016) )\
     or comp_file_name in os.listdir( '%s/plots/%s/ctageff/All_Runs/mass_discriminant/notag/notag' % (project, jobid_2017) ) ):
    print "The sf_comparison.root file doesn't exist in one or both of the directories."
    sys.exit()


### open root files
MCfile_2016 = root_open('%s/plots/%s/ctageff/All_Runs/mass_discriminant/notag/notag/sf_comparison.root' % (project, jobid_2016), 'read') 
MCfile_2017 = root_open('%s/plots/%s/ctageff/All_Runs/mass_discriminant/notag/notag/sf_comparison.root' % (project, jobid_2017), 'read') 

### define which variables to get
variables = [
    ('btag_sf', 'b-tag SF', 'Events', ''),
    ('btag_sf_up', 'b-tag SF, up sys', 'Events', ''),
    ('btag_sf_dw', 'b-tag SF, down sys', 'Events', ''),
    ('btag_sf_B', 'b-tag SF (b quarks)', 'Events', ''),
    ('btag_sf_B_up', 'b-tag SF (b quarks), up sys', 'Events', ''),
    ('btag_sf_B_dw', 'b-tag SF (b quarks), down sys', 'Events', ''),
    ('btag_sf_C', 'b-tag SF (c quarks)', 'Events', ''),
    ('btag_sf_C_up', 'b-tag SF (c quarks), up sys', 'Events', ''),
    ('btag_sf_C_dw', 'b-tag SF (c quarks), down sys', 'Events', ''),
    ('btag_sf_L', 'b-tag SF (l quarks)', 'Events', ''),
    ('btag_sf_L_up', 'b-tag SF (l quarks), up sys', 'Events', ''),
    ('btag_sf_L_dw', 'b-tag SF (l quarks), down sys', 'Events', ''),
    ('muon_sf_vs_mu_pt_etaL0p9', '$p_{T}(\mu)$ (GeV)', '$\mu$ SF', '|$\eta$|($\mu$) < 0.9'),
    ('muon_sf_vs_mu_pt_0p9eta1p2', '$p_{T}(\mu)$ (GeV)', '$\mu$ SF', '0.9 $\leq$ |$\eta$|($\mu$) < 1.2'),
    ('muon_sf_vs_mu_pt_1p2eta2p1', '$p_{T}(\mu)$ (GeV)', '$\mu$ SF', '1.2 $\leq$ |$\eta$|($\mu$) < 2.1')
]

#set_trace()
for var, xlabel, ylabel, title in variables:
    hist_2016 = asrootpy(MCfile_2016.Get(var)).Clone()
    hist_2017 = asrootpy(MCfile_2017.Get(var)).Clone()
    #set_trace()

    if hist_2016.DIM == 2:

            ### rebin hists, binning specific for each one
        xbins_2016 = np.array([30., 40., 50., 60., 120., 200., 500.])
        ybins_2016 = np.linspace(hist_2016.GetYaxis().GetBinLowEdge(1), hist_2016.GetYaxis().GetBinUpEdge(hist_2016.GetYaxis().GetNbins()), hist_2016.GetYaxis().GetNbins()+1 ) # ybinning remains unchanged
        hist_2016 = RebinView.newRebin2D(hist_2016, xbins_2016, ybins_2016)

        xbins_2017 = np.array([30., 32., 40., 50., 60., 120., 200., 500.])
        ybins_2017 = np.linspace(hist_2017.GetYaxis().GetBinLowEdge(1), hist_2017.GetYaxis().GetBinUpEdge(hist_2017.GetYaxis().GetNbins()), hist_2017.GetYaxis().GetNbins()+1 ) # ybinning remains unchanged
        hist_2017 = RebinView.newRebin2D(hist_2017, xbins_2017, ybins_2017)


        ## find x and y bin values for nonzero bins for 2016 2017 hists
        hist_2016_xbins = np.around( [hist_2016.GetXaxis().GetBinCenter(binx+1) for binx in range(hist_2016.GetNbinsX()) if hist_2016.Integral(binx+1, binx+1, 1, hist_2016.GetNbinsY()) != 0 ], decimals=3 )
        hist_2016_xbinlowedges = [hist_2016.GetXaxis().GetBinLowEdge(binx+1) for binx in range(hist_2016.GetNbinsX()) if hist_2016.Integral(binx+1, binx+1, 1, hist_2016.GetNbinsY()) != 0 ]
        hist_2016_xerrors = [ hist_2016_xbins[i]-hist_2016_xbinlowedges[i] for i in range(len(hist_2016_xbins)) ]

        hist_2016_nonzero_xbins = [binx+1 for binx in range(hist_2016.GetNbinsX()) if hist_2016.Integral(binx+1, binx+1, 1, hist_2016.GetNbinsY()) != 0 ]
        hist_2016_ybins = np.around( [hist_2016.GetYaxis().GetBinCenter(biny+1) for binx in hist_2016_nonzero_xbins for biny in range(hist_2016.GetNbinsY()) if hist_2016.Integral(binx, binx, biny+1, biny+1) != 0], decimals=3 )
            
        hist_2017_xbins = np.around( [hist_2017.GetXaxis().GetBinCenter(binx+1) for binx in range(hist_2017.GetNbinsX()) if hist_2017.Integral(binx+1, binx+1, 1, hist_2017.GetNbinsY()) != 0 ], decimals=3 )
        hist_2017_xbinlowedges = [hist_2017.GetXaxis().GetBinLowEdge(binx+1) for binx in range(hist_2017.GetNbinsX()) if hist_2017.Integral(binx+1, binx+1, 1, hist_2017.GetNbinsY()) != 0 ]
        hist_2017_xerrors = [ hist_2017_xbins[i]-hist_2017_xbinlowedges[i] for i in range(len(hist_2017_xbins)) ]

        hist_2017_nonzero_xbins = [binx+1 for binx in range(hist_2017.GetNbinsX()) if hist_2017.Integral(binx+1, binx+1, 1, hist_2017.GetNbinsY()) != 0 ]
        hist_2017_ybins = np.around( [hist_2017.GetYaxis().GetBinCenter(biny+1) for binx in hist_2017_nonzero_xbins for biny in range(hist_2017.GetNbinsY()) if hist_2017.Integral(binx, binx, biny+1, biny+1) != 0], decimals=3 )
            

        #set_trace()


        fig = plt.figure()
        plt.errorbar(hist_2016_xbins, hist_2016_ybins, xerr=hist_2016_xerrors, color='r', label='2016 MC', fmt='o', elinewidth=3, capthick=2)
        plt.errorbar(hist_2017_xbins, hist_2017_ybins, xerr=hist_2017_xerrors, color='b', label='2017 MC', fmt='o', elinewidth=3, capthick=2)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.tight_layout()
        plt.legend(numpoints=1)
        plt.xlim(0,500)
        plt.ylim(round(min(np.concatenate((hist_2016_ybins,hist_2017_ybins), axis=0))-0.01, 2), round(max(np.concatenate((hist_2016_ybins,hist_2017_ybins), axis=0))+0.01, 2) )
        plt.grid()
        fig.savefig('%s/%s_Comparison.png' % (plotter.outputdir, var))

        #set_trace()


    else:
        mean_2016 = hist_2016.GetMean()
        rms_2016 = hist_2016.GetRMS()
        mean_2017 = hist_2017.GetMean()
        rms_2017 = hist_2017.GetRMS()

        plotter.set_histo_style(hist_2016, fillstyle='hollow', color = 'r', title='2016 MC Mean=%.2f, RMS=%.2f' % (mean_2016, rms_2016) )
        plotter.set_histo_style(hist_2017, fillstyle='hollow', color = 'b', title='2017 MC Mean=%.2f, RMS=%.2f' % (mean_2017, rms_2017) )
        plotter.overlay([hist_2016, hist_2017], legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=ylabel, drawstyle='hist')
        plotter.save('%s_Comparison' % var )
    
