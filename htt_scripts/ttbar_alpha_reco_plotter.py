'''
ttbar_reco_3J Analyzer Plotter macro
'''

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
import ROOT

out_analyzer = 'ttbar_alpha_reco'

parser = argparse.ArgumentParser(description='Create plots using files from ttbar_reco_3J.')

jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']

parser.add_argument('--analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('--sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, AtoTT_M...).')
#parser.add_argument('--perm', help='Choose between matched_perms and 3-jet event perms (Event).')
parser.add_argument('--plot', help='Choose type of plots to generate (Gen, Reconstruction, Resolution, Discriminant, Alpha_Correction, All).')
args = parser.parse_args()


##### check analysis type
if not (args.analysis == "Test" or args.analysis == "Full"):
    print "You chose %s as your analysis type.\nYou must choose Full or Test!" % args.analysis
    sys.exit()

analyzer = 'ttbar_alpha_reco'


##### check sample type
results_files = []
results_file = '' 
if args.analysis == "Test":
    results_files = filter(lambda x: '%s.%s.test.root' % (args.sample, analyzer) in x, os.listdir(project))
    results_file = '%s/%s' % (project, results_files[0])    

    plotter = BasePlotter(
    	'%s/htt_scripts/plots/%s/%s/%s/%s' % (project, out_analyzer, jobid, args.analysis, args.sample),
    	#defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
    	defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}},
        styles = {
            'sample' : styles['ttJetsM0'] if args.sample=='ttJets' else styles[args.sample]
        }
    )
    #set_trace()

if args.analysis == "Full":
    results_files = filter(lambda x: '%s.root' % args.sample in x, os.listdir('%s/results/%s/%s' % (project, jobid, analyzer)))
    results_file  = '%s/results/%s/%s/%s' % (project, jobid, analyzer, results_files[0])

    plotter = BasePlotter(
    	'%s/plots/%s/%s/%s' % (project, jobid, out_analyzer, args.sample),
    	#defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
    	defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}},
        styles = {
            'sample' : styles['ttJetsM0'] if args.sample=='ttJets' else styles[args.sample]
        }
    )


allowed_plot_types = ['Gen', 'Reconstruction', 'Resolution', 'Discriminant', 'Alpha_Correction']
if not args.plot in allowed_plot_types+['All']:
    print "You chose %s as your plot type.\nYou must choose from the help list!" % args.plot
    sys.exit()



myfile = root_open(results_file, 'read')
normfile = views.NormalizeView( root_open(results_file, 'read') )
lumifile = open('%s/inputs/%s/%s.lumi' % (project, jobid, args.sample), 'read')


#set_trace()

##### Global Var. Definitions ####
defcol = 'black' # default color
defyax = 'A.U.' # yaxis title
if args.sample == "ttJetsM0":
    mass_min = 0
    mass_max = 2000
    nbins = 41
    mtt_bins = np.linspace(mass_min, mass_max, nbins)
if args.sample == "ttJetsM700":
    mass_min = 700
    mass_max = 1000
    nbins = 31
    mtt_bins = np.linspace(mass_min, mass_max, nbins)
if args.sample == "ttJetsM1000":
    mass_min = 1000
    mass_max = 2000
    nbins = 11
    mtt_bins = np.linspace(mass_min, mass_max, nbins)


m_range = plotter.styles['sample']['name']
decay = plotter.styles['sample']['decay']

##########
'''
This section creates dictionaries used for making plots based on different classifications
1. Categories- categorization of (MERGED or LOST) best perms based on their objects compared to the matched perm's
2. Objects- objects used in reconstruction of system variables and their labels
3. CUT- the label and color of each region from a cut made using the combined discriminant
4. kinvar- kinematic variables used to make plots and their axis labels
5. DiscPlots- the 3 different types of discriminants and their axis labels
6. var_types- the objects used for each kinematic variable and their axis labels
'''

#Perm_Categories = ['NO_MP', 'CORRECT_WJET_CORRECT_Bs', 'CORRECT_WJET_SWAPPED_Bs', 'CORRECT_WJET_CORRECT_BHAD', 'CORRECT_WJET_CORRECT_BLEP', 'CORRECT_WJET_WRONG_Bs',
#                    'WRONG_WJET_CORRECT_Bs', 'WRONG_WJET_SWAPPED_Bs', 'WRONG_WJET_CORRECT_BHAD', 'WRONG_WJET_CORRECT_BLEP', 'WRONG_WJET_WRONG_Bs'
#                ]
Correct_WJet_Categories = ['CORRECT_WJET_CORRECT_Bs', 'CORRECT_WJET_SWAPPED_Bs', 'CORRECT_WJET_CORRECT_BHAD', 'CORRECT_WJET_CORRECT_BLEP', 'CORRECT_WJET_WRONG_Bs']
Wrong_WJet_Categories = ['WRONG_WJET_CORRECT_Bs', 'WRONG_WJET_SWAPPED_Bs', 'WRONG_WJET_CORRECT_BHAD', 'WRONG_WJET_CORRECT_BLEP', 'WRONG_WJET_WRONG_Bs']


#Categories = {'MERGED' : ['RIGHT', 'MERGED_SWAP', 'MERGED', 'WRONG'],\
#              'LOST' : ['RIGHT', 'LOST_SWAP', 'LOST', 'WRONG']}
#Categories = {'MERGED' : [('RIGHT', 'red', 3345), ('MERGED_SWAP', 'orange', 0), ('MERGED', 'green', 3006), ('WRONG', 'blue', 3354)],\
#              'LOST' : [('RIGHT', 'red', 3345), ('LOST_SWAP', 'orange', 0), ('LOST', 'green', 3006), ('WRONG', 'blue', 3354)]}
Objects = {'THad' : 't_{h}', 'TLep' : 't_{l}', 'TTbar' : 't#bar{t}'}
CUT = {'LCut' : ['#lambda_{comb} < 2', 'black'], 'GCut' : ['#lambda_{comb} #geq 2', 'red']}
kinvar = {'Mass' : ' Mass [GeV]', 'Pt' : ' p_{T}', 'Eta' : ' #eta', 'Costh' : ' cos(#theta^{*})'}
DiscPlots = {'Massdisc' : ['#lambda_{mass}'], 'NSdisc' : ['#lambda_{NS}'], 'Totaldisc' : ['#lambda_{comb}']}

#var_types = {}
#if args.perm == 'Matched':
#    var_types = {'Costh' : {'THad' : 't_{h} cos(#theta^{*})', 'TLep' : 't_{l} cos(#theta^{*})'},\
#            'Mass' : {'THad' : 'M(t_{h}) [GeV]', 'TTbar' : 'M(t#bar{t}) [GeV]'},\
#            'Pt' : {'THad' : 'p_{T}(t_{h}) [GeV]', 'THad_PT_P' : 't_{h} sin(#theta)', 'THad_PZ_P' : 't_{h} cos(#theta)', 'TLep' : 'p_{T}(t_{l}) [GeV]', 'TTbar' : 'p_{T}(t#bar{t}) [GeV]'},\
#            'Eta' : {'THad' : '#eta(t_{h})', 'TLep' : '#eta(t_{l})'}}#, 'TTbar' : 't#bar{t}'} }
#else: #elif args.perm == 'Event':
var_types = {'Costh' : {'THad' : 't_{h} cos(#theta^{*})', 'TLep' : 't_{l} cos(#theta^{*})'},\
        'Mass' : {'THad' : 'M(t_{h}) [GeV]', 'TTbar' : 'M(t#bar{t}) [GeV]'},\
        'Pt' : {'THad' : 'p_{T}(t_{h}) [GeV]', 'TLep' : 'p_{T}(t_{l}) [GeV]', 'TTbar' : 'p_{T}(t#bar{t}) [GeV]'},\
        'Eta' : {'THad' : '#eta(t_{h})', 'TLep' : '#eta(t_{l})'}}#, 'TTbar' : 't#bar{t}'} }
##########

## beginning of plot function definitions


##############################################################################################################################
def alpha_corrections(directory, subdir):

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]

    print "These plots aren't currently supported.\nThe alpha_hists.root file to be used for the alpha correction is made from new_alpha_interp.py in my PYTHON_FILES/ttbar Toshiba directory.\n"
    return

    #fitvars = [ # variables for hists
    #    ('THad_E/Alpha_THad_E', '173.1/Reco M(t_{h})', '#alpha_{E} = Gen E(t_{h})/Reco E(t_{h})', 'M($t\\overline{t}$) < $\\infty$'),
    #    ('THad_P/Alpha_THad_P', '173.1/Reco M(t_{h})', '#alpha_{P} = Gen P(t_{h})/Reco P(t_{h})', 'M($t\\overline{t}$) < $\\infty$'),
    #    ('THad_M/Alpha_THad_M', '173.1/Reco M(t_{h})', '#alpha_{M} = Gen M(t_{h})/Reco M(t_{h})', 'M($t\\overline{t}$) < $\\infty$'),
    #]

    fitvars = [
        ('THad_E/Alpha_THad_E_Mtt_vs_Mthad_vs_Alpha'),# '173.1/Reco M(t_{h})', '#alpha_{E} = Gen E(t_{h})/Reco E(t_{h})', 'M($t\\overline{t}$) < $\\infty$'),
        ('THad_P/Alpha_THad_P_Mtt_vs_Mthad_vs_Alpha'),# '173.1/Reco M(t_{h})', '#alpha_{P} = Gen P(t_{h})/Reco P(t_{h})', 'M($t\\overline{t}$) < $\\infty$'),
        ('THad_M/Alpha_THad_M_Mtt_vs_Mthad_vs_Alpha'),# '173.1/Reco M(t_{h})', '#alpha_{M} = Gen M(t_{h})/Reco M(t_{h})', 'M($t\\overline{t}$) < $\\infty$'),
    ]

    #mtt_ranges = np.array([200, 350, 400, 500, 700, 1000, np.Inf])
    #for mtt in range(len(mtt_ranges)-1):
    #    if mtt_ranges[mtt+1] == np.Inf:
    #        txt_box = "M($t\\overline{t}$) > %.0f" % mtt_ranges[mtt]
    #        Evar = ('THad_E/Alpha_THad_E_Mttbar%.0ftoInf' % mtt_ranges[mtt], '173.1/Reco M(t_{h})', '#alpha_{E} = Gen E(t_{h})/Reco E(t_{h})', txt_box )
    #        Pvar = ('THad_P/Alpha_THad_P_Mttbar%.0ftoInf' % mtt_ranges[mtt], '173.1/Reco M(t_{h})', '#alpha_{P} = Gen P(t_{h})/Reco P(t_{h})', txt_box )
    #        Mvar = ('THad_M/Alpha_THad_M_Mttbar%.0ftoInf' % mtt_ranges[mtt], '173.1/Reco M(t_{h})', '#alpha_{M} = Gen M(t_{h})/Reco M(t_{h})', txt_box )
    #    else:
    #        txt_box = "%.0f $\leq$ M($t\\overline{t}$) < %.0f" % (mtt_ranges[mtt], mtt_ranges[mtt+1])
    #        Evar = ('THad_E/Alpha_THad_E_Mttbar%.0fto%.0f' % (mtt_ranges[mtt], mtt_ranges[mtt+1]), '173.1/Reco M(t_{h})', '#alpha_{E} = Gen E(t_{h})/Reco E(t_{h})', txt_box )
    #        Pvar = ('THad_P/Alpha_THad_P_Mttbar%.0fto%.0f' % (mtt_ranges[mtt], mtt_ranges[mtt+1]), '173.1/Reco M(t_{h})', '#alpha_{P} = Gen P(t_{h})/Reco P(t_{h})', txt_box )
    #        Mvar = ('THad_M/Alpha_THad_M_Mttbar%.0fto%.0f' % (mtt_ranges[mtt], mtt_ranges[mtt+1]), '173.1/Reco M(t_{h})', '#alpha_{M} = Gen M(t_{h})/Reco M(t_{h})', txt_box )
    #    fitvars.append(Evar)
    #    fitvars.append(Pvar)
    #    fitvars.append(Mvar)


    ### create dictionary for fit values for E and P for each mttbar range
    #alphas = list(set([ fitvars[i][0].split('/')[0] for i in range(len(fitvars))]))
    #mranges = list(set([ fitvars[i][0].split('/')[1].split('_')[-1] for i in range(len(fitvars)) if 'Mttbar' in fitvars[i][0].split('/')[1] ]))+['All']


        ## create file with directories for 'THad_E', 'THad_P', and 'THad_M'
    with root_open('%s/alpha_hists.root' % '/'.join([project, 'inputs', jobid, 'INPUT']), 'w') as out:# write to $jobid/INPUT directory
        for hname in fitvars:
            out.mkdir(hname.split('/')[0]).cd()

    #set_trace()

    for hvar in fitvars:
        hname = '/'.join([directory, 'Alpha_Correction', 'CORRECT_WJET_CORRECT_Bs', hvar])
        hist = asrootpy(myfile.Get(hname)).Clone()

        if hist.Integral() == 0:
            continue

        #plotter.plot(hist)

        colors = ['green', 'red', 'black', 'blue', 'orange', 'magenta', 'cyan', 'yellow']
        #colors = ['green', 'red', 'black', 'blue']
        i = 0

        #xbins = np.linspace(0.9, 2.5, 9) ## rebin xaxis so [0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
        ##xbins = np.linspace(0.9, 1.5, 7) 
        #ybins = np.linspace( hist.GetYaxis().GetBinLowEdge(1), hist.GetYaxis().GetBinUpEdge(hist.GetYaxis().GetNbins()), hist.GetYaxis().GetNbins()+1 ) # ybinning remains unchanged

            ## define bin edges for future rebinning
        mthad_bins = np.linspace(0.9, 2.5, 9)
        mtt_bins = np.array([200., 350., 400., 500., 700., 1000., 2000.])
        #mtt_bins = np.linspace( hist.GetYaxis().GetBinLowEdge(1), hist.GetYaxis().GetBinUpEdge(hist.GetNbinsY()), hist.GetNbinsY()+1 )
        alpha_bins = np.linspace( hist.GetZaxis().GetBinLowEdge(1), hist.GetZaxis().GetBinUpEdge(hist.GetNbinsZ()), hist.GetNbinsZ()+1 )

        #set_trace()
        hist = RebinView.newRebin3D(hist, mthad_bins, mtt_bins, alpha_bins)

            ## initialize dictionary to keep median values for different dists of mtt bins
        alphas_dict = {'All': {}, 'Mtt' : {} }


        #set_trace()
        #hist = RebinView.newRebin2D(hist, xbins, ybins)
        #hist.xaxis.range_user = min(xbins), max(xbins)
        ##hist.xaxis.range_user = 0.9, 2.5
        ##hist.yaxis.range_user = 0.0, 4.0

        #yprojections = []
        #yproj_norms = []

            ## get medians for single bin of mtt
        medians, median_errors, median_weights = plotter.median_from_3d_hist(hist, projection='zx', xbins=mthad_bins, ybins=alpha_bins)
        alphas_dict['All'] = {'Medians' : medians, 'Errors' : median_errors, 'Weights' : median_weights }


            ## get medians for all bins of mtt
        for mtt_bin in range(1, len(mtt_bins)):

            mtt_yslice = hist.Clone()
            mtt_yslice.GetYaxis().SetRange(mtt_bin, mtt_bin+1)

            #set_trace() 
            medians, median_errors, median_weights = plotter.median_from_3d_hist(mtt_yslice, projection='zx', xbins=mthad_bins, ybins=alpha_bins)
            alphas_dict['Mtt']['Bin%s' % mtt_bin] = {'Medians' : medians, 'Errors' : median_errors, 'Weights' : median_weights }

        #set_trace()

            ## fill hists with medians in root file
        with root_open('%s/alpha_hists.root' % '/'.join([project, 'inputs', jobid, 'INPUT']), 'update') as out:# write to $jobid/INPUT directory

                ## fill alphas for single mtt bin
            mttbar_bins = np.array([min(mtt_bins), max(mtt_bins)])
            all_mtt = Hist2D(mthad_bins, mttbar_bins, name='%s_All' % hvar.split('/')[0], title='Median #alpha')
            all_mtt.set_x_title(hist.xaxis.title)
            all_mtt.set_y_title(hist.yaxis.title)

            #set_trace()
            for binx in range( 1, len(alphas_dict['All']['Medians'])+1 ):
                if alphas_dict['All']['Medians'][binx-1] < 0: continue # skip invalid alphas that have negative values

                all_mtt[binx, 1].value = alphas_dict['All']['Medians'][binx-1]
                all_mtt[binx, 1].error = alphas_dict['All']['Errors'][binx-1]

            all_mtt.Write()

                ## fill alphas for different mtt bins
            mttbar_bins = mtt_bins
            binned_mtt = Hist2D(mthad_bins, mttbar_bins, name='%s_Mtt' % hvar.split('/')[0], title='Median #alpha')
            binned_mtt.set_x_title(hist.xaxis.title)
            binned_mtt.set_y_title(hist.yaxis.title)

            for biny in range( 1, len(alphas_dict['Mtt'].keys())+1 ):
                for binx in range( 1, len(alphas_dict['Mtt']['Bin%s' % biny]['Medians'])+1 ):
                    #set_trace()
                    if alphas_dict['Mtt']['Bin%s' % biny]['Medians'][binx-1] < 0: continue # skip invalid alphas that have negative values

                    binned_mtt[binx, biny].value = alphas_dict['Mtt']['Bin%s' % biny]['Medians'][binx-1]
                    binned_mtt[binx, biny].error = alphas_dict['Mtt']['Bin%s' % biny]['Errors'][binx-1]

            binned_mtt.Write()
        #set_trace()

        #    #plotter.set_histo_style(hist_yproj, color='black', title='%s-%s Med=%.2f, Max=%.2f' % (hist.GetXaxis().GetBinLowEdge(xbin), hist.GetXaxis().GetBinUpEdge(xbin), median, maxval))
        #    #plotter.plot(hist_yproj, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=ylabel, ytitle=defyax, drawstyle='hist', x_range=(0.0, 5.0))
        #    #plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', hvar.split('/')[0], 'y_proj']))
        #    #plotter.save('%s_%s' % (hvar.split('/')[1], hist.GetXaxis().GetBinLowEdge(xbin)) )
        #    #set_trace()

        #    plotter.set_histo_style(hist_yproj, color=colors[i], title='%s-%s Med=%.2f, Max=%.2f' % (hist.GetXaxis().GetBinLowEdge(xbin), hist.GetXaxis().GetBinUpEdge(xbin), median, maxval))
        #    plotter.plot(hist_yproj, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=ylabel, ytitle=defyax, drawstyle='hist', x_range=(0.0, 5.0))
        #    yprojections.append(hist_yproj)
        #    yproj_norms.append(hist_yproj/hist_yproj.Integral())
        #    i += 1
        #    #set_trace()

        #for xbin in range(hist.GetNbinsX() + 1):
        #    if hist.GetXaxis().GetBinLowEdge(xbin) >= hist.xaxis.range_user[0] and hist.GetXaxis().GetBinUpEdge(xbin) <= hist.xaxis.range_user[1]:
        #        hist_yproj = hist.ProjectionY("", xbin, xbin).Clone()
        #        hist_yproj = asrootpy(hist_yproj)
        #        if hist_yproj.Integral() == 0:
        #            continue

        #        probs = np.zeros(hist.GetYaxis().GetNbins()+1)
        #        binvals = np.zeros(hist.GetYaxis().GetNbins()+1)

        #        for bins in range(hist.GetYaxis().GetNbins()+1):
        #            probs[bins] = hist_yproj.Integral(1, bins+1)/hist_yproj.Integral()
        #            binvals[bins] = hist_yproj.GetBinContent(bins)

        #        median = hist_yproj.GetBinCenter(np.where(probs>= 0.5)[0][0])
        #        median_error = 1.2533*hist_yproj.GetMeanError() #standard error of median
        #        #set_trace()
        #        if median_error == 0:
        #            #set_trace()
        #            median_weight = 1000 # 1/(standard error of median) to be used in fit
        #        else:
        #            median_weight = 1/(1.2533*hist_yproj.GetMeanError()) # 1/(standard error of median) to be used in fit
        #        maxval = round(hist_yproj.GetBinCenter(np.where(binvals==binvals.max())[0][0]), 2)

        #        medians.append(median)
        #        median_weights.append(median_weight)
        #        median_errors.append(median_error)
        #        maxvals.append(maxval)

        #        #plotter.set_histo_style(hist_yproj, color='black', title='%s-%s Med=%.2f, Max=%.2f' % (hist.GetXaxis().GetBinLowEdge(xbin), hist.GetXaxis().GetBinUpEdge(xbin), median, maxval))
        #        #plotter.plot(hist_yproj, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=ylabel, ytitle=defyax, drawstyle='hist', x_range=(0.0, 5.0))
        #        #plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', hvar.split('/')[0], 'y_proj']))
        #        #plotter.save('%s_%s' % (hvar.split('/')[1], hist.GetXaxis().GetBinLowEdge(xbin)) )
        #        #set_trace()

        #        plotter.set_histo_style(hist_yproj, color=colors[i], title='%s-%s Med=%.2f, Max=%.2f' % (hist.GetXaxis().GetBinLowEdge(xbin), hist.GetXaxis().GetBinUpEdge(xbin), median, maxval))
        #        plotter.plot(hist_yproj, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=ylabel, ytitle=defyax, drawstyle='hist', x_range=(0.0, 5.0))
        #        yprojections.append(hist_yproj)
        #        yproj_norms.append(hist_yproj/hist_yproj.Integral())
        #        i += 1
        #        #set_trace()

        #plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', hvar.split('/')[0], 'y_proj']))
        #plotter.overlay(yprojections, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=ylabel, ytitle=defyax, drawstyle='hist', x_range=(0.0, 5.0))
        #plotter.save('%s_yprojections' % hvar.split('/')[1])
        #plotter.overlay(yproj_norms, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=ylabel, ytitle=defyax, drawstyle='hist', x_range=(0.0, 5.0))
        #plotter.save('%s_yprojections_Norm' % hvar.split('/')[1])
    
        #xbin_ints = [hist.Integral(xbinid+1, xbinid+1, 1, hist.GetNbinsY()+1) for xbinid in range(hist.GetNbinsX())]
        #xbin_ratios = [xbin_ints[i]/sum(xbin_ints) for i in range(len(xbin_ints))]

        #hist.yaxis.range_user = 0.0, 4.0

        ##### section to fit polynomials to median values from xaxis bins
        ##xfit_bins = np.linspace(hist.GetXaxis().GetXmin(),hist.GetXaxis().GetXmax(), 200) #makes points to be used for fit results
        #xfit_bins = np.linspace(0, 3.0, 200) #makes points to be used for fit results
        ##set_trace()
        #hist.GetXaxis().GetCenter(xbins) # change array of xbin values from edges to centers

        #fit_1d_groom = np.polyfit(xbins[:-1], medians, 1, full=True, w=median_weights) # fit medians to 1 degree polynomial for 0.9 < xbin < 2.5
        #fit_2d_groom = np.polyfit(xbins[:-1], medians, 2, full=True, w=median_weights) # fit medians to 2 degree polynomial for 0.9 < xbin < 2.5

        ##set_trace()
        #p1_groom = np.poly1d(fit_1d_groom[0]) #gets parameters for 1 degree fit getting rid of first 2 points
        #p2_groom = np.poly1d(fit_2d_groom[0]) #gets parameters for 2 degree fit getting rid of first 2 points

        #fig = plt.figure()
        #fit_medians = plt.errorbar(xbins[:-1], medians, yerr=median_errors, label='Medians', linestyle='None', color='black', marker='.') #plots median values with errors
        #    # 1 degree fits
        #p1_groom_fit = plt.plot(xfit_bins, p1_groom(xfit_bins), label='weighted deg=1 bins > 0.9', linestyle='--', color='red') #plots p1 for bins > 0.9
        #    # 2 degree fits
        #p2_groom_fit = plt.plot(xfit_bins, p2_groom(xfit_bins), label='weighted deg=2 bins > 0.9', linestyle='--', color='green') #plots p2 for bins > 0.9

        #### fill dictionary with fit values
        #mrange = 'All' if 'Mttbar' not in hvar else hvar.split('/')[1].split('_')[-1]
        #fits[hvar.split('/')[0]][mrange] = { '1d' : fit_1d_groom[0].tolist(), '2d' : fit_2d_groom[0].tolist() }
        ##fits[hvar.split('/')[0]][mrange]['1D_g0.9'] = {}
        ##fits[hvar.split('/')[0]][mrange]['2D_g0.9'] = {}

        ##fits[hvar.split('/')[0]][mrange]['1D_g0.9'] = { 'Pars' : fit_1d_groom[0].tolist() }
        ##fits[hvar.split('/')[0]][mrange]['2D_g0.9'] = { 'Pars' : fit_2d_groom[0].tolist() }
        ##fits[hvar.split('/')[0]][mrange]['1D_g0.9'] = { 'Pars' : fit_1d_groom[0].tolist(), 'Residual' : '%.4f' % float(fit_1d_groom[1]) }
        ##fits[hvar.split('/')[0]][mrange]['2D_g0.9'] = { 'Pars' : fit_2d_groom[0].tolist(), 'Residual' : '%.4f' % float(fit_2d_groom[1]) }

        ##set_trace()

        #plt.xlabel('173.1/Reco M($t_{h}$)')
        #if 'THad_E' in hvar:
        #    plt.ylabel('$\\alpha_{E}$ = Gen E($t_{h}$)/Reco E($t_{h}$)')
        #elif 'THad_P' in hvar:
        #    plt.ylabel('$\\alpha_{P}$ = Gen P($t_{h}$)/Reco P($t_{h}$)')
        #else: # 'THad_M'
        #    plt.ylabel('$\\alpha_{M}$ = Gen M($t_{h}$)/Reco M($t_{h}$)')
        #plt.title(txt_box_label)
        ##set_trace()
        #plt.ylim(0.0, 4.0)
        #plt.grid()
        #plt.legend(loc='upper right',fontsize=8, numpoints=1)
        #plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', hvar.split('/')[0]]))
        #plotting_dir = plotter.outputdir
        #fig.savefig('%s/%s_fits' % (plotting_dir, hvar.split('/')[1]))

        #plotter.set_histo_style(hist, xtitle=xlabel, ytitle=ylabel)
        #plotter.plot(hist, drawstyle='colz')

        #mean_x = hist.ProjectionX().GetMean()
        #rms_x = hist.ProjectionX().GetRMS()
        #mean_y = hist.ProjectionY().GetMean()
        #rms_y = hist.ProjectionY().GetRMS()

        #box2 = plotter.make_text_box('X Mean=%.2f, RMS=%.2f\nY Mean=%.2f, RMS=%.2f' % (mean_x, rms_x, mean_y, rms_y), position='NE')
        #box2.Draw()
        ##plotter.set_subdir('/'.join([subdir, 'Alpha_Correction']))

        #plotter.save(hvar.split('/')[1])

        ##set_trace()
        #median_dict[hvar.split('/')[0]][mrange] = { 'medians' : medians, 'mthad': xbins[:-1].tolist(), 'median_weights' : median_weights }#[ (round(medians[i], 2), round(xbins[i], 2)) for i in range(len(medians))]


    #xranges = np.linspace(0.9, 2.5, 9) ## rebin xaxis so [0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
    #hvars = []
    #for num in range(len(xranges)-1):
    #    txt_box = "%.1f < 173.1/Reco M(t_{h}) #leq %.1f" % (xranges[num], xranges[num+1])
    #    Evar = ('THad_E/Gen_vs_Reco_THadE_%.1fto%.1f' % (xranges[num], xranges[num+1]), 'Reco E(t_{h})', 'Gen E(t_{h})', txt_box )
    #    Pvar = ('THad_P/Gen_vs_Reco_THadP_%.1fto%.1f' % (xranges[num], xranges[num+1]), 'Reco E(t_{h})', 'Gen E(t_{h})', txt_box )
    #    hvars.append(Evar)
    #    hvars.append(Pvar)

    #for hvar, xlabel, ylabel, txt_box_label in hvars:
    #    hname = '/'.join([directory, 'Alpha_Correction', hvar])
    #    hist = asrootpy(myfile.Get(hname)).Clone()

    #    if hist.Integral() == 0:
    #        continue

    #    plotter.plot(hist)
    #    #if 'Gen_vs_Reco' in hvar:
    #    box1 = plotter.make_text_box( txt_box_label, position='NE')
    #    box1.Draw()
    #    if 'THadE' in hvar:
    #        plotter.set_subdir('/'.join([subdir, 'Alpha_Correction/THad_E/Gen_vs_Reco']))
    #    if 'THadP' in hvar:
    #        plotter.set_subdir('/'.join([subdir, 'Alpha_Correction/THad_P/Gen_vs_Reco']))
    #    #else:
    #    #    plotter.set_subdir('/'.join([subdir, 'Alpha_Correction']))

    #    plotter.save(hvar.split('/')[1])
        
        #set_trace()

        #### create hists based on perm category (RIGHT/WRONG, etc...)
        #if topology:
        #    for cat in Categories[topology]:
        #        plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', cat]))
        #        #evt_col=styles[cat]['fillcolor']
        #        evt_type = styles[cat]['name']
        #        cat_hname = '/'.join([directory, 'Alpha_Correction', cat, hvar])
        #        cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

        #        if cat_hist.Integral() == 0:
        #            continue

        #        plotter.plot(cat_hist)
    
        #        cat_hist.Draw('colz')
        #        cat_hist.xaxis.set_title(xlabel)
        #        cat_hist.yaxis.set_title(ylabel)

        #        #r = cat_hist.Fit("pol1", "Q") # fit hist with polynomial of degree 1 with no output

        #        box1 = plotter.make_text_box( cat, position='NW')
        #        box1.Draw()
    
        #        plotter.save('%s_%s' % (hvar, cat))


    #mttbar_ranges = [ 'Mttbar200to350', 'Mttbar350to400', 'Mttbar400to500', 'Mttbar500to700', 'Mttbar700to1000', 'Mttbar1000toInf']
    #mtt_bins = np.array([275, 375, 450, 600, 850, 1200])
    #mthad_vals = { '1.0' : '0.9 < 173.1/Reco M($t_{h}$) $\leq$ 1.1', '1.2' : '1.1 < 173.1/Reco M($t_{h}$) $\leq$ 1.3', '1.4' : '1.3 < 173.1/Reco M($t_{h}$) $\leq$ 1.5',
    #                '1.6' : '1.5 < 173.1/Reco M($t_{h}$) $\leq$ 1.7', '1.8' : '1.7 < 173.1/Reco M($t_{h}$) $\leq$ 1.9', '2.0' : '1.9 < 173.1/Reco M(t_{h}) $\leq$ 2.1',
    #                '2.2' : '2.1 < 173.1/Reco M($t_{h}$) $\leq$ 2.3', '2.4' : '2.3 < 173.1/Reco M($t_{h}$) $\leq$ 2.5'
    #            }

    #xmin, xmax, ymin, ymax = 200, 1300, 0.0, 3.0

    #for alpha_key in median_dict.keys():
    #    ## initialize fit results dict to include mthad
    #    for mthad_val in mthad_vals.keys():
    #        fits[alpha_key]['MTHad%s' % mthad_val] = {}


    #        ## make plot comparing median values for all mtt bins
    #    fig = plt.figure()
    #    i=0
    #    for mttbar_range in mttbar_ranges:

    #        meds = np.array(median_dict[alpha_key][mttbar_range]['medians'])
    #        mthad_bins = np.array(median_dict[alpha_key][mttbar_range]['mthad'])

    #        plt.plot(mthad_bins, meds, linestyle='-', color=colors[i], label=mttbar_range)
    #        i += 1

    #    plt.xlabel('173.1/Reco M($t_{h}$)')
    #    if alpha_key == 'THad_E':
    #        plt.ylabel('$\\alpha_{E}$ = Gen E($t_{h}$)/Reco E($t_{h}$)')
    #    elif alpha_key == 'THad_P':
    #        plt.ylabel('$\\alpha_{P}$ = Gen P($t_{h}$)/Reco P($t_{h}$)')
    #    else: # 'THad_M'
    #        plt.ylabel('$\\alpha_{M}$ = Gen M($t_{h}$)/Reco M($t_{h}$)')
    #    #plt.title(txt_box_label)
    #    #set_trace()
    #    plt.ylim(0.5, 2.5)
    #    plt.grid()
    #    plt.legend(loc='upper left',fontsize=8, numpoints=1)
    #    plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', alpha_key]))
    #    plotting_dir = plotter.outputdir
    #    #set_trace()
    #    #fig.savefig('test')
    #    fig.savefig('%s/Alpha_vs_mthad_compare_Mtt' % plotting_dir)
    #    print '%s/Alpha_vs_mthad_compare_Mtt.png has been created' % plotting_dir

    #    for idx, binval in enumerate(sorted(mthad_vals.keys())):
    #        #set_trace()
    #        alpha_medians = np.array([median_dict[alpha_key][key]['medians'][idx] for key in mttbar_ranges])
    #        alpha_median_weights = np.array([median_dict[alpha_key][key]['median_weights'][idx] for key in mttbar_ranges])
    #        alpha_median_errors = np.array([1/median_dict[alpha_key][key]['median_weights'][idx] for key in mttbar_ranges])

    #            ## get results of 1d and 2d fits
    #        xfit_bins = np.linspace(xmin, xmax, (xmax-xmin)*2+1) #makes points to be used for fit results

    #        fit_1d = np.polyfit(mtt_bins, alpha_medians, 1, full=True, w=alpha_median_weights) # fit medians to 1 degree polynomial
    #        fit_2d = np.polyfit(mtt_bins, alpha_medians, 2, full=True, w=alpha_median_weights) # fit medians to 2 degree polynomial

    #        poly1d = np.poly1d(fit_1d[0]) #gets parameters for 1 degree fit
    #        poly2d = np.poly1d(fit_2d[0]) #gets parameters for 2 degree fit

    #            ## plot medians, 1d and 2d fits
    #        fig, ax = plt.subplots()
    #        alpha_vs_mttbar = plt.errorbar(mtt_bins, alpha_medians, yerr=alpha_median_errors, label='Medians', linestyle='None', color='black', marker='.') #plots median values with errors
    #        p1_fit = plt.plot(xfit_bins, poly1d(xfit_bins), label='weighted deg=1', linestyle='-', color='red') #plots p1 for bins > 0.9
    #        p2_fit = plt.plot(xfit_bins, poly2d(xfit_bins), label='weighted deg=2', linestyle='-', color='blue') #plots p2 for bins > 0.9

    #            ## format hists
    #        plt.xlabel('M($t\\overline{t}$)')
    #        if alpha_key == 'THad_E':
    #            plt.ylabel('$\\alpha_{E}$ = Gen E($t_{h}$)/Reco E($t_{h}$)')
    #        elif alpha_key == 'THad_P':
    #            plt.ylabel('$\\alpha_{P}$ = Gen P($t_{h}$)/Reco P($t_{h}$)')
    #        else: # 'THad_M'
    #            plt.ylabel('$\\alpha_{M}$ = Gen M($t_{h}$)/Reco M($t_{h}$)')
    #        plt.title(mthad_vals[binval])
    #        plt.xlim(xmin, xmax)
    #        plt.ylim(ymin, ymax)
    #        ax.set_xticks([200, 350, 400, 500, 700, 1000])
    #        ax.xaxis.grid(True, which='major')
    #        ax.set_yticks(np.linspace(ymin,ymax,(ymax-ymin)*2+1).tolist())
    #        ax.yaxis.grid(True, which='major')
    #        #set_trace()
    #        plt.legend(loc='upper right',fontsize=8, numpoints=1)
    #        plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', alpha_key]))
    #        plotting_dir = plotter.outputdir
    #        #set_trace()
    #        #fig.savefig('test')
    #        fig.savefig('%s/alpha_vs_mttbar_mthad%s_fits' % (plotting_dir, binval.replace('.', 'p')))
    #        print '%s/alpha_vs_mttbar_mthad%s_fits.png has been created' % (plotting_dir, binval.replace('.', 'p'))

    #        fits[alpha_key]['MTHad%s' % binval] = { '2d' : fit_2d[0].tolist() }
    #        #set_trace()
    #

    #    ## write fit parameters dict to json files
    #with open('%s/median_Mtt_Mthad.json' % '/'.join([plotter.base_out_dir, subdir, 'Alpha_Correction']), 'w') as f: # write to same dir as plots
    #    f.write(prettyjson.dumps(median_dict))
    ##set_trace()

    ##fit_pars_root('%s/fit_parameters' % '/'.join([project, 'inputs', jobid, 'INPUT']), fits)

    #with open('%s/fit_parameters.json' % '/'.join([plotter.base_out_dir, subdir, 'Alpha_Correction']), 'w') as f: # write to same dir as plots
    #    f.write(prettyjson.dumps(fits))
    #with open('%s/fit_parameters.json' % '/'.join([project, 'inputs', jobid, 'INPUT']), 'w') as f: # write to $jobid/INPUT directory
    #    f.write(prettyjson.dumps(fits))



def fit_pars_root(fname, fit_dict):

    mtt_xbins = np.array([200, 350, 400, 500, 700, 1000, 2000]) #bins of mtt
    all_xbins = np.array([200, 2000]) #bins of mtt
    with root_open('%s.root' % fname, 'w') as out:
        for fit_var in fit_dict.keys():
            out.mkdir(fit_var).cd()
            for fit_degree in ['1D', '2D']:
                out.mkdir(fit_var+'/'+fit_degree).cd()
                out.cd(fit_var+'/'+fit_degree)
                if fit_degree == '1D':
                    Mtt_slope = Hist(mtt_xbins, name='Mtt_slope', title='Mtt slope')
                    for i in range(Mtt_slope.GetNbinsX()):
                        hname = [j for j in fit_dict[fit_var].keys() if 'Mttbar' in j if int(j.split('to')[0].split('Mttbar')[1]) == int(Mtt_slope.GetBinLowEdge(i+1))][0] ## get 'Mttbar...to...' hist corresponding to bin
                        Mtt_slope[i+1].value = fit_dict[fit_var][hname]['%s_g0.9' % fit_degree]['Pars'][0] # set bin value to slope

                    Mtt_slope.Write()

                    Mtt_yint = Hist(mtt_xbins, name='Mtt_yint', title='Mtt yint')
                    for i in range(Mtt_yint.GetNbinsX()):
                        hname = [j for j in fit_dict[fit_var].keys() if 'Mttbar' in j if int(j.split('to')[0].split('Mttbar')[1]) == int(Mtt_yint.GetBinLowEdge(i+1))][0] ## get 'Mttbar...to...' hist corresponding to bin
                        Mtt_yint[i+1].value = fit_dict[fit_var][hname]['%s_g0.9' % fit_degree]['Pars'][1] # set bin value to slope

                    Mtt_yint.Write()

                    All_slope = Hist(all_xbins, name='All_slope', title='All slope')
                    for i in range(All_slope.GetNbinsX()):
                        hname = [j for j in fit_dict[fit_var].keys() if 'All' in j][0] ## get 'All' hist corresponding to bin
                        All_slope[i+1].value = fit_dict[fit_var][hname]['%s_g0.9' % fit_degree]['Pars'][0] # set bin value to slope

                    All_slope.Write()

                    All_yint = Hist(mtt_xbins, name='All_yint', title='All yint')
                    for i in range(All_yint.GetNbinsX()):
                        hname = [j for j in fit_dict[fit_var].keys() if 'All' in j][0] ## get 'All' hist corresponding to bin
                        All_yint[i+1].value = fit_dict[fit_var][hname]['%s_g0.9' % fit_degree]['Pars'][1] # set bin value to slope

                    All_yint.Write()

                if fit_degree == '2D':
                    Mtt_c0 = Hist(mtt_xbins, name='Mtt_c0', title='Mtt c0')
                    for i in range(Mtt_c0.GetNbinsX()):
                        hname = [j for j in fit_dict[fit_var].keys() if 'Mttbar' in j if int(j.split('to')[0].split('Mttbar')[1]) == int(Mtt_c0.GetBinLowEdge(i+1))][0] ## get 'Mttbar...to...' hist corresponding to bin
                        Mtt_c0[i+1].value = fit_dict[fit_var][hname]['%s_g0.9' % fit_degree]['Pars'][0] # set bin value to c0

                    Mtt_c0.Write()

                    Mtt_c1 = Hist(mtt_xbins, name='Mtt_c1', title='Mtt c1')
                    for i in range(Mtt_c1.GetNbinsX()):
                        hname = [j for j in fit_dict[fit_var].keys() if 'Mttbar' in j if int(j.split('to')[0].split('Mttbar')[1]) == int(Mtt_c1.GetBinLowEdge(i+1))][0] ## get 'Mttbar...to...' hist corresponding to bin
                        Mtt_c1[i+1].value = fit_dict[fit_var][hname]['%s_g0.9' % fit_degree]['Pars'][1] # set bin value to slope

                    Mtt_c1.Write()

                    Mtt_c2 = Hist(mtt_xbins, name='Mtt_c2', title='Mtt c2')
                    for i in range(Mtt_c2.GetNbinsX()):
                        hname = [j for j in fit_dict[fit_var].keys() if 'Mttbar' in j if int(j.split('to')[0].split('Mttbar')[1]) == int(Mtt_c2.GetBinLowEdge(i+1))][0] ## get 'Mttbar...to...' hist corresponding to bin
                        Mtt_c2[i+1].value = fit_dict[fit_var][hname]['%s_g0.9' % fit_degree]['Pars'][2] # set bin value to slope

                    Mtt_c2.Write()

                    All_c0 = Hist(all_xbins, name='All_c0', title='All c0')
                    for i in range(All_c0.GetNbinsX()):
                        hname = [j for j in fit_dict[fit_var].keys() if 'All' in j][0] ## get 'All' hist corresponding to bin
                        All_c0[i+1].value = fit_dict[fit_var][hname]['%s_g0.9' % fit_degree]['Pars'][0] # set bin value to c0

                    All_c0.Write()

                    All_c1 = Hist(all_xbins, name='All_c1', title='All c1')
                    for i in range(All_c1.GetNbinsX()):
                        hname = [j for j in fit_dict[fit_var].keys() if 'All' in j][0] ## get 'All' hist corresponding to bin
                        All_c1[i+1].value = fit_dict[fit_var][hname]['%s_g0.9' % fit_degree]['Pars'][1] # set bin value to c1

                    All_c1.Write()

                    All_c2 = Hist(all_xbins, name='All_c2', title='All c2')
                    for i in range(All_c2.GetNbinsX()):
                        hname = [j for j in fit_dict[fit_var].keys() if 'All' in j][0] ## get 'All' hist corresponding to bin
                        All_c2[i+1].value = fit_dict[fit_var][hname]['%s_g0.9' % fit_degree]['Pars'][2] # set bin value to c2

                    All_c2.Write()
                #    set_trace()


##############################################################################################
'''
These are plots whose permutations are made from 3-jet events (at least 2 pass b-tagging reqs)
and which have at least one solution from using the merged- and lost-jet assumption discriminants.

solutions_dict-either only one solution exists (Only_Merged(Lost)) from the discriminants,
    both solutions exist but give the same classification (Both_BP/Class_Merged(Lost)),
    or all distributions in which a discriminant has a solution (Merged(Lost)_BP)

both_sols_dict- classifications for events in which both discriminants have solutions (some of the possibilities)
    1. disctributions for both classifications for a discriminant assumption (Merged(Lost)_BP/Class_Merged(Lost))
    2. discr solutions give opposite classifications (Opposite_Class/Merged(Lost)_BP/Class...)
'''



##############################################################################################

def Gen_Plots(directory, subdir):

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]
    rebin_hist = {'Costh' : {'THad' : (-1., 1., 26), 'TLep' : (-1., 1., 26)},\
                  'Eta' : {'TTbar' : (-2.5, 2.5, 26), 'THad' : (-2.5, 2.5, 26), 'TLep' : (-2.5, 2.5, 26)},\
                  'Mass' : {'TTbar' : (200., 2000., 46), 'THad' : (100., 250., 46)},\
                  'Pt' : {'TTbar' : (0., 1000., 51), 'THad' : (0., 1000., 51), 'TLep' : (0., 1000., 51)}}

        ### Gen plots for thad, tlep, and ttbar objects looking at kinvar dists (mass, costh, ...) based on likelihood value
    for kvar in var_types.keys():
        for obj in var_types[kvar].keys():
            plotter.set_subdir('/'.join([subdir, 'Gen', kvar]))
            xlabel = 'Gen %s' % var_types[kvar][obj]
            hname = '/'.join([directory, 'Gen', kvar, obj])

            hist = asrootpy(myfile.Get(hname)).Clone()

            if hist.Integral() == 0:
                continue

            #new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], rebin_hist[kvar][obj][2])
            new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], hist.nbins()/2+1)
            hist = RebinView.rebin(hist, new_bins)
            #hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

            hist_mean = hist.GetMean()
            hist_rms = hist.GetRMS()

            plotter.set_histo_style(hist, xtitle=xlabel, ytitle=defyax)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hist_mean, hist_rms), position='NE')
            box1.Draw()
            plotter.save('Gen_%s_%s' % (obj, kvar))


                ### create hists based on perm category (Correct b, wrong b, etc...)
            to_draw = []
            for cat in Correct_WJet_Categories+Wrong_WJet_Categories+['NO_MP']:
                plotter.set_subdir('/'.join([subdir, 'Gen', kvar, 'Perm_Categories']))
                #set_trace()
                evt_col=styles[cat]['fillcolor']
                evt_type = styles[cat]['name']
                cat_hname = '/'.join([directory, 'Gen', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                if cat_hist.Integral() == 0:
                    continue
    
                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)
    
            if not to_draw:
                continue
            stack, norm_stack, ratio = plotter.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Gen_%s_%s_Stack' % (obj, kvar))
    
            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Gen_%s_%s_Stack_Norm' % (obj, kvar))


                ### create hists based on correct wjet perm assignment
            to_draw = []
            for cat in Correct_WJet_Categories:
                plotter.set_subdir('/'.join([subdir, 'Gen', kvar, 'Perm_Categories', 'Correct_WJet']))
                #set_trace()
                evt_col=styles[cat]['fillcolor']
                evt_type = styles[cat]['name']
                cat_hname = '/'.join([directory, 'Gen', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

                if cat_hist.Integral() == 0:
                    continue

                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)

            if not to_draw:
                continue
            stack, norm_stack, ratio = plotter.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            #set_trace()
            hmean = sum([ stack[i] for i in range(len(stack))]).GetMean()
            hrms = sum([ stack[i] for i in range(len(stack))]).GetRMS()
            box2 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
            box2.Draw()
            plotter.save('Gen_%s_%s_Stack' % (obj, kvar))

            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box2.Draw()
            plotter.save('Gen_%s_%s_Stack_Norm' % (obj, kvar))

                ### create hists based on wrong wjet perm assignment
            to_draw = []
            for cat in Wrong_WJet_Categories:
                plotter.set_subdir('/'.join([subdir, 'Gen', kvar, 'Perm_Categories', 'Wrong_WJet']))
                #set_trace()
                evt_col=styles[cat]['fillcolor']
                evt_type = styles[cat]['name']
                cat_hname = '/'.join([directory, 'Gen', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

                if cat_hist.Integral() == 0:
                    continue

                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)

            if not to_draw:
                continue
            stack, norm_stack, ratio = plotter.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            #set_trace()
            hmean = sum([ stack[i] for i in range(len(stack))]).GetMean()
            hrms = sum([ stack[i] for i in range(len(stack))]).GetRMS()
            box2 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
            box2.Draw()
            plotter.save('Gen_%s_%s_Stack' % (obj, kvar))

            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box2.Draw()
            plotter.save('Gen_%s_%s_Stack_Norm' % (obj, kvar))



##############################################################################################


def Reco_Plots(directory, subdir):#, topology):

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]
    rebin_hist = {'Costh' : {'THad' : (-1., 1., 26), 'TLep' : (-1., 1., 26)},\
                  'Eta' : {'TTbar' : (-2.5, 2.5, 26), 'THad' : (-2.5, 2.5, 26), 'TLep' : (-2.5, 2.5, 26)},\
                  'Mass' : {'TTbar' : (200., 2000., 46), 'THad' : (50., 250., 46)},\
                  'Pt' : {'TTbar' : (0., 1000., 51), 'THad' : (0., 1000., 51), 'TLep' : (0., 1000., 51)}}

        ### Reco plots for thad, tlep, and ttbar objects looking at kinvar dists (mass, costh, ...) based on likelihood value
    for kvar in var_types.keys():
        for obj in var_types[kvar].keys():
            plotter.set_subdir('/'.join([subdir, 'Reconstruction', kvar]))
            xlabel = 'Reco %s' % var_types[kvar][obj]
            hname = '/'.join([directory, 'Reconstruction', kvar, obj])

            hist = asrootpy(myfile.Get(hname)).Clone()

            if hist.Integral() == 0:
                continue

            #if kvar == 'Mass':
            #    pre_rebin = [hist.Integral(i+1,i+1) for i in range(hist.nbins())]
            #    #new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], rebin_hist[kvar][obj][2])
            #    post_rebin = [hist.Integral(i+1,i+1) for i in range(hist.nbins())]
            #set_trace()
            new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], hist.nbins()/2+1)
            hist = RebinView.rebin(hist, new_bins)
            #hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

            hist_mean = hist.GetMean()
            hist_rms = hist.GetRMS()

            plotter.set_histo_style(hist, xtitle=xlabel, ytitle=defyax)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hist_mean, hist_rms), position='NE')
            box1.Draw()
            plotter.save('Reco_%s_%s' % (obj, kvar))


                ### create hists based on perm category (Correct b, wrong b, etc...)
            to_draw = []
            for cat in Correct_WJet_Categories+Wrong_WJet_Categories+['NO_MP']:
                plotter.set_subdir('/'.join([subdir, 'Reconstruction', kvar, 'Perm_Categories']))
                #set_trace()
                evt_col=styles[cat]['fillcolor']
                evt_type = styles[cat]['name']
                cat_hname = '/'.join([directory, 'Reconstruction', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                if cat_hist.Integral() == 0:
                    continue
    
                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)
    
            if not to_draw:
                continue
            stack, norm_stack, ratio = plotter.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Reco_%s_%s_Stack' % (obj, kvar))
    
            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Reco_%s_%s_Stack_Norm' % (obj, kvar))


                ### create hists based on correct wjet perm assignment
            to_draw = []
            for cat in Correct_WJet_Categories:
                plotter.set_subdir('/'.join([subdir, 'Reconstruction', kvar, 'Perm_Categories', 'Correct_WJet']))
                #set_trace()
                evt_col=styles[cat]['fillcolor']
                evt_type = styles[cat]['name']
                cat_hname = '/'.join([directory, 'Reconstruction', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

                if cat_hist.Integral() == 0:
                    continue

                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)

            if not to_draw:
                continue
            stack, norm_stack, ratio = plotter.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            #set_trace()
            hmean = sum([ stack[i] for i in range(len(stack))]).GetMean()
            hrms = sum([ stack[i] for i in range(len(stack))]).GetRMS()
            box2 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
            box2.Draw()
            plotter.save('Reco_%s_%s_Stack' % (obj, kvar))

            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box2.Draw()
            plotter.save('Reco_%s_%s_Stack_Norm' % (obj, kvar))


                ### create hists based on wrong wjet perm assignment
            to_draw = []
            for cat in Wrong_WJet_Categories:
                plotter.set_subdir('/'.join([subdir, 'Reconstruction', kvar, 'Perm_Categories', 'Wrong_WJet']))
                #set_trace()
                evt_col=styles[cat]['fillcolor']
                evt_type = styles[cat]['name']
                cat_hname = '/'.join([directory, 'Reconstruction', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

                if cat_hist.Integral() == 0:
                    continue

                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)

            if not to_draw:
                continue
            stack, norm_stack, ratio = plotter.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            #set_trace()
            hmean = sum([ stack[i] for i in range(len(stack))]).GetMean()
            hrms = sum([ stack[i] for i in range(len(stack))]).GetRMS()
            box2 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
            box2.Draw()
            plotter.save('Reco_%s_%s_Stack' % (obj, kvar))

            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box2.Draw()
            plotter.save('Reco_%s_%s_Stack_Norm' % (obj, kvar))


##############################################################################################


def Resolution_Plots(directory, subdir):#, topology):

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]
    rebin_hist = {'Costh' : {'TTbar' : (-2., 2., 26), 'THad' : (-2., 2., 26), 'TLep' : (-2., 2., 26)},\
                  'Eta' : {'TTbar' : (-5.0, 50, 26), 'THad' : (-5.0, 5.0, 26), 'TLep' : (-5.0, 5.0, 26)},\
                  'Mass' : {'TTbar' : (-1000., 2000., 46), 'THad' : (-1000., 500., 46), 'TLep' : (-500., 500., 46)},\
                  'Pt' : {'TTbar' : (-1000., 500., 51), 'THad' : (-500., 500., 51), 'TLep' : (-1000., 1000., 51)}}

        ### Resolution plots for thad, tlep, and ttbar objects looking at kinvar dists (mass, costh, ...) based on likelihood value
    for kvar in var_types.keys():
        for obj in var_types[kvar].keys():
            plotter.set_subdir('/'.join([subdir, 'Resolution', kvar]))
            xlabel = '%s Resolution (Gen-Reco)' % var_types[kvar][obj]
            hname = '/'.join([directory, 'Resolution', kvar, obj])

            hist = asrootpy(myfile.Get(hname)).Clone()

            if hist.Integral() == 0:
                continue

            #new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], rebin_hist[kvar][obj][2])
            new_bins = np.linspace(rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1], hist.nbins()/2+1)
            hist = RebinView.rebin(hist, new_bins)
            #hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

            hist_mean = hist.GetMean()
            hist_rms = hist.GetRMS()

            plotter.set_histo_style(hist, xtitle=xlabel, ytitle=defyax)
            plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hist_mean, hist_rms), position='NE')
            box1.Draw()
            plotter.save('Reso_%s_%s' % (obj, kvar))



                ### create hists based on perm category (Correct b, wrong b, etc...)
            to_draw = []
            for cat in Correct_WJet_Categories+Wrong_WJet_Categories+['NO_MP']:
                plotter.set_subdir('/'.join([subdir, 'Resolution', kvar, 'Perm_Categories']))
                #set_trace()
                evt_col=styles[cat]['fillcolor']
                evt_type = styles[cat]['name']
                cat_hname = '/'.join([directory, 'Resolution', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                if cat_hist.Integral() == 0:
                    continue
    
                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)
    
            if not to_draw:
                continue
            stack, norm_stack, ratio = plotter.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Reso_%s_%s_Stack' % (obj, kvar))
    
            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Reso_%s_%s_Stack_Norm' % (obj, kvar))


                ### create hists based on correct wjet perm assignment
            to_draw = []
            for cat in Correct_WJet_Categories:
                plotter.set_subdir('/'.join([subdir, 'Resolution', kvar, 'Perm_Categories', 'Correct_WJet']))
                #set_trace()
                evt_col=styles[cat]['fillcolor']
                evt_type = styles[cat]['name']
                cat_hname = '/'.join([directory, 'Resolution', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

                if cat_hist.Integral() == 0:
                    continue

                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)

            if not to_draw:
                continue
            stack, norm_stack, ratio = plotter.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            #set_trace()
            hmean = sum([ stack[i] for i in range(len(stack))]).GetMean()
            hrms = sum([ stack[i] for i in range(len(stack))]).GetRMS()
            box2 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
            box2.Draw()
            plotter.save('Reso_%s_%s_Stack' % (obj, kvar))

            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box2.Draw()
            plotter.save('Reso_%s_%s_Stack_Norm' % (obj, kvar))


                ### create hists based on wrong wjet perm assignment
            to_draw = []
            for cat in Wrong_WJet_Categories:
                plotter.set_subdir('/'.join([subdir, 'Resolution', kvar, 'Perm_Categories', 'Wrong_WJet']))
                #set_trace()
                evt_col=styles[cat]['fillcolor']
                evt_type = styles[cat]['name']
                cat_hname = '/'.join([directory, 'Resolution', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

                if cat_hist.Integral() == 0:
                    continue

                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)

            if not to_draw:
                continue
            stack, norm_stack, ratio = plotter.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            #set_trace()
            hmean = sum([ stack[i] for i in range(len(stack))]).GetMean()
            hrms = sum([ stack[i] for i in range(len(stack))]).GetRMS()
            box2 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
            box2.Draw()
            plotter.save('Reso_%s_%s_Stack' % (obj, kvar))

            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box2.Draw()
            plotter.save('Reso_%s_%s_Stack_Norm' % (obj, kvar))


    
#####################################################################################################


def Discriminant_Plots(directory, subdir):#, topology):

        ### Discriminant plots for Mass, NS, and Combined discs looking at evt_type (RIGHT, WRONG, ...)

    plotter.defaults['watermark'] = ['%s %s (13 TeV, 25ns)' % (decay, m_range), False]
    #plotter.set_subdir('/'.join([subdir, 'Discr']))
    for disc in DiscPlots:
        plotter.set_subdir('/'.join([subdir, 'Discr']))
        xlabel = DiscPlots[disc][0]
        hname = '/'.join([directory, 'Discr', disc])

        hist = asrootpy(myfile.Get(hname)).Clone()
        if hist.Integral() == 0:
            continue
        hist_mean = hist.GetMean()
        hist_rms = hist.GetRMS()

        plotter.set_histo_style(hist, xtitle=xlabel+' 3 jets', ytitle=defyax)
        plotter.plot(hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
        box1 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hist_mean, hist_rms), position='NE')
        box1.Draw()
        plotter.save(disc)


                ### create hists based on perm category (Correct b, wrong b, etc...)
        to_draw = []
        for cat in Correct_WJet_Categories+Wrong_WJet_Categories+['NO_MP']:
            plotter.set_subdir('/'.join([subdir, 'Discr', disc, 'Perm_Categories']))
            #set_trace()
            evt_col=styles[cat]['fillcolor']
            evt_type = styles[cat]['name']
            cat_hname = '/'.join([directory, 'Discr', cat, disc])
            cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
        
            #cat_hist = RebinView.rebin(cat_hist, new_bins)
            #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
        
            if cat_hist.Integral() == 0:
                continue
        
            plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel+' 3 jets', ytitle=defyax)
            cat_hist.SetFillStyle(styles[cat]['fillstyle'])
            plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(cat_hist)
        
        if not to_draw:
            continue
        stack, norm_stack, ratio = plotter.stack_plots(to_draw)
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        box1.Draw()
        plotter.save(disc+'_Stack')
        
        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        box1.Draw()
        plotter.save(disc+'_Stack_Norm')
        

            ### create hists based on correct wjet perm assignment
        to_draw = []
        for cat in Correct_WJet_Categories:
            plotter.set_subdir('/'.join([subdir, 'Discr', disc, 'Perm_Categories', 'Correct_WJet']))
            #set_trace()
            evt_col=styles[cat]['fillcolor']
            evt_type = styles[cat]['name']
            cat_hname = '/'.join([directory, 'Discr', cat, disc])
            cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

            #cat_hist = RebinView.rebin(cat_hist, new_bins)
            #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

            if cat_hist.Integral() == 0:
                continue

            plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel+' 3 jets', ytitle=defyax)
            cat_hist.SetFillStyle(styles[cat]['fillstyle'])
            plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(cat_hist)

        if not to_draw:
            continue
        stack, norm_stack, ratio = plotter.stack_plots(to_draw)
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        #set_trace()
        hmean = sum([ stack[i] for i in range(len(stack))]).GetMean()
        hrms = sum([ stack[i] for i in range(len(stack))]).GetRMS()
        box2 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
        box2.Draw()
        plotter.save(disc+'_Stack')

        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        box2.Draw()
        plotter.save(disc+'_Stack_Norm')

            ### create hists based on wrong wjet perm assignment
        to_draw = []
        for cat in Wrong_WJet_Categories:
            plotter.set_subdir('/'.join([subdir, 'Discr', disc, 'Perm_Categories', 'Wrong_WJet']))
            #set_trace()
            evt_col=styles[cat]['fillcolor']
            evt_type = styles[cat]['name']
            cat_hname = '/'.join([directory, 'Discr', cat, disc])
            cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()

            #cat_hist = RebinView.rebin(cat_hist, new_bins)
            #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]

            if cat_hist.Integral() == 0:
                continue

            plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel+' 3 jets', ytitle=defyax)
            cat_hist.SetFillStyle(styles[cat]['fillstyle'])
            plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(cat_hist)

        if not to_draw:
            continue
        stack, norm_stack, ratio = plotter.stack_plots(to_draw)
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        #set_trace()
        hmean = sum([ stack[i] for i in range(len(stack))]).GetMean()
        hrms = sum([ stack[i] for i in range(len(stack))]).GetRMS()
        box2 = plotter.make_text_box('Mean=%.2f\nRMS=%.2f' % (hmean, hrms), position='NE')
        box2.Draw()
        plotter.save(disc+'_Stack')

        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        box2.Draw()
        plotter.save(disc+'_Stack_Norm')
        
    

#####################################################################################################

def Final_Reco_Plots( plot ):
    #set_trace()
    reco_dir = '3J/nosys'
    reco_subdir = '3J_Event_Plots/Lost_BP'

    if plot == 'Discriminant':
        print '\nMaking Discr plots for 3-jet events\n\n'
        Discriminant_Plots( reco_dir , reco_subdir ) 

    if plot == 'Gen':
        print '\nMaking Gen plots for 3-jet events\n\n' 
        Gen_Plots( reco_dir , reco_subdir )

    if plot == 'Reconstruction':
        print '\nMaking reco plots for 3-jet events\n\n'
        Reco_Plots( reco_dir , reco_subdir)

    if plot == 'Resolution':
        print '\nMaking resolution plots for 3-jet events\n\n'
        Resolution_Plots( reco_dir , reco_subdir)

    if plot == 'Alpha_Correction':
        print '\nMaking alpha correction plots for 3-jet events\n\n'
        alpha_corrections( reco_dir , reco_subdir )



if args.plot == 'All':
    for plot_type in allowed_plot_types:
        Final_Reco_Plots(plot_type)
else:
    Final_Reco_Plots(args.plot)

