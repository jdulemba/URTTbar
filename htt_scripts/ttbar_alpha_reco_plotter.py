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
import functions as fncts
import numpy as np
from styles import styles
import URAnalysis.Utilities.prettyjson as prettyjson


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

hist_styles = {
            'CORRECT_B': {
                'legendstyle' : 'l',
                'drawstyle' : 'hist',
                'fillcolor' : 'red',
                'linecolor' : 'black',
                'linewidth' : 1,
                'name' : "Correct b's",
                'fillstyle': '3345',
            },
            'WRONG_B': {
                'legendstyle' : 'l',
                'drawstyle' : 'hist',
                'fillcolor' : 'blue',
                'linecolor' : 'black',
                'linewidth' : 1,
                'name' : "Wrong b's",
                'fillstyle': '3354',
            },
            'OTHER' : styles['*OTHER'],

            'CORRECT_BHAD': {
                'legendstyle' : 'l',
                'drawstyle' : 'hist',
                'fillcolor' : 'blue',
                'linecolor' : 'black',
                'linewidth' : 1,
                'name' : "Correct b_{h} match",
                'fillstyle': '3345',
            },
            'CORRECT_BLEP': {
                'legendstyle' : 'l',
                'drawstyle' : 'hist',
                'fillcolor' : 'red',
                'linecolor' : 'black',
                'linewidth' : 1,
                'name' : "Correct b_{l} match",
                'fillstyle': '3354',
            },
            'CORRECT_Bs': {
                'legendstyle' : 'l',
                'drawstyle' : 'hist',
                'fillcolor' : 'black',
                'linecolor' : 'black',
                'linewidth' : 1,
                'name' : "Correct b's",
                'fillstyle': '0',
            },
            'SWAPPED_Bs': styles['*SWAP'],
            'OTHER_MATCH' : styles['*OTHER'],
}



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
            'RIGHT' : styles['*RIGHT'],
            'MERGED_SWAP' : styles['*SWAP'],
            'MERGED' : styles['*OTHER'],
            'WRONG' : styles['*WRONG'],
            'LOST_SWAP' : styles['*SWAP'],
            'LOST' : styles['*OTHER'],
            'sample' : styles[args.sample]
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
            'RIGHT' : styles['*RIGHT'],
            'MERGED_SWAP' : styles['*SWAP'],
            'MERGED' : styles['*OTHER'],
            'WRONG' : styles['*WRONG'],
            'LOST_SWAP' : styles['*SWAP'],
            'LOST' : styles['*OTHER'],
            'sample' : styles['ttJetsM0']
            #'sample' : styles[args.sample]
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

Perm_Categories = ['CORRECT_B', 'WRONG_B', 'OTHER']
Gen_Categories = ['CORRECT_BHAD', 'CORRECT_BLEP', 'CORRECT_Bs', 'SWAPPED_Bs', 'OTHER_MATCH']


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


    hvars = [ # variables for hists
        ('THad_E/Alpha_THad_E', '173.1/Reco M(t_{h})', '#alpha_{E} = Gen E(t_{h})/Reco E(t_{h})', ''),
        ('THad_P/Alpha_THad_P', '173.1/Reco M(t_{h})', '#alpha_{P} = Gen P(t_{h})/Reco P(t_{h})', ''),
    ]

    xranges = np.linspace(0.9, 2.5, 9) ## rebin xaxis so [0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
    for num in range(len(xranges)-1):
        txt_box = "%.1f < 173.1/Reco M(t_{h}) #leq %.1f" % (xranges[num], xranges[num+1])
        Evar = ('THad_E/Gen_vs_Reco_THadE_%.1fto%.1f' % (xranges[num], xranges[num+1]), 'Reco E(t_{h})', 'Gen E(t_{h})', txt_box )
        Pvar = ('THad_P/Gen_vs_Reco_THadP_%.1fto%.1f' % (xranges[num], xranges[num+1]), 'Reco E(t_{h})', 'Gen E(t_{h})', txt_box )
        hvars.append(Evar)
        hvars.append(Pvar)

    fits = {
        'Alpha_THad_E' : {},
        'Alpha_THad_P' : {},
    }

    for hvar, xlabel, ylabel, txt_box_label in hvars:
        hname = '/'.join([directory, 'Alpha_Correction', hvar])
        hist = asrootpy(myfile.Get(hname)).Clone()

        if hist.Integral() == 0:
            continue

        plotter.plot(hist)
        #set_trace()
        if hvar == 'THad_E/Alpha_THad_E' or hvar == 'THad_P/Alpha_THad_P':

            xbins = np.linspace(0.9, 2.5, 9) ## rebin xaxis so [0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5]
            ybins = np.linspace( hist.GetYaxis().GetBinLowEdge(1), hist.GetYaxis().GetBinUpEdge(hist.GetYaxis().GetNbins()), hist.GetYaxis().GetNbins()+1 ) # ybinning remains unchanged

            hist = RebinView.newRebin2D(hist, xbins, ybins)
            hist.xaxis.range_user = 0.9, 2.5
            #hist.yaxis.range_user = 0.0, 4.0
            #set_trace()

            medians = []
            median_weights = []
            median_errors = []
            maxvals = []
            for xbin in range(hist.GetNbinsX() + 1):
                if hist.GetXaxis().GetBinLowEdge(xbin) >= hist.xaxis.range_user[0] and hist.GetXaxis().GetBinUpEdge(xbin) <= hist.xaxis.range_user[1]:
                    hist_yproj = hist.ProjectionY("", xbin, xbin)
                    if hist_yproj.Integral() == 0:
                        continue

                    hist_yproj.Draw()
                    hist_yproj.GetXaxis().SetRangeUser(0.0,15.0)
                    hist_yproj.GetXaxis().SetTitle(ylabel)

                    probs = np.zeros(len(ybins))
                    binvals = np.zeros(len(ybins))

                    for bins in range(len(ybins)):
                        probs[bins] = hist_yproj.Integral(1, bins+1)/hist_yproj.Integral()
                        binvals[bins] = hist_yproj.GetBinContent(bins)

                    median = hist_yproj.GetBinCenter(np.where(probs>= 0.5)[0][0])
                    median_error = 1.2533*hist_yproj.GetMeanError() #standard error of median
                    if median_error == 0:
                        #set_trace()
                        median_weight = 1000 # 1/(standard error of median) to be used in fit
                    else:
                        median_weight = 1/(1.2533*hist_yproj.GetMeanError()) # 1/(standard error of median) to be used in fit
                    maxval = round(hist_yproj.GetBinCenter(np.where(binvals==binvals.max())[0][0]), 2)

                    medians.append(median)
                    median_weights.append(median_weight)
                    median_errors.append(median_error)
                    maxvals.append(maxval)

                    box1 = plotter.make_text_box('%s #leq %s #leq %s' % (hist.GetXaxis().GetBinLowEdge(xbin), xlabel, hist.GetXaxis().GetBinUpEdge(xbin)), position='NE')
                    box1.Draw()
                    plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', hvar.split('/')[0], 'y_proj']))
                    plotter.save('%s_%s' % (hvar.split('/')[1], hist.GetXaxis().GetBinLowEdge(xbin)) )
                    #set_trace()

            xbin_ints = [hist.Integral(xbinid+1, xbinid+1, 1, hist.GetNbinsY()+1) for xbinid in range(hist.GetNbinsX())]
            xbin_ratios = [xbin_ints[i]/sum(xbin_ints) for i in range(len(xbin_ints))]

            hist.yaxis.range_user = 0.0, 4.0

            #### section to fit polynomials to median values from xaxis bins
            #xfit_bins = np.linspace(hist.GetXaxis().GetXmin(),hist.GetXaxis().GetXmax(), 200) #makes points to be used for fit results
            xfit_bins = np.linspace(0, 5.0, 200) #makes points to be used for fit results
            hist.GetXaxis().GetCenter(xbins) # change array of xbin values from edges to centers

            fit_1d_groom = np.polyfit(xbins[:-1], medians, 1, w=median_weights) # fit medians to 1 degree polynomial for 0.9 < xbin < 2.5
            fit_1d_1p1 = np.polyfit(xbins[1:-1], medians[1:], 1, w=median_weights[1:]) # weighted fit medians to 1 degree polynomial for 1.1 < xbin < 2.5
            fit_2d_groom = np.polyfit(xbins[:-1], medians, 2, w=median_weights) # fit medians to 2 degree polynomial for 0.9 < xbin < 2.5
            fit_2d_1p1 = np.polyfit(xbins[1:-1], medians[1:], 2, w=median_weights[1:]) # fit medians to 2 degree polynomial for 1.1 < xbin < 2.5

            p1_groom = np.poly1d(fit_1d_groom) #gets parameters for 1 degree fit getting rid of first 2 points
            p1_1p1 = np.poly1d(fit_1d_1p1) #gets parameters for weighted 1 degree fit getting rid of first 3 points
            p2_groom = np.poly1d(fit_2d_groom) #gets parameters for 2 degree fit getting rid of first 2 points
            p2_1p1 = np.poly1d(fit_2d_1p1) #gets parameters for weighted 2 degree fit getting rid of first 3 points

            fig = plt.figure()
            fit_medians = plt.errorbar(xbins[:-1], medians, yerr=median_errors, label='Medians', linestyle='None', color='black', marker='.') #plots median values with errors
                # 1 degree fits
            p1_1p1_weighted = plt.plot(xfit_bins, p1_1p1(xfit_bins), label='weighted deg=1 bins > 1.1', linestyle='-', color='red') #plots p1 for weighted bins > 1.1
            p1_groom_fit = plt.plot(xfit_bins, p1_groom(xfit_bins), label='weighted deg=1 bins > 0.9', linestyle='--', color='red') #plots p1 for bins > 0.9
                # 2 degree fits
            p2_1p1_weighted = plt.plot(xfit_bins, p2_1p1(xfit_bins), label='weighted deg=2 bins > 1.1', linestyle='-', color='green') #plots p2 for weighted bins > 1.1
            p2_groom_fit = plt.plot(xfit_bins, p2_groom(xfit_bins), label='weighted deg=2 bins > 0.9', linestyle='--', color='green') #plots p2 for bins > 0.9

            fits[hvar.split('/')[1]]['1D_g1.1'] = fit_1d_1p1.tolist()
            fits[hvar.split('/')[1]]['2D_g1.1'] = fit_2d_1p1.tolist()

            plt.xlabel('173.1/Reco M($t_{h}$)')
            if hvar == 'THad_E/Alpha_THad_E':
                plt.ylabel('$\\alpha_{E}$ = Gen E($t_{h}$)/Reco E($t_{h}$)')
            if hvar == 'THad_P/Alpha_THad_P':
               plt.ylabel('$\\alpha_{P}$ = Gen P($t_{h}$)/Reco P($t_{h}$)')
            plt.ylim(0.0, 4.0)
            plt.grid()
            plt.legend(loc='upper right',fontsize=8, numpoints=1)
            plotter.set_subdir('/'.join([subdir, 'Alpha_Correction']))
            plotting_dir = plotter.outputdir
            fig.savefig('%s/%s_fits' % (plotting_dir, hvar.split('/')[1]))
            #set_trace()

        plotter.set_histo_style(hist, xtitle=xlabel, ytitle=ylabel)
        plotter.plot(hist, drawstyle='colz')

        mean_x = hist.ProjectionX().GetMean()
        rms_x = hist.ProjectionX().GetRMS()
        mean_y = hist.ProjectionY().GetMean()
        rms_y = hist.ProjectionY().GetRMS()

        box2 = plotter.make_text_box('X Mean=%.2f, RMS=%.2f\nY Mean=%.2f, RMS=%.2f' % (mean_x, rms_x, mean_y, rms_y), position='NE')
        box2.Draw()

        if 'Gen_vs_Reco' in hvar:
            box1 = plotter.make_text_box( txt_box_label, position='NE')
            box1.Draw()
            if 'THadE' in hvar:
                plotter.set_subdir('/'.join([subdir, 'Alpha_Correction/THad_E/Gen_vs_Reco']))
            if 'THadP' in hvar:
                plotter.set_subdir('/'.join([subdir, 'Alpha_Correction/THad_P/Gen_vs_Reco']))
        else:
            plotter.set_subdir('/'.join([subdir, 'Alpha_Correction']))

        plotter.save(hvar.split('/')[1])
        
        #set_trace()

        #### create hists based on perm category (RIGHT/WRONG, etc...)
        #if topology:
        #    for cat in Categories[topology]:
        #        plotter.set_subdir('/'.join([subdir, 'Alpha_Correction', cat]))
        #        #evt_col=hist_styles[cat]['fillcolor']
        #        evt_type = hist_styles[cat]['name']
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

    with open('%s/fit_parameters.json' % '/'.join([plotter.base_out_dir, subdir, 'Alpha_Correction']), 'w') as f:
        #set_trace()
        f.write(prettyjson.dumps(fits))



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
            for cat in Perm_Categories:#[topology]:
                plotter.set_subdir('/'.join([subdir, 'Gen', kvar, 'Perm_Categories']))
                #set_trace()
                evt_col=hist_styles[cat]['fillcolor']
                evt_type = hist_styles[cat]['name']
                cat_hname = '/'.join([directory, 'Gen', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                if cat_hist.Integral() == 0:
                    continue
    
                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(hist_styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)
    
            if not to_draw:
                continue
            stack, norm_stack, ratio = fncts.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Gen_%s_%s_Stack' % (obj, kvar))
    
            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Gen_%s_%s_Stack_Norm' % (obj, kvar))


                ### create hists based on gen category (Correct bhad, correct blep, wrong b, etc...)
            to_draw = []
            for cat in Gen_Categories:#[topology]:
                plotter.set_subdir('/'.join([subdir, 'Gen', kvar, 'Gen_Categories']))
                #set_trace()
                evt_col=hist_styles[cat]['fillcolor']
                evt_type = hist_styles[cat]['name']
                cat_hname = '/'.join([directory, 'Gen', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                if cat_hist.Integral() == 0:
                    continue
    
                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(hist_styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)
    
            if not to_draw:
                continue
            stack, norm_stack, ratio = fncts.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Gen_%s_%s_Stack' % (obj, kvar))
    
            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
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
            for cat in Perm_Categories:#[topology]:
                plotter.set_subdir('/'.join([subdir, 'Reconstruction', kvar, 'Perm_Categories']))
                #set_trace()
                evt_col=hist_styles[cat]['fillcolor']
                evt_type = hist_styles[cat]['name']
                cat_hname = '/'.join([directory, 'Reconstruction', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                if cat_hist.Integral() == 0:
                    continue
    
                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(hist_styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)
    
            if not to_draw:
                continue
            stack, norm_stack, ratio = fncts.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Reco_%s_%s_Stack' % (obj, kvar))
    
            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Reco_%s_%s_Stack_Norm' % (obj, kvar))


                ### create hists based on gen category (Correct bhad, correct blep, wrong b, etc...)
            to_draw = []
            for cat in Gen_Categories:#[topology]:
                plotter.set_subdir('/'.join([subdir, 'Reconstruction', kvar, 'Gen_Categories']))
                #set_trace()
                evt_col=hist_styles[cat]['fillcolor']
                evt_type = hist_styles[cat]['name']
                cat_hname = '/'.join([directory, 'Reconstruction', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                if cat_hist.Integral() == 0:
                    continue
    
                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(hist_styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)
    
            if not to_draw:
                continue
            stack, norm_stack, ratio = fncts.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Reco_%s_%s_Stack' % (obj, kvar))
    
            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
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
            for cat in Perm_Categories:#[topology]:
                plotter.set_subdir('/'.join([subdir, 'Resolution', kvar, 'Perm_Categories']))
                #set_trace()
                evt_col=hist_styles[cat]['fillcolor']
                evt_type = hist_styles[cat]['name']
                cat_hname = '/'.join([directory, 'Resolution', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                if cat_hist.Integral() == 0:
                    continue
    
                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(hist_styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)
    
            if not to_draw:
                continue
            stack, norm_stack, ratio = fncts.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Reso_%s_%s_Stack' % (obj, kvar))
    
            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Reso_%s_%s_Stack_Norm' % (obj, kvar))


                ### create hists based on gen category (Correct bhad, correct blep, wrong b, etc...)
            to_draw = []
            for cat in Gen_Categories:#[topology]:
                plotter.set_subdir('/'.join([subdir, 'Resolution', kvar, 'Gen_Categories']))
                #set_trace()
                evt_col=hist_styles[cat]['fillcolor']
                evt_type = hist_styles[cat]['name']
                cat_hname = '/'.join([directory, 'Resolution', cat, kvar, obj])
                cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
    
                cat_hist = RebinView.rebin(cat_hist, new_bins)
                #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
    
                if cat_hist.Integral() == 0:
                    continue
    
                plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel, ytitle=defyax)
                cat_hist.SetFillStyle(hist_styles[cat]['fillstyle'])
                plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
                to_draw.append(cat_hist)
    
            if not to_draw:
                continue
            stack, norm_stack, ratio = fncts.stack_plots(to_draw)
            plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
            plotter.save('Reso_%s_%s_Stack' % (obj, kvar))
    
            plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel, ytitle=defyax, drawstyle='hist')
            box1.Draw()
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
        for cat in Perm_Categories:#[topology]:
            plotter.set_subdir('/'.join([subdir, 'Discr', disc, 'Perm_Categories']))
            #set_trace()
            evt_col=hist_styles[cat]['fillcolor']
            evt_type = hist_styles[cat]['name']
            cat_hname = '/'.join([directory, 'Discr', cat, disc])
            cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
        
            #cat_hist = RebinView.rebin(cat_hist, new_bins)
            #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
        
            if cat_hist.Integral() == 0:
                continue
        
            plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel+' 3 jets', ytitle=defyax)
            cat_hist.SetFillStyle(hist_styles[cat]['fillstyle'])
            plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(cat_hist)
        
        if not to_draw:
            continue
        stack, norm_stack, ratio = fncts.stack_plots(to_draw)
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        box1.Draw()
        plotter.save(disc+'_Stack')
        
        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        box1.Draw()
        plotter.save(disc+'_Stack_Norm')
        
        
                 ### create hists based on gen category (Correct bhad, correct blep, wrong b, etc...)
        to_draw = []
        for cat in Gen_Categories:#[topology]:
            plotter.set_subdir('/'.join([subdir, 'Discr', disc, 'Gen_Categories']))
            #set_trace()
            evt_col=hist_styles[cat]['fillcolor']
            evt_type = hist_styles[cat]['name']
            cat_hname = '/'.join([directory, 'Discr', cat, disc])
            cat_hist = asrootpy(myfile.Get(cat_hname)).Clone()
        
            #cat_hist = RebinView.rebin(cat_hist, new_bins)
            #cat_hist.xaxis.range_user = rebin_hist[kvar][obj][0], rebin_hist[kvar][obj][1]
        
            if cat_hist.Integral() == 0:
                continue
        
            plotter.set_histo_style(cat_hist, color=evt_col, title=evt_type, xtitle=xlabel+' 3 jets', ytitle=defyax)
            cat_hist.SetFillStyle(hist_styles[cat]['fillstyle'])
            plotter.plot(cat_hist, legend_def=LegendDefinition(position='NW'), legendstyle='l', drawstyle='hist')
            to_draw.append(cat_hist)
        
        if not to_draw:
            continue
        stack, norm_stack, ratio = fncts.stack_plots(to_draw)
        plotter.plot(stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        box1.Draw()
        plotter.save(disc+'_Stack')
        
        plotter.plot(norm_stack, legend_def=LegendDefinition(position='NW'), legendstyle='l', xtitle=xlabel+' 3 jets', ytitle=defyax, drawstyle='hist')
        box1.Draw()
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

