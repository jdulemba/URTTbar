import rootpy.plotting.views as views
import rootpy.io as io
import os, glob
import rootpy
from URAnalysis.PlotTools.data_views import extract_sample
from pdb import set_trace
import ROOT
rootpy.log["/"].setLevel(rootpy.log.INFO)
log = rootpy.log["/make_permutation_distros.py"]
from URAnalysis.PlotTools.BasePlotter import BasePlotter
from argparse import ArgumentParser
from rootpy import asrootpy
import numpy as np
from URAnalysis.PlotTools.views.RebinView import RebinView

parser = ArgumentParser()
parser.add_argument('out', help='output file name')
#parser.add_argument('--pdfs', action='store_true', help='make plots for the PDF uncertainties')
args = parser.parse_args()


#configuration
### Joseph added for 3J events (merged and lost jets)
Merged_shapes = [
    ('mbpjet_vs_maxmjet', 'm_{j}^{max} 3 jets [GeV]', 'm_{j_{i}+j_{j}} [GeV]', 'A.U.'),
    ('nusolver_chi2', '#chi^{2} 3 jets', 'A.U.', ''),
    ('nusolver_dist', 'D_{#nu, min} 3 jets', 'A.U.', '')
]

Lost_shapes = [
    ('mbpjet', 'm_{j_{i}+j_{j}} [Gev]', 'A.U', ''),
    ('nusolver_chi2', '#chi^{2} 3 jets', 'A.U.', ''),
    ('nusolver_dist', 'D_{#nu, min} 3 jets', 'A.U.', '')
]

shapes_3J = [
    ('mbpjet_vs_maxmjet', 'm_{j}^{max} 3 jets [GeV]', 'm_{j_{i}+j_{j}} [GeV]', 'A.U.'),
    ('mbpjet', 'm_{j_{i}+j_{j}} [Gev]', 'A.U', ''),
    ('nusolver_chi2', '#chi^{2} 3 jets', 'A.U.', ''),
    ('nusolver_dist', 'D_{#nu, min} 3 jets', 'A.U.', '')
]


right_3J = ['3J/MERGED/Correct_Plots', '3J/LOST/Correct_Plots']
wrong_3J = ['3J/MERGED/Wrong_Plots', '3J/LOST/Wrong_Plots']
###
### resolved jet cases (same as make_perm_distros
shapes = [
    ('mWhad_vs_mtophad', 'M(W_{had})', 'M(t_{had})', 'A.U.'),
    ('mWhad_vs_mtophad_4J', 'M(W_{had})', 'M(t_{had})', 'A.U.'),
    ('mWhad_vs_mtophad_5J', 'M(W_{had})', 'M(t_{had})', 'A.U.'),
    ('mWhad_vs_mtophad_6PJ', 'M(W_{had})', 'M(t_{had})', 'A.U.'),
    #('Rel_Delta_mWhad_vs_Rel_Delta_mtophad', '#Delta M(W_{had})/#sigma(2J)', '#Delta M(t_{had})/#sigma(3J)', 'A.U.'),
    #('Rel_Delta_mWhad_vs_Rel_Delta_mtophad_4J', '#Delta M(W_{had})/#sigma(2J)', '#Delta M(t_{had})/#sigma(3J)', 'A.U.'),
    #('Rel_Delta_mWhad_vs_Rel_Delta_mtophad_5J', '#Delta M(W_{had})/#sigma(2J)', '#Delta M(t_{had})/#sigma(3J)', 'A.U.'),
    #('Rel_Delta_mWhad_vs_Rel_Delta_mtophad_6PJ', '#Delta M(W_{had})/#sigma(2J)', '#Delta M(t_{had})/#sigma(3J)', 'A.U.'),
    ('nusolver_chi2', '\chi^2', 'A.U', ''),
    ('wjets_bcMVA_p11', 'cMVA^{11}', 'A.U', ''),
    ('wjets_wcMVA_p11', 'cMVA^{11}', 'A.U.', ''),
    #('wjets_bqgt', 'QG Tag', 'A.U.', ''),
    #('wjets_wqgt', 'QG Tag', 'A.U.', ''),
    ('lb_ratio', 'lb ratio', 'A.U.', ''),
    ('w2b_ratio', 'w2b ratio', 'A.U.', '')
]


jet_only_shapes = []
systematics = ['nosys']#, 'jes_up', 'jes_down', 'jer_up', 'met_up', 'met_down']

right = ['semilep_visible_right']
wrong = ['semilep_wrong', 'semilep_right_thad', 'semilep_right_tlep']
###


output_name_base = 'prob'

def merge_views(inview, subnames):
   subviews = [views.SubdirectoryView(inview, i) for i in subnames]
   return views.SumView(*subviews)

jobid = os.environ['jobid']
plotter = BasePlotter(
   'plots/%s/permutations' % jobid,
    defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}},
)
plotter.reset()
    

#################### scale ttbar xsec for combining 3 different files
ttJets_fnames = ['ttJetsM0', 'ttJetsM700', 'ttJetsM1000']
#ttJets_fnames = ['ttJets']
lumis = []
for fname in ttJets_fnames:
    lumis.append(float(open('inputs/%s/%s.lumi' % (jobid, fname), 'read').readline()))
max_lumi = max(lumis)
lumis[:] = [x/max_lumi for x in lumis]

#set_trace()
ttJets_files = [ #scales found in lumi files normalized by largest value
    (glob.glob('results/%s/permProbComputer/%s.root' % (jobid, ttJets_fnames[0])), lumis[0]),
    (glob.glob('results/%s/permProbComputer/%s.root' % (jobid, ttJets_fnames[1])), lumis[1]),
    (glob.glob('results/%s/permProbComputer/%s.root' % (jobid, ttJets_fnames[2])), lumis[2])
]

# rebin x and y axex
    # xaxis
#maxmjet_bins = np.linspace(0.,100.,26)
#xarray_bins = np.linspace(120., 200., 5)
maxmjet_bins = np.linspace(0.,100.,11)
xarray_bins = np.array([120., 160., 200.])

xbins = np.concatenate((maxmjet_bins,xarray_bins), axis=0)
#print xbins
    # yaxis
#mbpjet_bins = np.linspace(0.,200.,21)
#yarray_bins = np.linspace(300.,2000.,18)
mbpjet_bins = np.linspace(0.,200.,21)
yarray_bins = np.array([400., 800., 1200., 1600., 2000.])

ybins = np.concatenate((mbpjet_bins,yarray_bins), axis=0)
#print ybins

bins_mbpjet = np.concatenate((np.linspace(0., 250., 26), np.array([400., 600.])), axis=0)
#######


#write output file
outname = 'inputs/%s/INPUT/%s_%s.root' % (jobid,output_name_base, args.out)
with io.root_open(outname, 'w') as out:
    out.mkdir('nosys').cd()

    ### 3J event hists
    for i in range(len(right_3J)):
        if i == 0:
            evt_type = 'merged'
        else:
            evt_type = 'lost'
        for shape, xtit, ytit, ztit in shapes_3J: ### make dists for merged and lost events
            hright = 0 # initialize hright hist
            hwrong = 0 # initialize hwrong hist
            for fnames, scales in ttJets_files:
            ### combine same right/wrong shape hists from each ttjets file and normalize based on lumi scale
                for fname in fnames:
                    tfile = io.root_open(fname)
                    right_path = '/'.join([right_3J[i], shape])
                    wrong_path = '/'.join([wrong_3J[i], shape])
                    hright += scales*asrootpy(tfile.Get(right_path)).Clone()
                    hwrong += scales*asrootpy(tfile.Get(wrong_path)).Clone()

            #set_trace()    
            out.cd('nosys')

            #format hright hists
            if shape == 'mbpjet' and evt_type == 'lost':
                #print shape, bins_mbpjet
                hright = RebinView.rebin(hright, bins_mbpjet)
            if ztit:
               hright = RebinView.newRebin2D(hright, xbins, ybins)
               hright.drawstyle = 'colz'
               ROOT.gStyle.SetPalette(56)
               #plotter.canvas.SetRightMargin(0.18)
            hright.name = '3J_%s_%s' % (shape, evt_type+'_right')
            hright.xaxis.title = xtit            
            hright.yaxis.title = ytit
            hright.yaxis.SetTitleOffset(1.5)
            hright.zaxis.title = ztit            
            plotter.pad.SetLeftMargin(0.11)
            plotter.pad.SetRightMargin(0.15)
            hright.zaxis.SetLabelSize(0.8*hright.zaxis.GetLabelSize())
            hright.Draw()
            plotter.save(hright.name)
            hright.Write()
    
            #format hwrong hists
            #if shape == 'mbpjet' and evt_type == 'lost':
            #    #print shape, bins_mbpjet
            #    hwrong = RebinView.rebin(hwrong, bins_mbpjet)
            if ztit:
               hwrong = RebinView.newRebin2D(hwrong, xbins, ybins)
               hwrong.drawstyle = 'colz'
               ROOT.gStyle.SetPalette(56)
               #plotter.canvas.SetRightMargin(0.18)
            hwrong.name = '3J_%s_%s' % (shape, evt_type+'_wrong')
            hwrong.xaxis.title = xtit            
            hwrong.yaxis.title = ytit
            hwrong.yaxis.SetTitleOffset(1.5)
            hwrong.zaxis.title = ztit            
            plotter.pad.SetLeftMargin(0.11)
            plotter.pad.SetRightMargin(0.15)
            hwrong.zaxis.SetLabelSize(0.8*hwrong.zaxis.GetLabelSize())
            hwrong.Draw()
            plotter.save(hwrong.name)
            hwrong.Write()

    ### resolved jet hists
    for i in range(len(right)):
        #if i == 0:
        #    evt_type = 'merged'
        #else:
        #    evt_type = 'lost'
        for shift in systematics:
            for shape, xtit, ytit, ztit in shapes: ### make dists for merged and lost events
                hright = 0 # initialize hright hist
                hwrong = 0 # initialize hwrong hist
                for fnames, scales in ttJets_files:
                ### combine same right/wrong shape hists from each ttjets file and normalize based on lumi scale
                    for fname in fnames:
                        tfile = io.root_open(fname)
                        right_path = '/'.join([right[i], shift, shape])
                        wrong_path = '/'.join([wrong[i], shift, shape])
                        hright += scales*asrootpy(tfile.Get(right_path)).Clone()
                        hwrong += scales*asrootpy(tfile.Get(wrong_path)).Clone()
    
                #set_trace()    
                out.cd('nosys')
    
                #format hright hists
                if ztit:
                   hright.drawstyle = 'colz'
                   ROOT.gStyle.SetPalette(56)
                   #plotter.canvas.SetRightMargin(0.18)
                hright.name = '%s_%s' % (shape, 'right')
                hright.xaxis.title = xtit            
                hright.yaxis.title = ytit
                hright.yaxis.SetTitleOffset(1.5)
                hright.zaxis.title = ztit            
                plotter.pad.SetLeftMargin(0.11)
                plotter.pad.SetRightMargin(0.15)
                hright.zaxis.SetLabelSize(0.8*hright.zaxis.GetLabelSize())
                hright.Draw()
                plotter.save(hright.name)
                hright.Write()
        
                #format hwrong hists
                #if shape == 'mbpjet' and evt_type == 'lost':
                #    #print shape, bins_mbpjet
                #    hwrong = RebinView.rebin(hwrong, bins_mbpjet)
                if ztit:
                   hwrong.drawstyle = 'colz'
                   ROOT.gStyle.SetPalette(56)
                   #plotter.canvas.SetRightMargin(0.18)
                hwrong.name = '%s_%s' % (shape, 'wrong')
                hwrong.xaxis.title = xtit            
                hwrong.yaxis.title = ytit
                hwrong.yaxis.SetTitleOffset(1.5)
                hwrong.zaxis.title = ztit            
                plotter.pad.SetLeftMargin(0.11)
                plotter.pad.SetRightMargin(0.15)
                hwrong.zaxis.SetLabelSize(0.8*hwrong.zaxis.GetLabelSize())
                hwrong.Draw()
                plotter.save(hwrong.name)
                hwrong.Write()


    log.info("%s written." % outname)

