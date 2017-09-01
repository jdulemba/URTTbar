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
### Joseph added
right_permdisc_shapes = [
   ('3J_mbpjet_vs_maxmjet_correct', 'Max M_{jet} [GeV]', 'M_{b+jet} [GeV]', 'A.U.')
   #('4J_mbpjet_vs_maxmjet_correct', 'Max M_{jet}', 'M_{b+jet}', 'A.U.')
]
right_nschi_shapes = [
   ('3J_nusolver_chi2_correct', '\chi^2', 'A.U', '')
]   

wrong_permdisc_shapes = [
  ('3J_mbpjet_vs_maxmjet_wrong', 'Max M_{jet} [GeV]', 'M_{b+jet} [GeV]', 'A.U.')
  #('4J_mbpjet_vs_maxmjet_wrong', 'Max M_{jet}', 'M_{b+jet}', 'A.U.') 
]
wrong_nschi_shapes = [
   ('3J_nusolver_chi2_wrong', '\chi^2', 'A.U', '')
]

####
#shapes = [
#    ('mWhad_vs_mtophad', 'M(W_{had})', 'M(t_{had})', 'A.U.'),
#    ('mWhad_vs_mtophad_4J', 'M(W_{had})', 'M(t_{had})', 'A.U.'),
#    ('mWhad_vs_mtophad_5J', 'M(W_{had})', 'M(t_{had})', 'A.U.'),
#    ('mWhad_vs_mtophad_6PJ', 'M(W_{had})', 'M(t_{had})', 'A.U.'),
#    ('Rel_Delta_mWhad_vs_Rel_Delta_mtophad', '#Delta M(W_{had})/#sigma(2J)', '#Delta M(t_{had})/#sigma(3J)', 'A.U.'),
#    ('Rel_Delta_mWhad_vs_Rel_Delta_mtophad_4J', '#Delta M(W_{had})/#sigma(2J)', '#Delta M(t_{had})/#sigma(3J)', 'A.U.'),
#    ('Rel_Delta_mWhad_vs_Rel_Delta_mtophad_5J', '#Delta M(W_{had})/#sigma(2J)', '#Delta M(t_{had})/#sigma(3J)', 'A.U.'),
#    ('Rel_Delta_mWhad_vs_Rel_Delta_mtophad_6PJ', '#Delta M(W_{had})/#sigma(2J)', '#Delta M(t_{had})/#sigma(3J)', 'A.U.'),
#    ('nusolver_chi2', '\chi^2', 'A.U', ''),
#    ('wjets_bcMVA_p11', 'cMVA^{11}', 'A.U', ''),
#    ('wjets_wcMVA_p11', 'cMVA^{11}', 'A.U.', ''),
#    #('wjets_bqgt', 'QG Tag', 'A.U.', ''),
#    #('wjets_wqgt', 'QG Tag', 'A.U.', ''),
#    ('lb_ratio', 'lb ratio', 'A.U.', ''),
#    ('w2b_ratio', 'w2b ratio', 'A.U.', '')
#]



jet_only_shapes = []
systematics = ['nosys']#, 'jes_up', 'jes_down', 'jer_up', 'met_up', 'met_down']

### Joseph added
right = ['Correct_Disc_Plots']
wrong = ['Wrong_Disc_Plots']
###

#right = ['semilep_visible_right']
#wrong = ['semilep_wrong', 'semilep_right_thad', 'semilep_right_tlep']
output_name_base = 'prob'

def merge_views(inview, subnames):
   subviews = [views.SubdirectoryView(inview, i) for i in subnames]
   return views.SumView(*subviews)

jobid = os.environ['jobid']
##input_files = glob.glob('ttJetsM0.permProbComputer.test.root')
#input_files = glob.glob('results/%s/permProbComputer/ttJets.root' % jobid)
#
##input_files = glob.glob('results/%s/permProbComputer/ttJets.root' % jobid)
plotter = BasePlotter(
   'plots/%s/permutations' % jobid,
)
plotter.reset()
    

#log.info('found %d input files' % len(input_files))


########### Joseph added
# scale ttbar xsec for combining 3 different files

ttJets_files = [ #scales found in lumi files normalized by largest value
    (glob.glob('results/%s/jet_perm_disc/ttJetsM0.root' % jobid), 180591.247475/527103.981164),
    (glob.glob('results/%s/jet_perm_disc/ttJetsM700.root' % jobid), 373640.713257/527103.981164),
    (glob.glob('results/%s/jet_perm_disc/ttJetsM1000.root' % jobid), 527103.981164/527103.981164)
    ]

    ## perm disc arrays
scaled_3J_permdisc_right = []
scaled_3J_permdisc_wrong = []
#scaled_4J_permdisc_right = []
#scaled_4J_permdisc_wrong = []

    ## ns chi arrays
scaled_3J_nschi_right = []
scaled_3J_nschi_wrong = []
#scaled_4J_nschi_right = []
#scaled_4J_nschi_wrong = []

for fname, scale in ttJets_files: # filename, normalized scale
    for name in fname:
        tfile = io.root_open(name)

    ## perm disc files
#        Right_3J = asrootpy(tfile.Get(right[0]+'/3J_mbpjet_vs_maxmjet_correct')).Clone()
        Scale_3J_permdisc_Right = scale*asrootpy(tfile.Get(right[0]+'/3J_mbpjet_vs_maxmjet_correct')).Clone()
        scaled_3J_permdisc_right.append(Scale_3J_permdisc_Right)

#        Wrong_3J = asrootpy(tfile.Get(wrong[0]+'/3J_mbpjet_vs_maxmjet_wrong')).Clone()
        Scale_3J_permdisc_Wrong = scale*asrootpy(tfile.Get(wrong[0]+'/3J_mbpjet_vs_maxmjet_wrong')).Clone()
        scaled_3J_permdisc_wrong.append(Scale_3J_permdisc_Wrong)

##        Right_4J = asrootpy(tfile.Get(right[0]+'/4J_mbpjet_vs_maxmjet_correct')).Clone()
#        Scale_4J_permdisc_Right = scale*asrootpy(tfile.Get(right[0]+'/4J_mbpjet_vs_maxmjet_correct')).Clone()
#        scaled_4J_permdisc_right.append(Scale_4J_permdisc_Right)
#
##        Wrong_4J = asrootpy(tfile.Get(wrong[0]+'/4J_mbpjet_vs_maxmjet_wrong')).Clone()
#        Scale_4J_permdisc_Wrong = scale*asrootpy(tfile.Get(wrong[0]+'/4J_mbpjet_vs_maxmjet_wrong')).Clone()
#        scaled_4J_permdisc_wrong.append(Scale_4J_permdisc_Wrong)

    ## ns chi files
        nschi_Right_3J = asrootpy(tfile.Get(right[0]+'/3J_nusolver_chi2_correct')).Clone()
        Scale_3J_nschi_Right = scale*asrootpy(tfile.Get(right[0]+'/3J_nusolver_chi2_correct')).Clone()
        scaled_3J_nschi_right.append(Scale_3J_nschi_Right)
#        print '3J correct: ', nschi_Right_3J.Integral()

        nschi_Wrong_3J = asrootpy(tfile.Get(wrong[0]+'/3J_nusolver_chi2_wrong')).Clone()
        Scale_3J_nschi_Wrong = scale*asrootpy(tfile.Get(wrong[0]+'/3J_nusolver_chi2_wrong')).Clone()
        scaled_3J_nschi_wrong.append(Scale_3J_nschi_Wrong)
#        print '3J wrong: ', nschi_Wrong_3J.Integral()

#        Scale_4J_nschi_Right = scale*asrootpy(tfile.Get(right[0]+'/4J_mbpjet_vs_maxmjet_correct')).Clone()
#        scaled_4J_nschi_right.append(Scale_4J_permdisc_Right)
#
#        Scale_4J_nschi_Wrong = scale*asrootpy(tfile.Get(wrong[0]+'/4J_mbpjet_vs_maxmjet_wrong')).Clone()
#        scaled_4J_nschi_wrong.append(Scale_4J_permdisc_Wrong)


## combine perm disc files
Scaled_3J_permdisc_correct = scaled_3J_permdisc_right[0]+scaled_3J_permdisc_right[1]+scaled_3J_permdisc_right[2]
Scaled_3J_permdisc_wrong = scaled_3J_permdisc_wrong[0]+scaled_3J_permdisc_wrong[1]+scaled_3J_permdisc_wrong[2]
Scaled_permdisc_correct = [Scaled_3J_permdisc_correct]
Scaled_permdisc_wrong = [Scaled_3J_permdisc_wrong]
#Scaled_4J_permdisc_correct = scaled_4J_permdisc_right[0]+scaled_4J_permdisc_right[1]+scaled_4J_permdisc_right[2]
#Scaled_4J_permdisc_wrong = scaled_4J_permdisc_wrong[0]+scaled_4J_permdisc_wrong[1]+scaled_4J_permdisc_wrong[2]
##print '3J correct: ', Scaled_3J_permdisc_correct.Integral()
##print '4J corect: ', Scaled_4J_permdisc_correct.Integral()
##print '3J wrong: ', Scaled_3J_permdisc_wrong.Integral()
##print '4J wrong: ', Scaled_4J_permdisc_wrong.Integral()
#Scaled_permdisc_correct = [Scaled_3J_permdisc_correct, Scaled_4J_permdisc_correct]
#Scaled_permdisc_wrong = [Scaled_3J_permdisc_wrong, Scaled_4J_permdisc_wrong]

## combine ns chi files
Scaled_3J_nschi_correct = scaled_3J_nschi_right[0]+scaled_3J_nschi_right[1]+scaled_3J_nschi_right[2]
Scaled_3J_nschi_wrong = scaled_3J_nschi_wrong[0]+scaled_3J_nschi_wrong[1]+scaled_3J_nschi_wrong[2]
#Scaled_4J_nschi_correct = scaled_4J_nschi_right[0]+scaled_4J_nschi_right[1]+scaled_4J_nschi_right[2]
#Scaled_4J_nschi_wrong = scaled_4J_nschi_wrong[0]+scaled_4J_nschi_wrong[1]+scaled_4J_nschi_wrong[2]

#print '3J correct: ', Scaled_3J_permdisc_correct.Integral()
#print '4J corect: ', Scaled_4J_permdisc_correct.Integral()
#print '3J wrong: ', Scaled_3J_permdisc_wrong.Integral()
#print '4J wrong: ', Scaled_4J_permdisc_wrong.Integral()
Scaled_nschi_correct = [Scaled_3J_nschi_correct]#, Scaled_4J_permdisc_correct]
Scaled_nschi_wrong = [Scaled_3J_nschi_wrong]#, Scaled_4J_permdisc_wrong]


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

######


##### Joseph added
#write output file
outname = 'inputs/%s/INPUT/%s_%s.root' % (jobid,output_name_base, args.out)
with io.root_open(outname, 'w') as out:
    out.mkdir('nosys').cd()

## perm disc
    right_permdisc_idx = 0 
    wrong_permdisc_idx = 0 
    for shape, xtit, ytit, ztit in right_permdisc_shapes:
#        hright = Scaled_permdisc_correct[right_permdisc_idx]
        hright = RebinView.newRebin2D(Scaled_permdisc_correct[right_permdisc_idx],xbins,ybins)        
#        print "shape: %s" % shape
#        print "Right int: %f" % hright.Integral()
        
        hright.name = shape
        hright.xaxis.title = xtit            
        hright.yaxis.title = ytit
        hright.yaxis.SetTitleOffset(1.5)
        hright.zaxis.title = ztit            
        if ztit:
            hright.drawstyle = 'colz'
            ROOT.gStyle.SetPalette(56)
        plotter.pad.SetLeftMargin(0.15)
        plotter.pad.SetRightMargin(0.15)
        hright.zaxis.SetLabelSize(0.8*hright.zaxis.GetLabelSize())
        hright.Draw()
        plotter.save(shape)
        hright.Write()
        right_permdisc_idx +=1
    
    for shape, xtit, ytit, ztit in wrong_permdisc_shapes:
#        hwrong = Scaled_permdisc_wrong[wrong_permdisc_idx]
        hwrong = RebinView.newRebin2D(Scaled_permdisc_wrong[wrong_permdisc_idx],xbins,ybins)        
#        print "shape: %s" % shape
#        print "Wrong int: %f" % hwrong.Integral()
        
        hwrong.name = shape
        hwrong.xaxis.title = xtit            
        hwrong.yaxis.title = ytit
        hwrong.yaxis.SetTitleOffset(1.5)
        hwrong.zaxis.title = ztit            
        if ztit:
            hwrong.drawstyle = 'colz'
            ROOT.gStyle.SetPalette(56)
        plotter.pad.SetLeftMargin(0.15)
        plotter.pad.SetRightMargin(0.15)
        hwrong.zaxis.SetLabelSize(0.8*hwrong.zaxis.GetLabelSize())
        hwrong.Draw()
        plotter.save(shape)
        hwrong.Write()
        wrong_permdisc_idx +=1
  

## ns chi 
    right_nschi_idx = 0 
    wrong_nschi_idx = 0 
    for shape, xtit, ytit, ztit in right_nschi_shapes:
        hright = Scaled_nschi_correct[right_nschi_idx]
#        print "shape: %s" % shape
#        print "Right int: %f" % hright.Integral()
        
        hright.name = shape
        hright.xaxis.title = xtit            
        hright.yaxis.title = ytit
        hright.yaxis.SetTitleOffset(1.5)
        hright.zaxis.title = ztit            
        if ztit:
            hright.drawstyle = 'colz'
            ROOT.gStyle.SetPalette(56)
        plotter.pad.SetLeftMargin(0.15)
        plotter.pad.SetRightMargin(0.15)
        hright.zaxis.SetLabelSize(0.8*hright.zaxis.GetLabelSize())
        hright.Draw()
        plotter.save(shape)
        hright.Write()
        right_nschi_idx +=1
    
    for shape, xtit, ytit, ztit in wrong_nschi_shapes:
        hwrong = Scaled_nschi_wrong[wrong_nschi_idx]
#        print "shape: %s" % shape
#        print "Wrong int: %f" % hwrong.Integral()
        
        hwrong.name = shape
        hwrong.xaxis.title = xtit            
        hwrong.yaxis.title = ytit
        hwrong.yaxis.SetTitleOffset(1.5)
        hwrong.zaxis.title = ztit            
        if ztit:
            hwrong.drawstyle = 'colz'
            ROOT.gStyle.SetPalette(56)
        plotter.pad.SetLeftMargin(0.15)
        plotter.pad.SetRightMargin(0.15)
        hwrong.zaxis.SetLabelSize(0.8*hwrong.zaxis.GetLabelSize())
        hwrong.Draw()
        plotter.save(shape)
        hwrong.Write()
        wrong_nschi_idx +=1
  
log.info("%s written." % outname)

###


#for fname in input_files:
#    sample = extract_sample(fname)
#    tfile = io.root_open(fname)
#    test_dir = tfile.Get(right[0])
#    right_view = merge_views(tfile, right)
#    wrong_view = merge_views(tfile, wrong)
#   
#    #write output file
#    outname = 'inputs/%s/INPUT/%s_%s.root' % (jobid,output_name_base, args.out)
#    with io.root_open(outname, 'w') as out:
#        for shift in systematics:
#            if not hasattr(test_dir, shift):
#               log.warning('I could not find %s in %s, skipping systematic' % (shift, sample))
#               continue
#            out.mkdir(shift).cd()
#        #out.mkdir('nosys').cd()
#        for shape, xtit, ytit, ztit in shapes:
#            path = '/'.join([shift, shape])
#            hright = right_view.Get(path)
#            #hright = right_view.Get(shape)
#            
#            hright.name = '%s_%s' % (shape, 'right')
#            hright.xaxis.title = xtit            
#            hright.yaxis.title = ytit
#            hright.yaxis.SetTitleOffset(1.5)
#            hright.zaxis.title = ztit            
#            if ztit:
#               hright.drawstyle = 'colz'
#               ROOT.gStyle.SetPalette(56)
#               #plotter.canvas.SetRightMargin(0.18)
#            #plotter.plot(hright)
#            plotter.pad.SetLeftMargin(0.11)
#            plotter.pad.SetRightMargin(0.15)
#            hright.zaxis.SetLabelSize(0.8*hright.zaxis.GetLabelSize())
#            hright.Draw()
#            plotter.save(shape)
#            hright.Write()
#
#
#            hwrong = wrong_view.Get(path)
###          hwrong = wrong_view.Get(shape)
##
#            hwrong.name = '%s_%s' % (shape, 'wrong')
#            hwrong.xaxis.title = xtit            
#            hwrong.yaxis.title = ytit
#            hwrong.yaxis.SetTitleOffset(1.5)
#            hwrong.zaxis.title = ztit            
#            if ztit:
#               hwrong.drawstyle = 'colz'
#               ROOT.gStyle.SetPalette(56)
#               #plotter.canvas.SetRightMargin(0.18)
#            #plotter.plot(hwrong)
#            plotter.pad.SetLeftMargin(0.11)
#            plotter.pad.SetRightMargin(0.15)
#            hwrong.zaxis.SetLabelSize(0.8*hwrong.zaxis.GetLabelSize())
#            hwrong.Draw()
#            plotter.save(shape)
#            hwrong.Write()
##         
##       for shape in jet_only_shapes:
##          hright = tfile.Get('/'.join([right[0], shift, shape]))
##          hright.name = '%s_%s' % (shape, 'right')
##          hright.Write()
##
##          hwrong = tfile.Get('/'.join([wrong[0], shift, shape]))
##          hwrong.name = '%s_%s' % (shape, 'wrong')
##          hwrong.Write()
#    log.info("%s written." % outname)

