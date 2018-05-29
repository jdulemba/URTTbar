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
output_name_base = 'prob'

def merge_views(inview, subnames):
   subviews = [views.SubdirectoryView(inview, i) for i in subnames]
   return views.SumView(*subviews)

jobid = os.environ['jobid']
#input_files = glob.glob('ttJetsM0.permProbComputer.test.root')
input_files = glob.glob('results/%s/permProbComputer/ttJetsM0.root' % jobid)
#
##input_files = glob.glob('results/%s/permProbComputer/ttJets.root' % jobid)
plotter = BasePlotter(
   'plots/%s/permutations' % jobid,
)
plotter.reset()


#################### scale ttbar xsec for combining 3 different files
ttJets_fnames = ['ttJetsM0', 'ttJetsM700', 'ttJetsM1000']
lumis = []
for fname in ttJets_fnames:
    lumis.append(float(open('inputs/%s/%s.lumi' % (jobid, fname), 'read').readline()))
max_lumi = max(lumis)
lumis[:] = [x/max_lumi for x in lumis]

##scales found in lumi files normalized by largest value
ttJets_files = [ #scales found in lumi files normalized by largest value
    (glob.glob('results/%s/permProbComputer/%s.root' % (jobid, ttJets_fnames[0])), lumis[0]),
    (glob.glob('results/%s/permProbComputer/%s.root' % (jobid, ttJets_fnames[1])), lumis[1]),
    (glob.glob('results/%s/permProbComputer/%s.root' % (jobid, ttJets_fnames[2])), lumis[2])
]


log.info('found %d input files' % len(input_files))


for fname in input_files:
    sample = extract_sample(fname)
    tfile = io.root_open(fname)
    test_dir = tfile.Get(right[0])
    right_view = merge_views(tfile, right)
    wrong_view = merge_views(tfile, wrong)
   
    #write output file
    outname = 'inputs/%s/INPUT/%s_%s.root' % (jobid,output_name_base, args.out)
    with io.root_open(outname, 'w') as out:
        for shift in systematics:
            if not hasattr(test_dir, shift):
               log.warning('I could not find %s in %s, skipping systematic' % (shift, sample))
               continue
            out.mkdir(shift).cd()
        #out.mkdir('nosys').cd()
        for shape, xtit, ytit, ztit in shapes:
            path = '/'.join([shift, shape])
            hright = right_view.Get(path)
            #hright = right_view.Get(shape)
            
            hright.name = '%s_%s' % (shape, 'right')
            hright.xaxis.title = xtit            
            hright.yaxis.title = ytit
            hright.yaxis.SetTitleOffset(1.5)
            hright.zaxis.title = ztit            
            if ztit:
               hright.drawstyle = 'colz'
               ROOT.gStyle.SetPalette(56)
               #plotter.canvas.SetRightMargin(0.18)
            #plotter.plot(hright)
            plotter.pad.SetLeftMargin(0.11)
            plotter.pad.SetRightMargin(0.15)
            hright.zaxis.SetLabelSize(0.8*hright.zaxis.GetLabelSize())
            hright.Draw()
            plotter.save(shape)
            hright.Write()


            hwrong = wrong_view.Get(path)
##          hwrong = wrong_view.Get(shape)
#
            hwrong.name = '%s_%s' % (shape, 'wrong')
            hwrong.xaxis.title = xtit            
            hwrong.yaxis.title = ytit
            hwrong.yaxis.SetTitleOffset(1.5)
            hwrong.zaxis.title = ztit            
            if ztit:
               hwrong.drawstyle = 'colz'
               ROOT.gStyle.SetPalette(56)
               #plotter.canvas.SetRightMargin(0.18)
            #plotter.plot(hwrong)
            plotter.pad.SetLeftMargin(0.11)
            plotter.pad.SetRightMargin(0.15)
            hwrong.zaxis.SetLabelSize(0.8*hwrong.zaxis.GetLabelSize())
            hwrong.Draw()
            plotter.save(shape)
            hwrong.Write()
#         
#       for shape in jet_only_shapes:
#          hright = tfile.Get('/'.join([right[0], shift, shape]))
#          hright.name = '%s_%s' % (shape, 'right')
#          hright.Write()
#
#          hwrong = tfile.Get('/'.join([wrong[0], shift, shape]))
#          hwrong.name = '%s_%s' % (shape, 'wrong')
#          hwrong.Write()
    log.info("%s written." % outname)

