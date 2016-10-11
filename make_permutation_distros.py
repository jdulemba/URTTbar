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

parser = ArgumentParser()
parser.add_argument('out', help='output file name')
#parser.add_argument('--pdfs', action='store_true', help='make plots for the PDF uncertainties')
args = parser.parse_args()


#configuration
shapes = [
   ('mWhad_vs_mtophad', 'M(W_{had})', 'M(t_{had})', 'A.U.'),
   ('nusolver_chi2', '\chi^2', 'A.U', ''),
	 ('wjets_bcMVA_p11', 'cMVA^{11}', 'A.U', ''),
	 ('wjets_wcMVA_p11', 'cMVA^{11}', 'A.U.', ''),
	 ('wjets_bqgt', 'QG Tag', 'A.U.', ''),
	 ('wjets_wqgt', 'QG Tag', 'A.U.', ''),
	 ('lb_ratio', 'lb ratio', 'A.U.', ''),
	 ('w2b_ratio', 'w2b ratio', 'A.U.', ''),
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
input_files = glob.glob('results/%s/permProbComputer/ttJets*.root' % jobid)
plotter = BasePlotter(
   'plots/%s/permutations' % jobid,
)
plotter.reset()

log.info('found %d input files' % len(input_files))

for fname in input_files:
   sample = extract_sample(fname)
   tfile = io.root_open(fname)
   test_dir = tfile.Get(right[0])
   right_view = merge_views(tfile, right)
   wrong_view = merge_views(tfile, wrong)
   #write output file
   outname = 'inputs/%s/INPUT/%s.root' % (jobid, args.out)
   with io.root_open(outname, 'w') as out:
      for shift in systematics:
         if not hasattr(test_dir, shift):
            log.warning('I could not find %s in %s, skipping systematic' % (shift, sample))
            continue
         out.mkdir(shift).cd()
         for shape, xtit, ytit, ztit in shapes:
            path = '/'.join([shift, shape])
            hright = right_view.Get(path)
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
            hwrong.name = '%s_%s' % (shape, 'wrong')
            hwrong.Write()
         
         for shape in jet_only_shapes:
            hright = tfile.Get('/'.join([right[0], shift, shape]))
            hright.name = '%s_%s' % (shape, 'right')
            hright.Write()

            hwrong = tfile.Get('/'.join([wrong[0], shift, shape]))
            hwrong.name = '%s_%s' % (shape, 'wrong')
            hwrong.Write()
   log.info("%s written." % outname)
