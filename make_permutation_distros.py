import rootpy.plotting.views as views
import rootpy.io as io
import os, glob
import rootpy
from URAnalysis.PlotTools.data_views import extract_sample
from pdb import set_trace
rootpy.log["/"].setLevel(rootpy.log.INFO)
log = rootpy.log["/make_permutation_distros.py"]

#configuration
shapes = ['mWhad_vs_mtophad', 'nusolver_chi2']
jet_only_shapes = ['btag_value']
systematics = ['nosys', 'jes_up', 'jes_down', 'jer_up', 'met_up', 'met_down']
right = ['semilep_visible_right']
wrong = ['semilep_wrong', 'semilep_right_thad', 'semilep_right_tlep', 'other_tt_decay']
output_name_base = 'prob'

def merge_views(inview, subnames):
   subviews = [views.SubdirectoryView(inview, i) for i in subnames]
   return views.SumView(*subviews)

jobid = os.environ['jobid']
input_files = glob.glob('results/%s/permProbComputer/ttJets*.root' % jobid)

log.info('found %d input files' % len(input_files))

for fname in input_files:
   sample = extract_sample(fname)
   tfile = io.root_open(fname)
   test_dir = tfile.Get(right[0])
   right_view = merge_views(tfile, right)
   wrong_view = merge_views(tfile, wrong)
   #write output file
   outname = 'inputs/%s/INPUT/%s_%s.root' % (jobid, output_name_base, sample)
   with io.root_open(outname, 'w') as out:
      for shift in systematics:
         if not hasattr(test_dir, shift):
            log.warning('I could not find %s in %s, skipping systematic' % (shift, sample))
            continue
         out.mkdir(shift).cd()
         for shape in shapes:
            path = '/'.join([shift, shape])
            hright = right_view.Get(path)
            hright.name = '%s_%s' % (shape, 'right')
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
