import rootpy.plotting.views as views
import rootpy.io as io
import os, glob
import rootpy
from URAnalysis.PlotTools.data_views import data_views 
from URAnalysis.PlotTools.data_views import extract_sample 
from pdb import set_trace
import ROOT
rootpy.log["/"].setLevel(rootpy.log.INFO)
log = rootpy.log["/make_permutation_distros.py"]
from URAnalysis.PlotTools.BasePlotter import BasePlotter
from URAnalysis.PlotTools.Plotter import Plotter, BasePlotter
from styles import styles
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
	 #('wjets_bqgt', 'QG Tag', 'A.U.', ''),
	 #('wjets_wqgt', 'QG Tag', 'A.U.', ''),
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


class PERMPlotter(Plotter):
    def __init__(self, lumi=None):
        'Inits the plotter and identifies what will be plotted'
        lumi = lumi if lumi > 0 else None
        jobid = os.environ['jobid']
        project = os.environ['URA_PROJECT']

        ## get files, lumis
        ttJetsfiles = glob.glob('results/%s/permProbComputer/*.root' % jobid )
        data_files =  glob.glob('results/%s/ctag_eff/*.root' % jobid )
        data_files = [ fname for fname in data_files if (os.path.basename(fname).startswith('data') ) ]
        files = ttJetsfiles+data_files
        
        lumis = glob.glob('inputs/%s/*.lumi' % jobid)
        data_lumis = [ fname for fname in lumis if (os.path.basename(fname).startswith('data') ) ]
        ttJets_lumis = [ fname for fname in lumis if os.path.basename(fname).startswith('ttJets') and 'up' not in fname and 'down' not in fname ]
        lumis = data_lumis+ttJets_lumis

        outdir = 'plots/%s/permutations' % jobid

        super(PERMPlotter, self).__init__(
            files, lumis, outdir, styles, None, lumi
        )

        self.defaults = {
            'blurb': [13, self.views['data']['intlumi']],
            'save' : {'png' : True, 'pdf' : False}
        }
        self.jobid = jobid
        self.project = project

plotter = PERMPlotter()

ttJets_view = plotter.views['ttJetsAll']['view']
right_view = merge_views(ttJets_view, right)
wrong_view = merge_views(ttJets_view, wrong)
#write output file
outname = 'inputs/%s/INPUT/%s.root' % (plotter.jobid, args.out)
with io.root_open(outname, 'w') as out:
   for shift in systematics:
      #if not hasattr(test_dir, shift):
      #   log.warning('I could not find %s in %s, skipping systematic' % (shift, sample))
      #   continue
      out.mkdir(shift).cd()
      for shape, xtit, ytit, ztit in shapes:
         path = '/'.join([shift, shape])
         #set_trace()
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
