import rootpy.plotting.views as views
import rootpy.io as io
import os, glob
import rootpy
from URAnalysis.PlotTools.data_views import extract_sample
from pdb import set_trace
rootpy.log["/"].setLevel(rootpy.log.INFO)
log = rootpy.log["/make_permutation_distros.py"]

jet_types = [('bjet', 'bottom'), ('cjet', 'charm'), ('ljet', 'light')]
cut_types = ['loose', 'tight']

jobid = os.environ['jobid']
input_file = 'results/%s/permProbComputer/ttJets.root' % jobid
tfile = io.root_open(input_file)

def make_efficiency(hpass, hall):
   eff = hpass.Clone()
   eff.Reset()
   for bin_eff, bin_pass, bin_all in zip(eff, hpass, hall):
      bin_eff.value = bin_pass.value/bin_all.value if bin_all.value else 0
   return eff

with io.root_open('inputs/%s/INPUT/btag_efficiencies.root' % jobid, 'recreate') as outfile:
   for jtype, dname in jet_types:
      jdir = outfile.mkdir(dname)
      jdir.cd()
      for cut_type in cut_types:
         h_all  = tfile.Get('semilep_visible_right/nosys/btag_%s_%s_all'  % (cut_type, jtype)) 
         h_pass = tfile.Get('semilep_visible_right/nosys/btag_%s_%s_pass' % (cut_type, jtype))
         ## h_all.Rebin2D()
         ## h_pass.Rebin2D()
         eff = make_efficiency(h_pass, h_all)
         jdir.WriteTObject(eff, '%s_eff' % cut_type)
