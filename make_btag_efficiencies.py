import rootpy.plotting.views as views
import rootpy.io as io
import os, glob
import rootpy
import URAnalysis.PlotTools.views as urviews
from pdb import set_trace
rootpy.log["/"].setLevel(rootpy.log.INFO)
log = rootpy.log["/make_permutation_distros.py"]

jet_types = [('bjet', 'bottom'), ('cjet', 'charm'), ('ljet', 'light')]
cut_types = ['loose', 'tight']

jobid = os.environ['jobid']
input_file = 'results/%s/permProbComputer/ttJets.root' % jobid
tfile = io.root_open(input_file)

pt_bins  = [0.]+[10.*i+50 for i in range(25)]+[20.*i+300. for i in range(5)]+[30.*i+400. for i in range(10)]+[700, 750, 800, 900, 1000]
eta_bins = [-2.4]+[0.2*i-2. for i in range(20)]+[2., 2.4]
hview = urviews.RebinView(tfile, [pt_bins, eta_bins])

def make_efficiency(hpass, hall):
   eff = hpass.Clone()
   eff.Reset()
   for bin_eff, bin_pass, bin_all in zip(eff, hpass, hall):
      if not bin_eff.overflow and bin_pass.value < 20:
         log.warning('bin (%.0f, %.2f) has %i entries' % (bin_eff.x.center, bin_eff.y.center, bin_pass.value))
      bin_eff.value = bin_pass.value/bin_all.value if bin_all.value else 0
      if bin_eff.value == 0 and not bin_eff.overflow:
         log.error('bin (%.0f, %.2f) has 0 efficiency' % (bin_eff.x.center, bin_eff.y.center))
   return eff

with io.root_open('inputs/%s/INPUT/btag_efficiencies.root' % jobid, 'recreate') as outfile:   
   for jtype, dname in jet_types:
      jdir = outfile.mkdir(dname)
      jdir.cd()
      for cut_type in cut_types:
         log.info("computing efficiency for %s, %s jets" % (cut_type, jtype))
         h_all  = hview.Get('semilep_visible_right/nosys/btag_%s_%s_all'  % (cut_type, jtype))
         h_pass = hview.Get('semilep_visible_right/nosys/btag_%s_%s_pass' % (cut_type, jtype))
         ## h_all.Rebin2D()
         ## h_pass.Rebin2D()
         eff = make_efficiency(h_pass, h_all)
         jdir.WriteTObject(eff, '%s_eff' % cut_type)
