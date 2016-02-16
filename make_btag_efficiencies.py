import rootpy.plotting.views as views
import rootpy.io as io
import os, glob
import rootpy
import URAnalysis.PlotTools.views as urviews
from pdb import set_trace
rootpy.log["/"].setLevel(rootpy.log.INFO)
log = rootpy.log["/make_permutation_distros.py"]
log.setLevel(rootpy.log.ERROR)

jet_types = [('bjet', 'bottom'), ('cjet', 'charm'), ('ljet', 'light')]

jobid = os.environ['jobid']
input_file = 'results/%s/btag_topology_effs/ttJets.root' % jobid
tfile = io.root_open(input_file)

pt_bins  = [0., 30, 60, 100, 150, 200, 1000]
eta_bins = [-2.4, -1.4, -0.8, 0.0, 0.8, 1.4, 2.4]
hview = urviews.RebinView(tfile, [pt_bins, eta_bins])

alljet_cut_types = set([i.name.split('_')[1] for i in tfile.nosys.alljets.keys()])
wjet_cut_types   = set([i.name.split('_')[1] for i in tfile.nosys.Wjets.keys()])
wjet_cut_types.discard('NONE')

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

with io.root_open('inputs/%s/INPUT/btag_general_efficiencies.root' % jobid, 'recreate') as outfile:   
   for jtype, dname in jet_types:
      jdir = outfile.mkdir(dname)
      jdir.cd()
      for cut_type in alljet_cut_types:
         log.info("computing efficiency for %s, %s jets" % (cut_type, jtype))
         h_all  = hview.Get('nosys/alljets/btag_%s_%s_all'  % (cut_type, jtype))
         h_pass = hview.Get('nosys/alljets/btag_%s_%s_pass' % (cut_type, jtype))
         ## h_all.Rebin2D()
         ## h_pass.Rebin2D()
         eff = make_efficiency(h_pass, h_all)
         jdir.WriteTObject(eff, '%s_eff' % cut_type)


with io.root_open('inputs/%s/INPUT/btag_wjets_efficiencies.root' % jobid, 'recreate') as outfile:   
   for jtype, dname in jet_types:
      jdir = outfile.mkdir(dname)
      jdir.cd()
      for cut_type in wjet_cut_types:
         log.info("computing efficiency for %s, %s jets" % (cut_type, jtype))
         h_all  = hview.Get('nosys/Wjets/WjetTag_%s_%s_all'  % (cut_type, jtype))
         h_pass = hview.Get('nosys/Wjets/WjetTag_%s_%s_pass' % (cut_type, jtype))
         ## h_all.Rebin2D()
         ## h_pass.Rebin2D()
         eff = make_efficiency(h_pass, h_all)
         jdir.WriteTObject(eff, '%s_eff' % cut_type)
