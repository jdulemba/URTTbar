import glob
import os
import rootpy.io as io
import rootpy.plotting as plotting
from URAnalysis.Utilities.roottools import Envelope
from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
from pdb import set_trace
import rootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)

from URAnalysis.Utilities.quad import quad

plotter = BasePlotter()
unf_dir = 'plots/2015Jul03/ttxsec/ptthad_fullcov_recotoy_notscaled_Curvature_15Jul21_Staggered'
spath =    os.path.join(
      unf_dir,
      'toy_*',
      'Handmade',
      'scan_info.root'
      )
print spath
tfiles  = glob.glob(
   spath
   )
tfiles = [io.root_open(i) for i in tfiles]
scan_points = [(i.name, int(i.name.split('_')[1])) for i in tfiles[0].keys() if i.name.startswith('pt_')]
scan_points.sort(key=lambda x: x[1])
truth = tfiles[0].truth

biases = plotting.Graph(len(scan_points), title='Bias Scan')
uncs = plotting.Graph(len(scan_points), title='Uncertainty Scan')

for spoint, idx in scan_points:
   tau = None
   envelope = Envelope()
   avg_uncs = []
   for tfile in tfiles:
      unfolded = tfile.Get('%s/unfolded' % spoint)
      envelope += unfolded
      if tau is None:
         tau = tfile.Get('%s/tau' % spoint).value
      rel_unc = [bin.error/abs(bin.value)
                 for bin in unfolded
                 if not bin.overflow]
      avg_uncs.append(sum(rel_unc)/len(rel_unc))
   
   bias = []
   for bu, bt in zip(envelope.median, truth):
      if bu.overflow: continue
      bias.append(
         abs(bu.value - bt.value)/bt.value
         )

   biases.SetPoint(idx, tau, sum(bias)/len(bias))
   uncs.SetPoint(idx, tau, sum(avg_uncs)/len(avg_uncs))


canvas = plotter.create_and_write_graph_canvas('unc_scan',[0],[1],True,True,[uncs],write=False)
uncs.yaxis.SetMoreLogLabels(True)
canvas.SaveAs('%s/unc_scan.png' % unf_dir)

canvas = plotter.create_and_write_graph_canvas('bias_scan',[0],[1],True,True,[biases],write=False)
biases.yaxis.SetMoreLogLabels(True)
canvas.SaveAs('%s/bias_scan.png' % unf_dir)

bias_points = [(biases.GetX()[i], biases.GetY()[i])
               for i in xrange(biases.GetN())]
unc_points = [(uncs.GetX()[i], uncs.GetY()[i])
               for i in xrange(uncs.GetN())]
fom_scan = plotting.Graph(uncs.GetN())
for idx, info in enumerate(zip(bias_points, unc_points)):
    binfo, uinfo = info
    tau, bias = binfo
    _, unc = uinfo
    fom_scan.SetPoint(idx, tau, quad(bias, unc))
fom_scan.name = 'Handmade'
fom_scan.title = 'Figure of merit - Handmade'
canvas = plotter.create_and_write_graph_canvas('fom_scan',[0],[1],True,True,[fom_scan],write=False)
fom_scan.yaxis.SetMoreLogLabels(True)
canvas.SaveAs('%s/fom_scan.png' % unf_dir)

