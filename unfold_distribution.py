from pdb import set_trace
import os, glob, sys, logging, rootpy, itertools, math
import rootpy.plotting as plotting
import rootpy.io as io
import ROOT
import math
from URAnalysis.AnalysisTools.unfolding.urunfolding import URUnfolding
rootpy.log["/"].setLevel(rootpy.log.INFO)
log = rootpy.log["/URUnfolding"]
rootpy.log.basic_config_colorized()
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('var', type=str, help='varible to unfold')
parser.add_argument('fit_file', type=str, help='file where to find the fitting info')
parser.add_argument('truth_file', type=str, help='file where to find the truth info (migration matrix and true distribution)')
parser.add_argument('-o', type=str, dest='out', default='result_unfolding.root', help='output file')
parser.add_argument('-d', type=str, dest='dir', default='', help='output directory')

## parser.add_argument('--noplots', dest='noplots', action='store_true',
##                     help='skip plot making')
## parser.add_argument('--noshapes', dest='noshapes', action='store_true',
##                     help='skip shape making')
opts = parser.parse_args()

def linearize(graph, axes="X"):
    'axes X or XY'
    axes = axes.upper()
    npts = graph.get_n()
    retval = plotting.Graph(npts)
    for i in xrange(npts):
        x = graph.get_x()[i]
        y = graph.get_y()[i] 
        if 'X' in axes:
            x = 10**x
        if 'Y' in axes:
            y = 10**y
        retval.set_point(i, x, y)
    #retval.get_xaxis().SetTitle(graph.get_xaxis().GetTitle())
    #retval.get_yaxis().SetTitle(graph.get_yaxis().GetTitle())
    return retval


def make_cov_matrix(full_cov, h_input):
    '''pick from the fit covariance matrix 
    only the values we care about'''
    nbins = h_input.GetNbinsX()
    full_nbins = full_cov.GetNbinsX()
    matrix = plotting.Hist2D(
        nbins, h_input.GetBinLowEdge(1), h_input.GetBinLowEdge(nbins+1),
        nbins, h_input.GetBinLowEdge(1), h_input.GetBinLowEdge(nbins+1)
        )
    bin_map_x = dict((full_cov.xaxis.GetBinLabel(i), i) 
                   for i in range(1, full_nbins+1))
    bin_map_y = dict((full_cov.yaxis.GetBinLabel(i), i) 
                   for i in range(1, full_nbins+1))
    to_save = [h_input.xaxis.GetBinLabel(i) 
               for i in range(1, nbins+1)]
    bin_map_new = dict((j, i+1) for i,j in enumerate(to_save))
    for i, j in itertools.product(to_save, to_save):
        full_i = bin_map_x[i]
        full_j = bin_map_y[j]

        new_i = bin_map_new[i]
        new_j = bin_map_new[j]
        
        cov_ij = full_cov.GetBinContent(full_i,full_j)*\
           h_input.GetBinError(new_i) *\
           h_input.GetBinError(new_j)

        if i == j and cov_ij < 0:
            set_trace()
        matrix.SetBinContent(
            new_i,
            new_j,
            cov_ij
            )
    return matrix

def save(canvas, sub, name):
    folder = os.path.join(opts.dir, sub)
    if not os.path.isdir(folder):
        os.mkdir(folder)
    canvas.SaveAs(
        '%s.png' % os.path.join(folder, name))
    canvas.SaveAs(
        '%s.pdf' % os.path.join(folder, name))

def overlay(reference, target):
    canvas = plotting.Canvas(name='adsf', title='asdf')
    reference.linecolor = 'red'
    reference.Draw('hist')
    
    target.markerstyle = 20
    target.Draw('E1P same')
    
    maxreference = reference.GetBinContent(reference.GetMaximumBin())
    maxtarget = target.GetBinContent(target.GetMaximumBin())
    minreference = reference.GetBinContent(reference.GetMinimumBin())
    mintarget = target.GetBinContent(target.GetMinimumBin())

    maxy = max(maxreference, maxtarget)*1.2
    miny = min(minreference, mintarget)
    miny = 0 if miny>0 else miny*1.2
    reference.GetYaxis().SetRangeUser(miny, maxy)
    return canvas


canvas = plotting.Canvas(name='adsf', title='asdf')
resp_file = io.root_open(opts.truth_file)
data_file = io.root_open(opts.fit_file)
scale = 1.
myunfolding = URUnfolding()
myunfolding.matrix   = getattr(resp_file, opts.var).migration_matrix
#log.warning("Using the MC reco distribution for the unfolding!")
myunfolding.measured = getattr(data_file, opts.var).tt_right
#myunfolding.measured = getattr(resp_file, opts.var).reco_distribution
myunfolding.truth    = getattr(resp_file, opts.var).true_distribution
myunfolding.cov_matrix = make_cov_matrix(
    data_file.correlation_matrix,
    getattr(data_file, opts.var).tt_right
    )
myunfolding.InitUnfolder()
hdata = myunfolding.measured # Duplicate. Remove!

#optimize
best_taus = {}
t_min = 0.00001
t_max = 7
best_l, l_curve, graph_x, graph_y  = myunfolding.DoScanLcurve(100, t_min, t_max)
best_taus['L_curve'] = best_l
l_curve.SetName('lcurve')
l_curve.name = 'lcurve'
graph_x.name = 'l_scan_x'
graph_y.name = 'l_scan_y'
l_tau = math.log10(best_l)
points = [(graph_x.GetX()[i], graph_x.GetY()[i], graph_y.GetY()[i]) 
          for i in xrange(graph_x.GetN())]
best = [(x,y) for i, x, y in points if l_tau == i]
graph_best = plotting.Graph(1)
graph_best.SetPoint(0, *best[0])
graph_best.SetMarkerStyle(29)
graph_best.SetMarkerSize(3)
graph_best.SetMarkerColor(2)
l_curve.Draw('ALP')
graph_best.Draw('P SAME')
info = ROOT.TPaveText(.9,0.9,.999,0.999, "brNDC")
info.SetFillColor(0)
info.SetBorderSize(0)
info.SetMargin(0.)
info.AddText('#tau = %.5f' % best_l)
info.Draw()
save(canvas, 'L_curve', 'L_curve')

modes = ['RhoMax', 'RhoSquareAvg', 'RhoAvg']
for mode in modes:
    best_tau, tau_curve = myunfolding.DoScanTau(100, t_min, t_max, mode)
    best_taus[mode] = best_tau
    tau_curve.SetName(mode)
    tau_curve.SetMarkerStyle(1)
    points = [(tau_curve.GetX()[i], tau_curve.GetY()[i])
              for i in xrange(tau_curve.GetN())]
    best = [(x,y) for x, y in points if math.log10(best_tau) == x]
    graph_best = plotting.Graph(1)
    graph_best.SetPoint(0, *best[0])
    graph_best.SetMarkerStyle(29)
    graph_best.SetMarkerSize(3)
    graph_best.SetMarkerColor(2)
    tau_curve.Draw('ALP')
    graph_best.Draw('P SAME')
    info = ROOT.TPaveText(.9,0.9,.999,0.999, "brNDC")
    info.SetFillColor(0)
    info.SetBorderSize(0)
    info.SetMargin(0.)
    info.AddText('#tau = %.5f' % best_tau)
    info.Draw()
    save(canvas, mode, 'L_curve')

for name, best_tau in best_taus.iteritems():
    log.warning('best tau option for %s: %.3f' % (name, best_tau))
    
to_save = []
for name, best_tau in best_taus.iteritems():
    myunfolding.tau = best_tau

    hdata_unfolded = myunfolding.unfolded
    hdata_unfolded.name = 'hdata_unfolded_%s' % name
    to_save.append(hdata_unfolded)
    hdata_refolded = myunfolding.refolded
    hdata_refolded.name = 'hdata_refolded_%s' % name
    to_save.append(hdata_refolded)
    error_matrix = myunfolding.ematrix_total
    error_matrix.name = 'error_matrix_%s' % name
    to_save.append(error_matrix)

    hcorrelations = myunfolding.rhoI_total
    hcorrelations.name = 'hcorrelations_%s' % name
    to_save.append(error_matrix)
    hbias = myunfolding.bias
    hbias.name = 'bias_%s' % name
    to_save.append(hbias)
    canvas = overlay(myunfolding.truth, hdata_unfolded)
    save(canvas, name, 'unfolding')
    canvas = overlay(myunfolding.measured, hdata_refolded)
    save(canvas, name, 'refolded')

htruth = myunfolding.truth
hmatrix = myunfolding.matrix
hmeasured = myunfolding.measured

numexp = 10
myunfolding.unfoldingparam = 10
#uncertainty10 = myunfolding.StatTest(numexp)

with rootpy.io.root_open(os.path.join(opts.dir, opts.out),'recreate') as outfile:
    to_save.extend([
        myunfolding.truth,     ## 4
        myunfolding.measured,  ## 5
        myunfolding.matrix,    ## 6
        l_curve,               ## 9
        tau_curve,             ## 10
        graph_x,
        graph_y,
        ])

    for i, j in enumerate(to_save):
        log.debug('Saving %s as %s' % (j.name, j.GetName()))
        j.Write()
    getattr(resp_file, opts.var).reco_distribution.Write()
    getattr(resp_file, opts.var).prefit_distribution.Write()
    myunfolding.write_to(outfile, 'urunfolder')
