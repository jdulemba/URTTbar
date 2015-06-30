from pdb import set_trace
import os, glob, sys, logging, rootpy, itertools, math, time
import rootpy.plotting as plotting
import rootpy.io as io
import ROOT
import math
from URAnalysis.AnalysisTools.unfolding.urunfolding import URUnfolding
from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
from unfolding_toy_diagnostics import unfolding_toy_diagnostics
rootpy.log["/"].setLevel(rootpy.log.ERROR)
log = rootpy.log["/URUnfolding"]
log.setLevel(rootpy.log.INFO)
rootpy.log.basic_config_colorized()
#ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(False)
ROOT.gStyle.SetObjectStat(False)
ROOT.gROOT.SetBatch()
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('var', type=str, help='varible to unfold')
parser.add_argument('fit_file', type=str, help='file where to find the fitting info')
parser.add_argument('truth_file', type=str, help='file where to find the truth info (migration matrix and true distribution)')
parser.add_argument('-o', type=str, dest='out', default='result_unfolding.root', help='output file')
parser.add_argument('-d', type=str, dest='dir', default='', help='output directory')
parser.add_argument('--cov_matrix', type=str, dest='cov_matrix', default='full', help='Covariance matrix to use: full (diagonal+off-diagonal), diag (diagonal only), none (let TUnfold build a diagonal one).')
parser.add_argument('--use_reco_truth', action='store_true', dest='use_reco_truth', help='Use the reco from migration matrix')
parser.add_argument('--reg_mode', type=str, dest='reg_mode', default='Curvature', help='Regularization mode to use: None, Size, Derivative, Curvature (default), Mixed.')
#parser.add_argument('--tau_range', type=str, dest='tau_range', default='(0.0000001,7)', help='Tau range to scan')
parser.add_argument('--tau_range', type=str, dest='tau_range', default='(0,20)', help='Tau range to scan')
#parser.add_argument('--bias_id', type=str, dest='bias_id', default='', help='bias applied to the true distribution')

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
        if new_i != new_j and opts.cov_matrix == 'diag':
            continue
        matrix.SetBinContent(
            new_i,
            new_j,
            cov_ij
            )
    return matrix
 
def make_corr_matrix(full_cov, h_input):
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
        
        corr_ij = full_cov.GetBinContent(full_i,full_j)

        if i == j and corr_ij < 0:
            set_trace()
        if new_i != new_j and opts.cov_matrix == 'diag':
            continue
        matrix.SetBinContent(
            new_i,
            new_j,
            corr_ij
            )
    return matrix

def save(canvas, sub, name, outdir):
    folder = os.path.join(outdir, sub)
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

def set_pretty_label(variable):
    if 'ptthad' in variable:
        return 'p_{T}(t_{had}) [GeV]'
    elif 'pttlep' in variable:
        return 'p_{T}(t_{lep}) [GeV]'
    elif 'etathad' in variable:
        return '|#eta(t_{had})|'
    elif 'etatlep' in variable:
        return '|#eta(t_{lep})|'
    else:
        return variable
    return ''

def run_unfolder(itoy = 0, outdir = opts.dir):
    
    plotter = BasePlotter()
    
    #canvas = plotting.Canvas(name='adsf', title='asdf')
    if "toy" in opts.fit_file:
        data_file_basedir = 'toy_' + str(itoy)
        data_file_dir = data_file_basedir + '/' + opts.var
    else:
        data_file_dir = opts.var
    xaxislabel = set_pretty_label(opts.var)
    scale = 1.
    myunfolding = URUnfolding(regmode = opts.reg_mode)
    myunfolding.matrix   = getattr(resp_file, opts.var).migration_matrix
    if opts.use_reco_truth:
        log.warning("Using the MC reco distribution for the unfolding!")
        myunfolding.measured = getattr(resp_file, opts.var).reco_distribution
    else:
        myunfolding.measured = getattr(data_file, data_file_dir).tt_right
    #if opts.bias_id == '':
        #myunfolding.truth    = getattr(resp_file, opts.var).true_distribution
    #else:
        #myunfolding.truth    = getattr(resp_file, opts.var + opts.bias_id).true_distribution
    myunfolding.truth    = getattr(resp_file, opts.var).true_distribution
    if opts.cov_matrix != 'none':
        if 'toy' in opts.fit_file:
            input_cov_matrix = make_cov_matrix(
                getattr(data_file, data_file_basedir).correlation_matrix,
                getattr(data_file, data_file_dir).tt_right
                )
            input_corr_matrix = make_corr_matrix(
                getattr(data_file, data_file_basedir).correlation_matrix,
                getattr(data_file, data_file_dir).tt_right
                )
        else:
            input_cov_matrix = make_cov_matrix(
                data_file.correlation_matrix,
                getattr(data_file, data_file_dir).tt_right
                )
            input_corr_matrix = make_corr_matrix(
                data_file.correlation_matrix,
                getattr(data_file, data_file_dir).tt_right
                )
        input_cov_matrix.name = 'input_cov_matrix'
        input_corr_matrix.name = 'input_corr_matrix'
        myunfolding.cov_matrix = input_cov_matrix
    myunfolding.InitUnfolder()
    hdata = myunfolding.measured # Duplicate. Remove!

    #optimize
    best_taus = {}
    t_min, t_max = eval(opts.tau_range)
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
    plotter.set_graph_style(graph_best,29,2)
    graph_best.SetPoint(0, *best[0])
    graph_best.SetMarkerSize(3)
    #l_curve.Draw('ALP')
    #graph_best.Draw('P SAME')
    canvas = plotter.create_and_write_graph_canvas('c'+l_curve.GetName(),[0],[1],False,False,[l_curve],write=False)
    graph_best.Draw('P SAME')
    
    info = ROOT.TPaveText(0.65,1-canvas.GetTopMargin(),1-canvas.GetRightMargin(),0.999, "brNDC")
    info.UseCurrentStyle()
    info.SetTextAlign(32)
    info.SetFillColor(0)
    info.SetBorderSize(0)
    info.SetMargin(0.)
    info.SetTextSize(0.037)
    info.AddText('#tau = %.5f' % best_l)
    info.Draw()
    canvas.Update()
    save(canvas, 'L_curve', 'L_curve', outdir)

    modes = ['RhoMax', 'RhoSquareAvg', 'RhoAvg']
    for mode in modes:
        best_tau, tau_curve, index_best = myunfolding.DoScanTau(100, t_min, t_max, mode)
        best_taus[mode] = best_tau
        tau_curve.SetName(mode)
        tau_curve.SetMarkerStyle(1)
        points = [(tau_curve.GetX()[i], tau_curve.GetY()[i])
                  for i in xrange(tau_curve.GetN())]
        best = [points[index_best]] 
        #best = [(x,y) for x, y in points if math.log10(best_tau) == x]
        #log.info('best = %s and best_from_index = %s' % (best,best_from_index))
        #if  best_from_index != best[0]:
            #log.error("Pair found by DoScanTau is different from pair found by this code!")
            #os.abort()
        #if itoy == 16 and mode == 'RhoSquareAvg':
            #set_trace()
        graph_best = plotting.Graph(1)
        graph_best.SetPoint(0, *best[0])
        graph_best.SetMarkerStyle(29)
        graph_best.SetMarkerSize(3)
        graph_best.SetMarkerColor(2)
        
        canvas = plotter.create_and_write_graph_canvas('c'+tau_curve.GetName(),[0],[1],False,False,[tau_curve],write=False)
    
        graph_best.Draw('P SAME')
        info = ROOT.TPaveText(0.65,1-canvas.GetTopMargin(),1-canvas.GetRightMargin(),0.999, "brNDC")
        info.UseCurrentStyle()
        info.SetFillColor(0)
        info.SetBorderSize(0)
        info.SetMargin(0.)
        info.AddText('#tau = %.5f' % best_tau)
        info.Draw()
        canvas.Update()
        save(canvas, mode, 'Tau_curve', outdir)

    #force running without regularization
    best_taus['NoReg'] = 0
    for name, best_tau in best_taus.iteritems():
        log.info('best tau option for %s: %.3f' % (name, best_tau))
    
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
        #canvas = overlay(myunfolding.truth, hdata_unfolded)
        leg = LegendDefinition(title=name,labels=['Truth','Unfolded'],position='ne')
        myunfolding.truth.xaxis.title = xaxislabel
        hdata_unfolded.xaxis.title = xaxislabel
        canvas = plotter.create_and_write_canvas_with_comparison('unfolding_pull', [1,0],[0,21], [2,1], leg, False, False, [myunfolding.truth, hdata_unfolded], write=False, comparison='pull')
        save(canvas, name, 'unfolding_pull', outdir)
        canvas = plotter.create_and_write_canvas_with_comparison('unfolding_ratio', [1,0],[0,21], [2,1], leg, False, False, [myunfolding.truth, hdata_unfolded], write=False, comparison='ratio')
        save(canvas, name, 'unfolding_ratio', outdir)
    
        nbins = myunfolding.measured.GetNbinsX()
        #for i in range(1, nbins+1):
            #myunfolding.measured.GetXaxis().SetBinLabel(
                #i, 
                #'%.0f' % myunfolding.measured.GetXaxis().GetBinLowEdge(i)
                #)
        leg = LegendDefinition(title=name,labels=['Reco','Refolded'],position='ne')
        myunfolding.measured.xaxis.title = xaxislabel
        hdata_refolded.xaxis.title = xaxislabel
        canvas = plotter.create_and_write_canvas_with_comparison('refolded_pull', [1,0],[0,21], [2,1], leg, False, False, [myunfolding.measured, hdata_refolded], write=False, comparison='pull')
        save(canvas, name, 'refolded_pull', outdir)
        canvas = plotter.create_and_write_canvas_with_comparison('refolded_ratio', [1,0],[0,21], [2,1], leg, False, False, [myunfolding.measured, hdata_refolded], write=False, comparison='ratio')
        save(canvas, name, 'refolded_ratio', outdir)

    htruth = myunfolding.truth
    hmatrix = myunfolding.matrix
    hmeasured = myunfolding.measured

    numexp = 10
    myunfolding.unfoldingparam = 10
    #uncertainty10 = myunfolding.StatTest(numexp)

    with rootpy.io.root_open(os.path.join(outdir, opts.out),'recreate') as outfile:
        to_save.extend([
            myunfolding.truth,     ## 4
            myunfolding.measured,  ## 5
            myunfolding.matrix,    ## 6
            l_curve,               ## 9
            tau_curve,             ## 10
            graph_x,
            graph_y
            ])

        if opts.cov_matrix != 'none':
           to_save.extend([input_cov_matrix])
           to_save.extend([input_corr_matrix])

        for i, j in enumerate(to_save):
            log.debug('Saving %s as %s' % (j.name, j.GetName()))
            j.Write()
        getattr(resp_file, opts.var).reco_distribution.Write()
        getattr(resp_file, opts.var).prefit_distribution.Write()
        myunfolding.write_to(outfile, 'urunfolder')


resp_file = io.root_open(opts.truth_file)
data_file = io.root_open(opts.fit_file)

if 'toy' in opts.fit_file:
    itoy = 0
    while True:
        data_file_basedir = 'toy_' + str(itoy)
        try:
            getattr(data_file, data_file_basedir)
        #except rootpy.io.DoesNotExist as err:
        except AttributeError as err:
            log.warning('Directory %s not found in the rootfile! Breaking the loop.' % data_file_basedir)
            break
        log.info('Found directory %s in the rootfile' % data_file_basedir)

        try:
            os.mkdir(os.path.join(opts.dir,data_file_basedir))
        except OSError as e:
            if e.errno == 17:
                log.warning('Directory %s already exists.' % data_file_basedir)
            else:
                raise e

        outdir = os.path.join(opts.dir,data_file_basedir)
        
        run_unfolder(itoy, outdir)
                
        itoy = itoy + 1
        #if itoy > 10:
            #break
else:
    run_unfolder()
    
if 'toy' in opts.fit_file:
    unfolding_toy_diagnostics(opts.dir, opts.var)



