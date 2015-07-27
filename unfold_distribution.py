from pdb import set_trace
import os, glob, sys, logging, rootpy, itertools, math, time
import rootpy.plotting as plotting
import rootpy.io as io
import ROOT
import math
from math import sqrt
from URAnalysis.AnalysisTools.unfolding.urunfolding import URUnfolding
from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
from URAnalysis.Utilities.quad import quad
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
import URAnalysis.Utilities.prettyjson as prettyjson

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
parser.add_argument('--tau', type=float, default=-1., help='Fix tau')
parser.add_argument('--mintoy', type=int, default=0, help='Limit toy analysis (>=)')
parser.add_argument('--maxtoy', type=int, default=-1., help='Limit toy analysis (<)')
parser.add_argument('--no_area_constraint', action='store_true', dest='no_area_constraint', help='Do not use the area constraint in the unfolding')
parser.add_argument('--runHandmade', action='store_true', help='Run the handmade tau scan')
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
    if 'thadpt' in variable:
        return 'p_{T}(t_{had}) [GeV]'
    elif 'tleppt' in variable:
        return 'p_{T}(t_{lep}) [GeV]'
    elif 'thadeta' in variable:
        return '|#eta(t_{had})|'
    elif 'tlepeta' in variable:
        return '|#eta(t_{lep})|'
    else:
        return variable
    return ''

def run_unfolder(itoy = 0, outdir = opts.dir, tau = opts.tau):
    
    plotter = BasePlotter()
    
    #canvas = plotting.Canvas(name='adsf', title='asdf')
    if "toy" in opts.fit_file:
        data_file_basedir = 'toy_' + str(itoy)
        data_file_dir = data_file_basedir + '/' + opts.var
    else:
        data_file_dir = opts.var
    xaxislabel = set_pretty_label(opts.var)
    scale = 1.
    if opts.no_area_constraint:
        area_constraint='None'
    else:
        area_constraint='Area'
    myunfolding = URUnfolding(regmode = opts.reg_mode, constraint = area_constraint)

    ## Migration matrix preprocessing
    ## remove oflow bins
    var_dir = getattr(resp_file, opts.var)
    migration_matrix = var_dir.migration_matrix
    for bin in migration_matrix: 
        if bin.overflow:
            bin.value = 0 
            bin.error = 0
    myunfolding.matrix = migration_matrix
    thruth_unscaled = var_dir.thruth_unscaled
    reco_unscaled = var_dir.reco_unscaled
    project_reco = 'X' if myunfolding.orientation == 'Vertical' else 'Y'
    project_gen = 'Y' if myunfolding.orientation == 'Vertical' else 'X'
    reco_project = rootpy.asrootpy(
        getattr(migration_matrix, 'Projection%s' % project_reco)()
        )
    gen_project = rootpy.asrootpy(
        getattr(migration_matrix, 'Projection%s' % project_gen)()
        )
    eff_correction = ROOT.TGraphAsymmErrors(gen_project, thruth_unscaled)
    purity_correction = ROOT.TGraphAsymmErrors(reco_project, reco_unscaled)

    #flush graphs into histograms (easier to handle)
    eff_hist = gen_project.Clone()
    eff_hist.reset()
    eff_hist.name = 'eff_hist'
    for idx in range(eff_correction.GetN()):
        eff_hist[idx+1].value = eff_correction.GetY()[idx]
        eff_hist[idx+1].error = max(
            eff_correction.GetEYhigh()[idx],
            eff_correction.GetEYlow()[idx]
            )

    purity_hist = reco_project.Clone()
    purity_hist.reset()
    purity_hist.name = 'purity_hist'
    for idx in range(purity_correction.GetN()):
        purity_hist[idx+1].value = purity_correction.GetY()[idx]
        purity_hist[idx+1].error = max(
            purity_correction.GetEYhigh()[idx],
            purity_correction.GetEYlow()[idx]
            )

    #Get measured histogram
    measured = None
    if opts.use_reco_truth:
        log.warning("Using the MC reco distribution for the unfolding!")
        measured = getattr(resp_file, opts.var).reco_distribution
    else:
        measured = getattr(data_file, data_file_dir).tt_right

    measured_no_correction = measured.Clone()
    measured_no_correction.name = 'measured_no_correction'
    measured.name = 'measured'
    measured.multiply(purity_hist)
    myunfolding.measured = measured

    #get gen-level distribution
    gen_distro = getattr(resp_file, opts.var).true_distribution.Clone()
    full_true  = gen_distro.Clone()
    full_true.name = 'complete_true_distro'
    gen_distro.multiply(eff_hist)
    gen_distro.name = 'true_distribution'    
    myunfolding.truth = gen_distro
    
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

    #plot covariance matrix
    canvas = plotting.Canvas()
    BasePlotter.set_canvas_style(canvas)
    input_corr_matrix.SetStats(False)
    input_corr_matrix.Draw('colz')
    canvas.Update()
    canvas.SetLogz(True)
    canvas.SaveAs('%s/correlation_matrix.png' % outdir)

    #optimize
    best_taus = {}
    if tau >= 0:
        best_taus['External'] = tau
    else:
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
            tau_curve.SetName('%s_scan' % mode)
            tau_curve.SetMarkerStyle(1)
            points = [(tau_curve.GetX()[i], tau_curve.GetY()[i])
                      for i in xrange(tau_curve.GetN())]
            best = [points[index_best]] 

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

        if opts.runHandmade:
            #hand-made tau scan
            unc_scan, bias_scan = myunfolding.scan_tau(
                200, 10**-6, 50, os.path.join(outdir, 'Handmade', 'scan_info.root'))

            bias_scan.name = 'Handmade'
            bias_scan.title = 'Avg. Bias - Handmade'
            bias_scan.markerstyle = 20
            bias_scan.markersize = 2
            canvas = plotter.create_and_write_graph_canvas('bias_scan',[0],[1],True,True,[bias_scan],write=False)
            save(canvas, 'Handmade', 'bias_scan', outdir)

            unc_scan.name = 'Handmade'
            unc_scan.title = 'Avg. Unc. - Handmade'
            unc_scan.markerstyle = 20
            unc_scan.markersize = 2
            canvas = plotter.create_and_write_graph_canvas('unc_scan',[0],[1],True,True,[unc_scan],write=False)
            save(canvas, 'Handmade', 'unc_scan', outdir)
        
            bias_points = [(bias_scan.GetX()[i], bias_scan.GetY()[i])
                           for i in xrange(bias_scan.GetN())]
            unc_points = [(unc_scan.GetX()[i], unc_scan.GetY()[i])
                           for i in xrange(unc_scan.GetN())]
            fom_scan = plotting.Graph(unc_scan.GetN())
            for idx, info in enumerate(zip(bias_points, unc_points)):
                binfo, uinfo = info
                tau, bias = binfo
                _, unc = uinfo
                fom_scan.SetPoint(idx, tau, quad(bias, unc))
            fom_scan.name = 'Handmade'
            fom_scan.title = 'Figure of merit - Handmade'
            fom_scan.markerstyle = 20
            fom_scan.markersize = 2
            canvas = plotter.create_and_write_graph_canvas('unc_scan',[0],[1],True,True,[fom_scan],write=False)
            save(canvas, 'Handmade', 'fom_scan', outdir)

    to_save = []
    outfile = rootpy.io.root_open(os.path.join(outdir, opts.out),'recreate')
    for name, best_tau in best_taus.iteritems():
        method_dir = outfile.mkdir(name)
        myunfolding.tau = best_tau

        hdata_unfolded = myunfolding.unfolded
        #apply phase space efficiency corrections
        hdata_unfolded_ps_corrected = hdata_unfolded.Clone()
        hdata_unfolded_ps_corrected.Divide(eff_hist)

        hdata_refolded = myunfolding.refolded
        #apply purity corrections
        hdata_refolded_wpurity = hdata_refolded.Clone()

        error_matrix = myunfolding.ematrix_total

        hcorrelations = myunfolding.rhoI_total
        hbias = myunfolding.bias
        #canvas = overlay(myunfolding.truth, hdata_unfolded)
        leg = LegendDefinition(
            title=name,
            labels=['Truth','Unfolded'],
            position='ne'
            )
        myunfolding.truth.xaxis.title = xaxislabel
        hdata_unfolded.xaxis.title = xaxislabel
        n_neg_bins = 0
        for ibin in range(1,hdata_unfolded.GetNbinsX()+1):
            if hdata_unfolded.GetBinContent(ibin) < 0:
                n_neg_bins = n_neg_bins + 1
        hn_neg_bins = plotting.Hist(2,-1, 1, name = 'nneg_bins', title = 'Negative bins in ' + hdata_unfolded.GetName()+ ';Bin sign; N_{bins}')
        hn_neg_bins.SetBinContent(1,n_neg_bins)
        hn_neg_bins.SetBinContent(2,hdata_unfolded.GetNbinsX()-n_neg_bins)
        canvas = plotter.create_and_write_canvas_single(1,0,1,False,False,hn_neg_bins, write=False)
        save(canvas, name, 'unfolding_bins_sign', outdir)


        sumofpulls = 0
        sumofratios = 0
        for ibin in range(1,myunfolding.truth.GetNbinsX()+1):
            binContent1 = myunfolding.truth.GetBinContent(ibin)
            binContent2 = hdata_unfolded.GetBinContent(ibin)
            binError1 = myunfolding.truth.GetBinError(ibin)
            binError2 = hdata_unfolded.GetBinError(ibin)
            error = sqrt(binError1*binError1 + binError2*binError2)
            if error != 0:
                pull = (binContent2-binContent1)/error
            else:
                pull = 9999
            if binContent1 != 0:
                ratio = binContent2/binContent1
            sumofpulls = sumofpulls + pull
            sumofratios = sumofratios + ratio
        sumofpulls = sumofpulls / myunfolding.truth.GetNbinsX()
        sumofratios = sumofratios / myunfolding.truth.GetNbinsX()
        
        hsum_of_pulls = plotting.Hist(1,0, 1, name = 'sum_of_pulls_' + hdata_unfolded.GetName(), title = 'Sum of pulls wrt truth for ' + hdata_unfolded.GetName()+ ';None; #Sigma(pulls) / N_{bins}')
        hsum_of_pulls.SetBinContent(1, sumofpulls)
        canvas = plotter.create_and_write_canvas_single(1,0,1,False,False,hsum_of_pulls, write=False)
        save(canvas, name, 'unfolding_sum_of_pulls', outdir)
        
        hsum_of_ratios = plotting.Hist(1,0, 1, name = 'sum_of_ratios_' + hdata_unfolded.GetName(), title = 'Sum of ratios wrt truth for ' + hdata_unfolded.GetName()+ ';None; #Sigma(ratios) / N_{bins}')
        hsum_of_ratios.SetBinContent(1, sumofratios)
        canvas = plotter.create_and_write_canvas_single(1,0,1,False,False,hsum_of_ratios, write=False)
        save(canvas, name, 'unfolding_sum_of_ratios', outdir)

        canvas = plotter.create_and_write_canvas_with_comparison('unfolding_pull', [1,0],[0,21], [2,1], leg, False, False, [myunfolding.truth, hdata_unfolded], write=False, comparison='pull')
        save(canvas, name, 'unfolding_pull', outdir)
        canvas = plotter.create_and_write_canvas_with_comparison('unfolding_ratio', [1,0],[0,21], [2,1], leg, False, False, [myunfolding.truth, hdata_unfolded], write=False, comparison='ratio')
        save(canvas, name, 'unfolding_ratio', outdir)

        canvas = plotter.create_and_write_canvas_with_comparison('unfolding_weff_pull', [1,0],[0,21], [2,1], leg, False, False, [full_true, hdata_unfolded_ps_corrected], write=False, comparison='pull')
        save(canvas, name, 'unfolding_weff_pull', outdir)
        canvas = plotter.create_and_write_canvas_with_comparison('unfolding_weff_ratio', [1,0],[0,21], [2,1], leg, False, False, [full_true, hdata_unfolded_ps_corrected], write=False, comparison='ratio')
        save(canvas, name, 'unfolding_weff_ratio', outdir)
    
        nbins = myunfolding.measured.GetNbinsX()
        #for i in range(1, nbins+1):
            #myunfolding.measured.GetXaxis().SetBinLabel(
                #i, 
                #'%.0f' % myunfolding.measured.GetXaxis().GetBinLowEdge(i)
                #)
        input_distro = getattr(resp_file, opts.var).prefit_distribution
        leg = LegendDefinition(title=name,labels=['Reco','Refolded','Input'],position='ne')
        myunfolding.measured.xaxis.title = xaxislabel
        hdata_refolded.xaxis.title = xaxislabel
        myunfolding.measured.drawstyle = 'e1'
        canvas = plotter.create_and_write_canvas_with_comparison('refolded_pull', [1,0,1],[0,21,0], [2,1,4], leg, False, False, [myunfolding.measured, hdata_refolded], write=False, comparison='pull')
        save(canvas, name, 'refolded_pull', outdir)
        canvas = plotter.create_and_write_canvas_with_comparison('refolded_ratio', [1,0,1],[0,21,0], [2,1,4], leg, False, False, [myunfolding.measured, hdata_refolded], write=False, comparison='ratio')
        save(canvas, name, 'refolded_ratio', outdir)

        measured_no_correction.drawstyle = 'e1'
        canvas = plotter.create_and_write_canvas_with_comparison('refolded_wpurity_pull', [1,0,1],[0,21,0], [2,1,4], leg, False, False, [measured_no_correction, hdata_refolded_wpurity, input_distro], write=False, comparison='pull')
        save(canvas, name, 'refolded_wpurity_pull', outdir)
        canvas = plotter.create_and_write_canvas_with_comparison('refolded_wpurity_ratio', [1,0,1],[0,21,0], [2,1,4], leg, False, False, [measured_no_correction, hdata_refolded_wpurity, input_distro], write=False, comparison='ratio')
        save(canvas, name, 'refolded_wpurity_ratio', outdir)

        method_dir.WriteTObject(hdata_unfolded, 'hdata_unfolded')
        method_dir.WriteTObject(hdata_unfolded_ps_corrected, 'hdata_unfolded_ps_corrected')
        method_dir.WriteTObject(hdata_refolded, 'hdata_refolded')
        method_dir.WriteTObject(hdata_refolded_wpurity, 'hdata_refolded_wpurity')
        method_dir.WriteTObject(error_matrix, 'error_matrix')
        method_dir.WriteTObject(hbias, 'bias')
        method_dir.WriteTObject(hn_neg_bins, 'hn_neg_bins')
        method_dir.WriteTObject(hsum_of_pulls, 'hsum_of_pulls')
        method_dir.WriteTObject(hsum_of_ratios, 'hsum_of_ratios')


    htruth = myunfolding.truth
    hmatrix = myunfolding.matrix
    hmeasured = myunfolding.measured

    #with rootpy.io.root_open(os.path.join(outdir, opts.out),'recreate') as outfile:
    outfile.cd()
    to_save.extend([
        measured_no_correction,
        eff_hist,
        purity_hist,
        full_true,
        myunfolding.truth,     ## 4
        myunfolding.measured,  ## 5
        myunfolding.matrix,])  ## 6
    if opts.tau < 0:
        to_save.extend([
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
    json = ROOT.TText(0., 0., prettyjson.dumps(best_taus))
    outfile.WriteTObject(json, 'best_taus')
    myunfolding.write_to(outfile, 'urunfolder')
    outfile.Close()

resp_file = io.root_open(opts.truth_file)
data_file = io.root_open(opts.fit_file)

#Draw migration matrix and its errors (w/ Under/Over-flow)
migration_matrix   = getattr(resp_file, opts.var).migration_matrix
xbins = set(i.x.high for i in migration_matrix) 
xbins.update(
    set(i.x.low for i in migration_matrix)
    )
xbins = list(xbins)
xbins.sort()

ybins = set(i.y.high for i in migration_matrix)
ybins.update(
    set(i.y.low for i in migration_matrix)
    )
ybins = list(ybins)
ybins.sort()

full_matrix = plotting.Hist2D(xbins, ybins)
for bin in migration_matrix:
    ibin = full_matrix.FindFixBin(
        bin.x.center, bin.y.center
        )
    full_matrix[ibin].value = bin.value
    full_matrix[ibin].error = bin.error

canvas = plotting.Canvas(1000, 800)
BasePlotter.set_canvas_style(canvas)
full_matrix.SetStats(False)
full_matrix.Draw('colz')
canvas.Update()
canvas.SetLogz(True)
canvas.SaveAs('%s/migration_matrix.png' % opts.dir)

#make relative uncertainty
for bin in full_matrix:
    if not bin.value: continue
    bin.value = bin.error / bin.value 
full_matrix.Draw('colz')
canvas.Update()
canvas.SaveAs('%s/migration_matrix_unc.png' % opts.dir)

if 'toy' in opts.fit_file:
    toys = [i.name for i in data_file.keys()]
    toys = [(i, int(i.replace('toy_',''))) for i in toys if i.startswith('toy_')]
    toys = [(i, j) for i, j in toys if j >= opts.mintoy]
    if opts.maxtoy > 0:
        if opts.mintoy >= opts.maxtoy:
            raise ValueError(
                'You ran with option --maxtoy=%i, but the option'
                ' --mintoy=%i, which makes no sense!' % (opts.maxtoy, opts.mintoy)
                )
        toys = [(i, j) for i, j in toys if j < opts.maxtoy]
    for data_file_basedir, itoy in toys:
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
else:
    run_unfolder()
    
if 'toy' in opts.fit_file:
    print "unfolding_toy_diagnostics %s %s" % (opts.dir, opts.var)
    unfolding_toy_diagnostics(opts.dir, opts.var)



