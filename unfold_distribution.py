from pdb import set_trace
import os, glob, sys, logging, rootpy, itertools
from styles import styles
import rootpy.plotting as plotting
import rootpy.plotting.views as views
from array import array
import rootpy.io as io
from fnmatch import fnmatch
import ROOT
import URAnalysis.Utilities.prettyjson as prettyjson
from URAnalysis.AnalysisTools.unfolding.urunfolding import URUnfolding
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('var', type=str, help='varible to unfold')
parser.add_argument('fit_file', type=str, help='file where to find the fitting info')
parser.add_argument('truth_file', type=str, help='file where to find the truth info (migration matrix and true distribution)')
parser.add_argument('-o', type=str, dest='out', default='result_unfolding.root', help='output file')
## parser.add_argument('--noplots', dest='noplots', action='store_true',
##                     help='skip plot making')
## parser.add_argument('--noshapes', dest='noshapes', action='store_true',
##                     help='skip shape making')
opts = parser.parse_args()

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

resp_file = io.root_open(opts.truth_file)
data_file = io.root_open(opts.fit_file)
scale = 1.
myunfolding = URUnfolding()
myunfolding.matrix   = getattr(resp_file, opts.var).migration_matrix
myunfolding.measured = getattr(data_file, opts.var).tt_right
myunfolding.truth    = getattr(resp_file, opts.var).true_distribution
myunfolding.cov_matrix = make_cov_matrix(
    data_file.correlation_matrix,
    myunfolding.measured
    )
myunfolding.InitUnfolder()
hdata = myunfolding.measured # Duplicate. Remove!

#optimize
best_l, l_curve   = myunfolding.DoScanLcurve(100)
l_curve.name = 'lcurve'
best_tau, tau_curve = myunfolding.DoScanTau(100)
tau_curve.SetName('tau_curve')
logging.info('best_tau_option for L: %.2f for scan: %.2f' % (best_l, best_tau))
myunfolding.tau = best_l

hdata_unfolded = myunfolding.unfolded
hdata_unfolded.name = 'hdata_unfolded'
hdata_refolded = myunfolding.refolded
hdata_refolded.name = 'hdata_refolded'
error_matrix = myunfolding.ematrix_total
error_matrix.name = "error_matrix"

hcorrelations = myunfolding.rhoI_total
hcorrelations.name = "hcorrelations"
htruth = myunfolding.truth
hmatrix = myunfolding.matrix
hmeasured = myunfolding.measured
hbias = myunfolding.bias
hbias.name = 'bias'

hgenerated = rootpy.asrootpy(hmatrix.ProjectionY())
hgenerated.name = 'hgenerated'
hreconstructed = rootpy.asrootpy(hmatrix.ProjectionX())
hreconstructed.name = 'hreconstructed'

#hdifference = hdata_unfolded - hgenerated
#hdifference.SetName('hdifference')

#hdifference_refolded = hdata_refolded - hreconstructed
#hdifference_refolded.SetName('hdifference_refolded')

numexp = 10
myunfolding.unfoldingparam = 10
#uncertainty10 = myunfolding.StatTest(numexp)

with rootpy.io.root_open(opts.out,'recreate') as outfile:
    to_write = [
        hdata_unfolded,        ## 0  
        hdata_refolded,        ## 1
        error_matrix,          ## 2
        hcorrelations,         ## 3
        myunfolding.truth,     ## 4
        myunfolding.measured,  ## 5
        myunfolding.matrix,    ## 6
        hgenerated,            ## 7
        hreconstructed,        ## 8
        l_curve,               ## 9
        tau_curve,             ## 10
        hbias,                 ## 11
        ]
    for i, j in enumerate(to_write):
        j.Write()
    myunfolding.write_to(outfile, 'urunfolder')
