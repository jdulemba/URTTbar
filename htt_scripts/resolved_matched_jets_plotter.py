'''
Resolved_Matched_Jets Analyzer Plotter macro
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os, glob, sys
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
import rootpy.plotting as plotting
from rootpy.plotting import views, Graph, Hist
from argparse import ArgumentParser
import URAnalysis.Utilities.roottools as urtools
import ROOT
from URAnalysis.Utilities.tables import latex_table
import argparse
#import matplotlib.pyplot as plt
from collections import OrderedDict


analyzer = 'resolved_matched_jets'

parser = argparse.ArgumentParser(description='Create plots using files from jet_perm_disc.')

jobid = jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000).')
parser.add_argument('plot', help='Choose type of plots to generate:\n Discriminants(Discs, Matched_Best_Perm)\n Tagger_Vars(TP_Tagger, MP_Tagger)\n Everything')
#parser.add_argument('plot', help='Choose type of plots to generate (Had_Comp, Delta_Plots, Discs, NS_Chi, Combine_Discs, TP_Tagger, MP_Tagger, Best_Perm, Best_Perm_Combined).')
args = parser.parse_args()


##### check analysis type
if not (args.analysis == 'Test' or args.analysis == 'Full'):
    print "You chose %s as your analysis type.\n You must choose Full or Test!" % args.analysis
    sys.exit()


##### check sample type
results_files = []
if args.analysis == "Test":
    for f in os.listdir('../'):
        if not '%s.test.root' % analyzer in f:
            continue
        results_files.append(f.replace(".root", ""))

    if not '%s.%s.test' % (args.sample, analyzer) in results_files:
        print "You chose %s as your sample file.\nYou must choose from the %s files in URTTbar!" % (args.sample, analyzer)
        sys.exit()

if args.analysis == "Full":
    for f in os.listdir('../results/%s/%s' % (jobid, analyzer)):
        if not '.root' in f:
            continue
        results_files.append(f.replace(".root", ""))

    if not args.sample in results_files:
        print "You chose %s as your sample file.\nYou must choose from the files in results/%s/%sJ!" % (args.sample, jobid, analyzer)
        sys.exit()


##### check plot type
if not (args.plot == "Everything" or args.plot == "Discs" or args.plot == "Gen_Plots":
    print "You chose %s as your plot type.\nYou must choose from the help list!"
    sys.exit()


print( 'Analysis: %s\nSample: %s\nPlot: %s' % (args.analysis, args.sample, args.plot) )


if args.analysis == "Test":
    myfile = root_open('../%s.%s.test.root' % (args.sample, analyzer), 'read')
    normfile = views.NormalizeView(root_open('../%s.%s.test.root' % (args.sample, analyzer), 'read'))#normalized file

elif args.analysis == "Full":
    myfile = root_open('../results/%s/%s/%s.root' % (jobid, analyzer, args.sample), 'read')
    normfile = views.NormalizeView(root_open('../results/%s/%s/%s.root' % (jobid, analyzer, args.sample), 'read'))

plotter = BasePlotter(
    'plots/%s/%s/%s/%s' % (analyzer, jobid, args.analysis, args.sample),
    defaults = {'show_title': False, 'save' : {'png' : True, 'pdf' : False}}
    #defaults = {'show_title': True, 'save' : {'png' : True, 'pdf' : False}, 'watermark': ['(13 TeV, 25ns)', False]}
)

def plot_dir(plot):
    print '\ncp -r /uscms/home/jdulemba/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/htt_scripts/plots/%s/%s/%s/%s/%s .\n' % (analyzer, jobid, args.analysis, args.sample, plot)

##### Global Var. Definitions ####
defcol = 'black' # default color
defyax = 'A.U.' # yaxis title
if args.sample == "ttJetsM0":
    mass_min = 0
    mass_max = 2000
if args.sample == "ttJetsM700":
    mass_min = 700
    mass_max = 1000
if args.sample == "ttJetsM1000":
    mass_min = 1000
    mass_max = 2000




######################################################################################
def Gen_Plots():
    plots = 'Gen_Plots'

    disc_type_dict = {'Correct' : '4PJ_M3j_vs_M2j_correct', 'Wrong' : '4PJ_M3j_vs_M2j_wrong'}
    for d_type in disc_type_dict:
        disc_dir = '%s_4PJ_Disc_Plots' % disc_type_dict[d_type]
        



######################################################################################

