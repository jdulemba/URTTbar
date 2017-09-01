'''
ttbar_reco_3J Analyzer Plotter macro
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os, glob
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
from rootpy.plotting import views, Graph, Hist
from argparse import ArgumentParser
import URAnalysis.Utilities.roottools as urtools
import ROOT
from URAnalysis.Utilities.tables import latex_table
import argparse
import matplotlib.pyplot as plt
import sys

parser = argparse.ArgumentParser(description='Create plots using files from ttbar_reco_3J or htt_simple.')

jobid = os.environ['jobid']

parser.add_argument('analyzer', help='Choose type of analyzer (ttbar_reco_3J or htt_simple).')
parser.add_argument('sample', help='Choose sample (ttJetsM700 or ttJetsM1000).')
args = parser.parse_args()

if args.analyzer == "ttbar":
    analyzer = "ttbar_reco_3J"

elif args.analyzer == "htt":
    analyzer = "htt_simple"

if not (args.sample == "ttJetsM700" or args.sample == "ttJetsM1000"):
    print "You chose %s as your sample." % args.sample
    print "Must choose ttJetsM700 or ttJetsM1000!"
    sys.exit()

if args.sample == "ttJetsM700":
    bin_range = ["700-800","800-900","900-1000"]

if args.sample == "ttJetsM1000":
    bin_range = ["1000-1200","1200-1500","> 1500"]

rows = []
rows.append(("ttbar Mass Range [GeV]", "Expected Events 3J (4+J)", "Expected Partial Merge", "Expected Other", "Expected 3J/4+J Gain", "Partial Merge Fraction", "Other Fraction"))

def print_table(lines, fname, separate_head=True):
    """Prints a formatted table given a 2 dimensional array"""
    #Count the column width
    widths = []
    for line in lines:
        for i,size in enumerate([len(x) for x in line]):
            while i >= len(widths):
                widths.append(0)
            if size > widths[i]:
                widths[i] = size

    #Generate the format string to pad the columns
    print_string = ""
    for i,width in enumerate(widths):
        print_string += "{" + str(i) + ":" + str(width) + "} | "
    if (len(print_string) == 0):
        return
    print_string = print_string[:-3]

    with open(fname, 'w') as f:
        #Print the actual data
        for i,line in enumerate(lines):
            print >> f, print_string.format(*line)
            if (i == 0 and separate_head):
                print >> f, "-"*(sum(widths)+3*(len(widths)-1))




#TTJets_files = []
#
#for f in os.listdir('../results/%s/%s' % (jobid, analyzer)):
#    if f.endswith('M700.root') or f.endswith('M1000.root'):
#        TTJets_files.append(f.replace(".root", ""))

#for fname in TTJets_files:
myfile = root_open('../results/%s/%s/%s.root' % (jobid, analyzer, args.sample), 'read')
lumifile = open('../inputs/%s/%s.lumi' % (jobid, args.sample), 'read')

lumi = lumifile.readline() #luminosity from lumi file
scale = 50000/float(lumi)


tot_vals = [1,5,9]

hist_3J = asrootpy(myfile.Get('Expected_Plots/Expected_Event_Categories_3J_massbins')).Clone()
hist_4PJ = asrootpy(myfile.Get('Expected_Plots/Expected_Event_Categories_4PJ_massbins')).Clone()

for i in range(len(tot_vals)):
    # 3J
    tot_3J = hist_3J.Integral(tot_vals[i], tot_vals[i])
    Exp_tot_3J = tot_3J*scale
    merge_3J = hist_3J.Integral(tot_vals[i]+1, tot_vals[i]+1)
    Exp_merge_3J = merge_3J*scale
    other_3J = hist_3J.Integral(tot_vals[i]+2, tot_vals[i]+2)
    Exp_other_3J = other_3J*scale
    
    tot_4PJ = hist_4PJ.Integral(tot_vals[i], tot_vals[i])
    Exp_tot_4PJ = tot_4PJ*scale
    merge_4PJ = hist_4PJ.Integral(tot_vals[i]+1, tot_vals[i]+1)
    Exp_merge_4PJ = merge_4PJ*scale
    other_4PJ = hist_4PJ.Integral(tot_vals[i]+2, tot_vals[i]+2)
    Exp_other_4PJ = other_4PJ*scale
    
    rows.append((bin_range[i], format(Exp_tot_3J, '.1f')+' ('+format(Exp_tot_4PJ, '.1f')+')',\
        format(Exp_merge_3J, '.1f')+' ('+format(Exp_merge_4PJ, '.1f')+')', \
        format(Exp_other_3J, '.1f')+' ('+format(Exp_other_4PJ, '.1f')+')',\
        format(Exp_tot_3J/Exp_tot_4PJ, '.3f'),\
        format(Exp_merge_3J/Exp_tot_3J, '.3f')+' ('+format(Exp_merge_4PJ/Exp_tot_4PJ, '.3f')+')',\
        format(Exp_other_3J/Exp_tot_3J, '.3f')+' ('+format(Exp_other_4PJ/Exp_tot_4PJ, '.3f')+')'))

print_table(rows, fname='../plots/%s/%s/%s_Expected_Event_Categories.raw_txt' % (jobid, args.analyzer, args.sample))
