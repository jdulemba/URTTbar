'''
This file makes tables for the expected number of partially merged events for samples.
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os, glob, sys
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

parser = argparse.ArgumentParser(description='Create tables using files from htt_simple or ttbar_reco_3J.')

jobid = os.environ['jobid']

parser.add_argument('analyzer', help='Choose type of analyzer (ttbar or htt).')
parser.add_argument('sample', help='Choose type of sample (AtoTT, HtoTT, ttJets, or All).')
args = parser.parse_args()



##### check analyzer type
if not (args.analyzer == "ttbar" or args.analyzer == "htt"):
    print "You chose %s as your analyzer.\nYou must choose from the help menu!" % args.analyzer
    sys.exit()


##### check sample type
if not (args.sample == "AtoTT" or args.sample == "HtoTT" or args.sample == "ttJets" or args.sample == "All"):
    print "You chose %s as your sample type.\nYou must choose from the help menu!" % args.sample
    sys.exit()
    


if args.analyzer == "ttbar":
    analyzer = "ttbar_reco_3J"
    print "This is a check for the ttJetsM files against the performance in htt_simple."

elif args.analyzer == "htt":
    analyzer = "htt_simple"

rows = []
rows.append(("Sample", "Expected Events 3J (4+J)", "Expected Partial Merge", "Expected Other", "Expected 3J/4+J Gain", "Partial Merge Fraction", "Other Fraction"))

def print_table(lines, filename, separate_head=True):
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

    with open(filename, 'w') as f:
        #Print the actual data
        for i,line in enumerate(lines):
            print >> f, print_string.format(*line)
            if (i == 0 and separate_head):
                print >> f, "-"*(sum(widths)+3*(len(widths)-1))


def Expected_Events(fname):
        myfile = root_open('../../results/%s/%s/%s.root' % (jobid, analyzer, fname), 'read')
        #lumifile = open('../inputs/%s/ttJetsM1000.lumi' % jobid, 'read')
        lumifile = open('../../inputs/%s/%s.lumi' % (jobid, fname), 'read')

        lumi = lumifile.readline() #luminosity from lumi file
        scale = 50000/float(lumi)

        # 3J
        hist_3J = asrootpy(myfile.Get('Expected_Plots/Expected_Event_Categories_3J')).Clone()
        tot_3J = hist_3J.Integral(0,1)
        Exp_tot_3J = tot_3J*scale
        merge_3J = hist_3J.Integral(2,3)
        Exp_merge_3J = merge_3J*scale
        other_3J = hist_3J.Integral(4,5)
        Exp_other_3J = other_3J*scale

        # 4PJ    
        hist_4PJ = asrootpy(myfile.Get('Expected_Plots/Expected_Event_Categories_4PJ')).Clone()
        tot_4PJ = hist_4PJ.Integral(0,1)
        Exp_tot_4PJ = tot_4PJ*scale
        merge_4PJ = hist_4PJ.Integral(2,3)
        Exp_merge_4PJ = merge_4PJ*scale
        other_4PJ = hist_4PJ.Integral(4,5)
        Exp_other_4PJ = other_4PJ*scale
    
        rows.append((fname, format(Exp_tot_3J, '.1f')+' ('+format(Exp_tot_4PJ, '.1f')+')',\
            format(Exp_merge_3J, '.1f')+' ('+format(Exp_merge_4PJ, '.1f')+')', \
            format(Exp_other_3J, '.1f')+' ('+format(Exp_other_4PJ, '.1f')+')',\
            format(Exp_tot_3J/Exp_tot_4PJ, '.3f'),\
            format(Exp_merge_3J/Exp_tot_3J, '.3f')+' ('+format(Exp_merge_4PJ/Exp_tot_4PJ, '.3f')+')',\
            format(Exp_other_3J/Exp_tot_3J, '.3f')+' ('+format(Exp_other_4PJ/Exp_tot_4PJ, '.3f')+')'))



#### initialize lists for different file
all_files = []
AtoTT_files = []
HtoTT_files = []
ttJets_files = []

for f in os.listdir('../../results/%s/%s' % (jobid, analyzer)):
    if f.startswith('AtoTT') or f.startswith('HtoTT') or f.startswith('ttJetsM'):
        all_files.append(f.replace(".root", ""))

    if f.startswith('AtoTT'):
        AtoTT_files.append(f.replace(".root",""))

    if f.startswith('HtoTT'):
        HtoTT_files.append(f.replace(".root",""))

    if f.startswith('ttJetsM'):
        ttJets_files.append(f.replace(".root", ""))


for fname in all_files:
    if 'Int' in fname:
        continue


########## all files

    if args.sample == "All":
        Expected_Events(fname)

########## AtoTT files

    if args.sample == "AtoTT":
        if 'AtoTT' in fname:
            Expected_Events(fname)
    
########## HtoTT files

    if args.sample == "HtoTT":
        if 'HtoTT' in fname:
            Expected_Events(fname)

########## ttJetsM files

    if args.sample == "ttJets":
        if 'ttJetsM' in fname:
            Expected_Events(fname)

print_table(rows, filename='../../plots/%s/%s/%s_Expected_Event_Categories.raw_txt' % (jobid, args.analyzer, args.sample))

print "\n%s_Expected_Event_Categories.raw_txt created.\n" % args.sample

print "The tables can be found in the directory:\n \
/uscms/home/jdulemba/nobackup/CMSSW_7_4_7/src/Analyses/URTTbar/plots/%s/%s" % (jobid, args.analyzer)
