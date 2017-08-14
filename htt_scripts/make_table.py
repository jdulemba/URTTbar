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

parser = argparse.ArgumentParser(description='Create plots using files from ttbar_reco_3J.')

jobid = os.environ['jobid']

parser.add_argument('analysis', help='Choose type of analysis (Test or Full).')
parser.add_argument('analyzer', help='Choose type of analyzer (ttbar_reco_3J or htt_simple).')
#parser.add_argument('sample', help='Choose a file (ttJetsM0, ttJetsM700, ttJetsM1000, AtoTT_M...).')
#parser.add_argument('plot', help='Choose type of plots to generate (Gen_Plots, Reco_Plots, Resolution_Plots).')
args = parser.parse_args()

if args.analyzer == "ttbar":
    analyzer = "ttbar_reco_3J"

elif args.analyzer == "htt":
    analyzer = "htt_simple"

rows = []
rows.append(("Sample", "Expected Events", "Expected Merged", "Expected Other", "Merged Fraction", "Other Fraction"))

def print_table(lines, separate_head=True):
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

    #Print the actual data
    for i,line in enumerate(lines):
        print(print_string.format(*line))
        if (i == 0 and separate_head):
            print("-"*(sum(widths)+3*(len(widths)-1)))



AtoTT_files = []

if args.analysis == "Test":
    for f in os.listdir('..'):
        if f.startswith('AtoTT'):
            AtoTT_files.append(f.replace(".root", ""))

#    for fname in AtoTT_files:
#        myfile = root_open('../%s.ttbar_reco_3J.test.root' % fname, 'read')
#        lumifile = open('../inputs/%s/%s.lumi' % (jobid, fname), 'read')
#
#        lumi = lumifile.readline() #luminosity from lumi file
#        scale = 50000/float(lumi)
##        print 'scale: ', scale
#        
#        hist = asrootpy(myfile.Get('Expected_Plots/Expected_Event_Categories')).Clone()
#
#        print fname
#        for i in range(10):
#            print hist.Integral(i,i+1)
#        tot = hist.Integral(0,1)
#        Exp_tot = tot*scale
#        merge = hist.Integral(2,3)
#        Exp_merge = merge*scale
#        other = hist.Integral(4,5)
#        Exp_other = other*scale
#
#        rows.append((fname, format(Exp_tot, '.1f'), format(Exp_merge, '.1f'), format(Exp_other, '.1f')))

    fname = 'AtoTT_M750_5pc_Peak'
    myfile = root_open('../%s.%s.test.root' % (fname, analyzer), 'read')
    lumifile = open('../inputs/%s/%s.lumi' % (jobid, fname), 'read')
    
    lumi = lumifile.readline() #luminosity from lumi file
    scale = 50000/float(lumi)
#     print 'scale: ', scale
    
    hist = asrootpy(myfile.Get('Expected_Plots/Expected_Event_Categories')).Clone()
    
    print fname
    print hist.Integral(0,1)
    print hist.Integral(2,3)
    print hist.Integral(4,5)
#    for i in range(10):
#        print hist.Integral(i,i+1)


elif args.analysis == "Full":
    for f in os.listdir('../'):
    #for f in os.listdir('../results/%s/%s' % (jobid, analyzer)):
        if f.startswith('AtoTT') or f.startswith('ttJetsM'):
            AtoTT_files.append(f.replace(".root", ""))
#        if f.startswith('ttJetsM'):
#            TTJets_files.append(f.replace(".root", ""))

    for fname in AtoTT_files:
        myfile = root_open('../%s.root' % fname, 'read')
        #myfile = root_open('../results/%s/%s/%s.root' % (jobid, analyzer, fname), 'read')
        lumifile = open('../inputs/%s/ttJetsM1000.lumi' % jobid, 'read')
        #lumifile = open('../inputs/%s/%s.lumi' % (jobid, fname), 'read')

        lumi = lumifile.readline() #luminosity from lumi file
        scale = 50000/float(lumi)
#        print 'scale: ', scale
        
        hist_3J = asrootpy(myfile.Get('Expected_Plots/Expected_Event_Categories_3J')).Clone()

        tot_3J = hist_3J.Integral(0,1)
        Exp_tot_3J = tot_3J*scale
        merge_3J = hist_3J.Integral(2,3)
        Exp_merge_3J = merge_3J*scale
        other_3J = hist_3J.Integral(4,5)
        Exp_other_3J = other_3J*scale

        print "merged fraction: ", Exp_merge_3J/Exp_other_3J

        rows.append((fname, format(Exp_tot_3J, '.1f'), format(Exp_merge_3J, '.1f'), format(Exp_other_3J, '.1f'), format(Exp_merge_3J/Exp_tot_3J, '.3f'), format(Exp_other_3J/Exp_tot_3J, '.3f')))

print_table(rows)

