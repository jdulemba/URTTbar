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

jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']

#parser = argparse.ArgumentParser()#description='Create tables using files from htt_simple or ttbar_reco_3J.')
#parser.add_argument('analyzer', help='Choose type of analyzer (ttbar or htt).')
#parser.add_argument('sample', help='Choose type of sample (AtoTT, HtoTT, ttJets, or All).')
#args = parser.parse_args()

rows = []
rows.append(("Run", "ditag", "notag", "subtag", "leadtag", "Total"))

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

def datacard_events(era):
    myfile = root_open('%s/plots/%s/ctageff/%s/mass_discriminant/DeepCsvTight/datacard.root' % (project, jobid, era), 'read')

    tag_names = ['subtag', 'notag', 'leadtag', 'ditag']
    tags ={}
    total = 0
    for tag in tag_names:
        hist = asrootpy(myfile.Get('%s/data_obs' % tag)).Clone()
        tags[tag] = hist.Integral()
        total += hist.Integral()
    tags['Total'] = total
    
    rows.append((era, format(tags['ditag'], '.0f'),\
                format(tags['notag'], '.0f'),\
                format(tags['subtag'], '.0f'),\
                format(tags['leadtag'], '.0f'),\
                format(tags['Total'], '.0f')))


eras = ['All_Runs', 'Run_B', 'Run_CtoE', 'Run_EtoF']

for era in eras:
    datacard_events(era)


print_table(rows, filename='%s/plots/%s/ctageff/Datacard_Yields.raw_txt' % (project, jobid))
print "\nDatacard_Yields.raw_txt and\n"


def ctageff_data_yields(fname):
    myfile = root_open('%s/results/%s/ctag_eff/%s.root' % (project, jobid, fname), 'read')

    tag_names = ['lead_tagged', 'sublead_tagged', 'both_tagged', 'both_untagged']
    tags={}
    total = 0
    for tag in tag_names:
        hist = asrootpy(myfile.Get('nosys/mass_discriminant/DeepCsvTight/%s/mass_discriminant' % tag)).Clone()
        tags[tag] = hist.Integral()
        total += hist.Integral()
    tags['Total'] = total

    rows.append((fname, format(tags['sublead_tagged'], '.0f'),\
                format(tags['lead_tagged'], '.0f'),\
                format(tags['both_tagged'], '.0f'),\
                format(tags['both_untagged'], '.0f'),\
                format(tags['Total'], '.0f')))


rows = []
rows.append(("File", "sublead tagged", "lead tagged", "both tagged", "both untagged", "Total"))

for f in os.listdir('%s/results/%s/ctag_eff' % (project, jobid)):
    if 'data' in f:
        fname = f.replace('.root', '')
        ctageff_data_yields(fname)

print_table(rows, filename='%s/plots/%s/ctageff/ctageff_data_Yields.raw_txt' % (project, jobid))
print "\nctageff_data_Yields.raw_txt were created in\n"


print "%s/plots/%s/ctageff" % (project, jobid)
