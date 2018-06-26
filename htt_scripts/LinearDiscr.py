import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colors
from scipy import linalg
import root_numpy as rtnp
import numpy as np
import os, glob, sys
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

#import argparse
#parser = argparse.ArgumentParser(description='Choose a ttJets sample.')
#parser.add_argument('Sample', help='Choose between ttJetsM700 and ttJetsM1000.')
#args = parser.parse_args()

project = os.environ['URA_PROJECT']
jobid = os.environ['jobid']

Merged_TreatMerged_THadPt = rtnp.root2array('%s/ttJetsM1000.ttbar_reco_3J.test.root' % project, treename="Merged_TreatMerged", branches="M_TM_THad_Pt")

ThadPt = plt.figure()
plt.hist(Merged_TreatMerged_THadPt)

