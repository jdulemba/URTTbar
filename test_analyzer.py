#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from pdb import set_trace

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from python_scripts.IDMuon import IDMuon
from python_scripts.IDElectron import IDElectron
from python_scripts.IDJet import IDJet


import ast
from argparse import ArgumentParser
parser = ArgumentParser(description=__doc__)
#parser.add_argument('--flist', type=str, help='list containing the'
#                            ' location of root files to be processed')
parser.add_argument('--output', type=str,
                          help='output file name')

parser.add_argument('--isMC', type=str,
                          help='input files are data (false) or MC (true)')
args = parser.parse_args()

isMC = True if args.isMC == 'True' else False

class ExampleAnalysis(Module):
    def __init__(self):
	self.writeHistFile=True

    def beginJob(self,histFile=None,histDirName=None):
	Module.beginJob(self,histFile,histDirName)

	self.event_hist=ROOT.TH1F('Nevents', 'Nevents', 5, 0, 5)
        self.addObject(self.event_hist)
	self.genweight_hist=ROOT.TH1F('genWeight', 'genWeight', 100, 0., 5.)
        self.addObject(self.genweight_hist)

	if isMC:
	    self.pu_hist=ROOT.TH1F('pu', 'pu', 100, 0, 100)
            self.addObject(self.pu_hist)

    def analyze(self, event):

        muons = IDMuon(event)
        electrons = IDElectron(event)
        jets = IDJet(event)
        if event.nMuon > 0:
            print 'nMuons = %i' % event.nMuon
        if event.nElectron > 0:
            print 'nElectrons = %i' % event.nElectron

        set_trace()

        self.event_hist.Fill(1.)
        if isMC:
            self.genweight_hist.Fill(event.genWeight)

            pileup = Object(event, "Pileup")
            self.pu_hist.Fill(pileup.nTrueInt)

        else:
            self.genweight_hist.Fill(1.)


        return True


#set_trace()
#files = ast.literal_eval(args.flist)
##if len(files) > 1:
##    files = files[0:2]
files=["root://cmseos.fnal.gov//store/group/lpcbtag/jdulemba/2019Jan07/ttJets/007B405B-74C1-E811-B4D7-44A84225C911.root"]


p=PostProcessor(".",files,branchsel=None,modules=[ExampleAnalysis()],noOut=True,histFileName=args.output, histDirName="plots")
p.run()
