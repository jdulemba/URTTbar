import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from pdb import set_trace

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from IDMuon import IDMuon
from IDElectron import IDElectron


class IDJet(Collection):

    def __init__(self, event):
        self.event = event
        self.jets = Collection(event, "Jet")

        self.nJets = event.nJet
        self.rho = event.fixedGridRhoFastjetAll

    def PFIsoDb(self):
        return [self.jets[i].pfRelIso04_all*self.jets[i].pt for i in range(self.nJets)]

    def RelPFIsoDb(self):
        return [self.jets[i].pfRelIso04_all for i in range(self.nJets)]

    def isTight(self):
        return [self.jets[i].tightId for i in range(self.nJets)]

    def isLoose(self):
        return [self.jets[i].softId for i in range(self.nJets)]
    #set_trace()

    def CorPFIsolation2015(self):
        rho = self.rho

        corpfIsos = []
        for mu in self.jets:
            eta = abs(mu.eta)
            effarea = 0.

            if eta < 0.8: effarea = 0.0913
            elif eta < 1.3: effarea = 0.0765
            elif eta < 2.0: effarea = 0.0546
            elif eta < 2.2: effarea = 0.0728
            elif eta < 2.5: effarea = 0.1177

            effarea *= max(rho, 0.);
            corpfiso2015 = ( mu.pfRelIso03_chg*mu.pt + max( (mu.pfRelIso03_all - mu.pfRelIso03_chg)*mu.pt - effarea, 0. ) )/mu.pt
            corpfIsos.append(corpfiso2015)

        return corpfIsos
