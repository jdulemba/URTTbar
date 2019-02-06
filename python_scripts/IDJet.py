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

    def area(self):
        return [self.jets[i].area for i in range(self.nJets)]

    def bRegCorr(self):
        return [self.jets[i].bRegCorr for i in range(self.nJets)]

    def bRegRes(self):
        return [self.jets[i].bRegRes for i in range(self.nJets)]

    def btagCMVA(self):
        return [self.jets[i].btagCMVA for i in range(self.nJets)]

    def btagCSVV2(self):
        return [self.jets[i].btagCSVV2 for i in range(self.nJets)]

    def btagDeepB(self):
        return [self.jets[i].btagDeepB for i in range(self.nJets)]

    def btagDeepC(self):
        return [self.jets[i].btagDeepC for i in range(self.nJets)]

    def btagDeepFlavB(self):
        return [self.jets[i].btagDeepFlavB for i in range(self.nJets)]

    def chEmEF(self):
        return [self.jets[i].chEmEF for i in range(self.nJets)]

    def chHEF(self):
        return [self.jets[i].chHEF for i in range(self.nJets)]

    def cleanmask(self):
        return [self.jets[i].cleanmask for i in range(self.nJets)]

    def electronIdx1(self):
        return [self.jets[i].electronIdx1 for i in range(self.nJets)]

    def electronIdx2(self):
        return [self.jets[i].electronIdx2 for i in range(self.nJets)]

    def eta(self):
        return [self.jets[i].eta for i in range(self.nJets)]

    def genJetIdx(self):
        return [self.jets[i].genJetIdx for i in range(self.nJets)]

    def hadronFlavour(self):
        return [self.jets[i].hadronFlavour for i in range(self.nJets)]

    def jetId(self):
        return [self.jets[i].jetId for i in range(self.nJets)]

    def mass(self):
        return [self.jets[i].mass for i in range(self.nJets)]

    def muEF(self):
        return [self.jets[i].muEF for i in range(self.nJets)]

    def muonIdx1(self):
        return [self.jets[i].muonIdx1 for i in range(self.nJets)]

    def muonIdx2(self):
        return [self.jets[i].muonIdx2 for i in range(self.nJets)]

    def nConstituents(self):
        return [self.jets[i].nConstituents for i in range(self.nJets)]

    def nElectrons(self):
        return [self.jets[i].nElectrons for i in range(self.nJets)]

    def nMuons(self):
        return [self.jets[i].nMuons for i in range(self.nJets)]

    def neEmEF(self):
        return [self.jets[i].neEmEF for i in range(self.nJets)]

    def neHEF(self):
        return [self.jets[i].neHEF for i in range(self.nJets)]

    def partonFlavour(self):
        return [self.jets[i].partonFlavour for i in range(self.nJets)]

    def phi(self):
        return [self.jets[i].phi for i in range(self.nJets)]

    def pt(self):
        return [self.jets[i].pt for i in range(self.nJets)]

    def puId(self):
        return [self.jets[i].puId for i in range(self.nJets)]

    def qgl(self):
        return [self.jets[i].qgl for i in range(self.nJets)]

    def rawFactor(self):
        return [self.jets[i].rawFactor for i in range(self.nJets)]

    #def RelPFIsoDb(self):
    #    return [self.jets[i].pfRelIso04_all for i in range(self.nJets)]

    #def isTight(self):
    #    return [self.jets[i].tightId for i in range(self.nJets)]

    #def isLoose(self):
    #    return [self.jets[i].softId for i in range(self.nJets)]
    ##set_trace()

    #def CorPFIsolation2015(self):
    #    rho = self.rho

    #    corpfIsos = []
    #    for mu in self.jets:
    #        eta = abs(mu.eta)
    #        effarea = 0.

    #        if eta < 0.8: effarea = 0.0913
    #        elif eta < 1.3: effarea = 0.0765
    #        elif eta < 2.0: effarea = 0.0546
    #        elif eta < 2.2: effarea = 0.0728
    #        elif eta < 2.5: effarea = 0.1177

    #        effarea *= max(rho, 0.);
    #        corpfiso2015 = ( mu.pfRelIso03_chg*mu.pt + max( (mu.pfRelIso03_all - mu.pfRelIso03_chg)*mu.pt - effarea, 0. ) )/mu.pt
    #        corpfIsos.append(corpfiso2015)

    #    return corpfIsos
