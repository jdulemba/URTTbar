import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from pdb import set_trace

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class IDMuon(Collection):

    def __init__(self, event):
        self.event = event
        self.muons = Collection(event, "Muon")

        self.nMuons = event.nMuon
        self.rho = event.fixedGridRhoFastjetAll


# create functions corresponding to all properties found
# https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc94X_doc.html

    def charge(self):
        return [self.muons[i].charge for i in range(self.nMuons)]

    def cleanmask(self):
        return [self.muons[i].cleanmask for i in range(self.nMuons)]

    def dxy(self):
        return [self.muons[i].dxy for i in range(self.nMuons)]

    def dxyErr(self):
        return [self.muons[i].dxyErr for i in range(self.nMuons)]

    def dz(self):
        return [self.muons[i].dz for i in range(self.nMuons)]

    def dzErr(self):
        return [self.muons[i].dzErr for i in range(self.nMuons)]

    def eta(self):
        return [self.muons[i].eta for i in range(self.nMuons)]

    def genPartFlav(self):
        return [self.muons[i].genPartFlav for i in range(self.nMuons)]

    def genPartIdx(self):
        return [self.muons[i].genPartIdx for i in range(self.nMuons)]

    def highPtId(self):
        return [self.muons[i].highPtId for i in range(self.nMuons)]

    def inTimeMuon(self):
        return [self.muons[i].inTimeMuon for i in range(self.nMuons)]

    def ip3d(self):
        return [self.muons[i].ip3d for i in range(self.nMuons)]

    def isGlobal(self):
        return [self.muons[i].isGlobal for i in range(self.nMuons)]

    def isPFcand(self):
        return [self.muons[i].isPFcand for i in range(self.nMuons)]

    def isTracker(self):
        return [self.muons[i].isTracker for i in range(self.nMuons)]

    def jetIdx(self):
        return [self.muons[i].jetIdx for i in range(self.nMuons)]

    def jetRelIso(self):
        return [self.muons[i].jetRelIso for i in range(self.nMuons)]

    def mass(self):
        return [self.muons[i].mass for i in range(self.nMuons)]

    def mediumId(self):
        return [self.muons[i].mediumId for i in range(self.nMuons)]

    def mediumPromptId(self):
        return [self.muons[i].mediumPromptId for i in range(self.nMuons)]

    def miniIsoId(self):
        return [self.muons[i].miniIsoId for i in range(self.nMuons)]

    def miniPFRelIso_all(self):
        return [self.muons[i].miniPFRelIso_all for i in range(self.nMuons)]

    def miniPFRelIso_chg(self):
        return [self.muons[i].miniPFRelIso_chg for i in range(self.nMuons)]

    def multiIsoId(self):
        return [self.muons[i].multiIsoId for i in range(self.nMuons)]

    def mvaId(self):
        return [self.muons[i].mvaId for i in range(self.nMuons)]

    def mvaTTH(self):
        return [self.muons[i].mvaTTH for i in range(self.nMuons)]

    def nStations(self):
        return [self.muons[i].nStations for i in range(self.nMuons)]

    def nTrackerLayers(self):
        return [self.muons[i].nTrackerLayers for i in range(self.nMuons)]

    def pdgId(self):
        return [self.muons[i].pdgId for i in range(self.nMuons)]

    def pfIsoId(self):
        return [self.muons[i].pfIsoId for i in range(self.nMuons)]

    def pfRelIso03_all(self):
        return [self.muons[i].pfRelIso03_all for i in range(self.nMuons)]

    def pfRelIso03_chg(self):
        return [self.muons[i].pfRelIso03_chg for i in range(self.nMuons)]

    def pfRelIso04_all(self):
        return [self.muons[i].pfRelIso04_all for i in range(self.nMuons)]

    def phi(self):
        return [self.muons[i].phi for i in range(self.nMuons)]

    def pt(self):
        return [self.muons[i].pt for i in range(self.nMuons)]

    def ptErr(self):
        return [self.muons[i].ptErr for i in range(self.nMuons)]

    def segmentComp(self):
        return [self.muons[i].segmentComp for i in range(self.nMuons)]

    def sip3d(self):
        return [self.muons[i].sip3d for i in range(self.nMuons)]

    def softId(self):
        return [self.muons[i].softId for i in range(self.nMuons)]

    def softMvaId(self):
        return [self.muons[i].softMvaId for i in range(self.nMuons)]

    def tightCharge(self):
        return [self.muons[i].tightCharge for i in range(self.nMuons)]

    def tightId(self):
        return [self.muons[i].tightId for i in range(self.nMuons)]

    def tkIsoId(self):
        return [self.muons[i].tkIsoId for i in range(self.nMuons)]

    def triggerIdLoose(self):
        return [self.muons[i].triggerIdLoose for i in range(self.nMuons)]


# create functions not corresponding to properties found
# https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc94X_doc.html

    def PFIsoDb(self):
        return [ pfRelIso04_all*pt for pfRelIso04_all, pt in zip( self.pfRelIso04_all(), self.pt() ) ]

    def RelPFIsoDb(self):
        return self.pfRelIso04_all()

    def isTight(self):
        return self.tightId()

    def isLoose(self):
        return self.softId()
    #set_trace()

    def CorPFIsolation2015(self):
        rho = self.rho

        corpfIsos = []
        for mu in self.muons:
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
