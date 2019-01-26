import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from pdb import set_trace

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
import math
from enum import Enum

class IDElectron(Collection):

    def __init__(self, event):
    #def __init__(self, event, prefix, lenVar=None):
    #    Collection.__init__(self, event, prefix, lenVar=None)
        self.event = event
        self.electrons = Collection(event, "Electron")

        self.nElectrons = event.nElectron
        self.rho = event.fixedGridRhoFastjetAll

        self.IDS = Enum('IDS', 'FAIL TIGHT_15 MEDIUM_15 VETO_15 TIGHT_15_NoECAL_Gap NOTVETO_15 FAKES')
        self.id_names = { i.name : i.value for i in self.IDS }
        self.USEISO = True

    def etaSC(self):
        return [self.electrons[i].deltaEtaSC+self.electrons[i].eta for i in range(self.nElectrons)]

    def IPCuts(self):
        #etaSC = self.etaSC()
        IPCuts = []
        for el in self.electrons:
            etaSC =  el.deltaEtaSC+el.eta
            dxy = math.fabs(el.dxy)
            dz = math.fabs(el.dz)

            #print 'etaSC=%f, dxy=%f, dz=%f' % (etaSC, dxy, dz)

            if etaSC < 1.479:
                dxy_dz_pass = True if (dxy < 0.05 and dz < 0.10) else False
                #print 'etaSC < 1.479, dxy < 0.05 and dz < 0.10: ', dxy_dz_pass 
                IPCuts.append(dxy_dz_pass)
            else:
                dxy_dz_pass = True if (dxy < 0.10 and dz < 0.20) else False
                #print 'etaSC < 1.479, dxy < 0.10 and dz < 0.20: ', dxy_dz_pass 
                IPCuts.append(dxy_dz_pass)

        return IPCuts

    #def PFIsoDb(self):
    #    return [self.electrons[i].pfRelIso04_all*self.electrons[i].pt for i in range(self.nElectrons)]

    #def RelPFIsoDb(self):
    #    return [self.electrons[i].pfRelIso04_all for i in range(self.nElectrons)]

    def VetoID25ns(self):
        return [self.electrons[i].cutBased == 1 for i in range(self.nElectrons)]

    def LooseID25ns(self):
        IPCuts = self.IPCuts()
        return [IPCuts[i] and self.electrons[i].cutBased == 2 for i in range(self.nElectrons)]

    def MediumID25ns(self):
        IPCuts = self.IPCuts()
        return [IPCuts[i] and self.electrons[i].cutBased == 3 for i in range(self.nElectrons)]

    def TightID25ns(self):
        IPCuts = self.IPCuts()
        return [IPCuts[i] and self.electrons[i].cutBased == 4 for i in range(self.nElectrons)]

    def TightID25ns_noECAL_Gap(self):
        tightID = self.TightID25ns()
        absetaSC = [math.fabs(etasc) for etasc in self.etaSC()]
        return [tightID[i] and (absetaSC[i] <= 1.4442 or absetaSC[i] >= 1.5660) for i in range(self.nElectrons)]

    def FakeID(self):
        etaSC = self.etaSC()
        abseta = [math.fabs(i) for i in self.etaSC()]
        ecalgap = [(i <= 1.4442 or i >= 1.5660) for i in abseta]
        IPCuts = self.IPCuts()
        PFIsoRho2015 = self.PFIsolationRho2015()

        Iso = []
        for i in range(self.nElectrons):
            if etaSC[i] < 1.479:
                iso = True if PFIsoRho2015[i] >= 0.0588 else False
                Iso.append(iso)
            else:
                iso = True if PFIsoRho2015[i] >= 0.0571 else False
                Iso.append(iso)

        return [IPCuts[i] and ecalgap[i] and Iso[i] for i in range(self.nElectrons)]

    def PFIsolationRho2015(self):
    # Isolation is calculated following an example in [1], as recommended in [2]
    #[1] https://github.com/cms-sw/cmssw/blob/CMSSW_7_6_2/RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleEffAreaPFIsoCut.cc#L75-L90
    #[2] https://hypernews.cern.ch/HyperNews/CMS/get/egamma/1664/1.html
    # The following values refer to EA for cone 0.3 and fixedGridRhoFastjetAll. 
    # They are valid for electrons only, different EA are available for muons.
    #EA Values available here: https://github.com/ikrav/cmssw/blob/egm_id_80X_v1/RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt

        rho = self.rho
        pfIsos = []

        for el in self.electrons:
            abseta = math.fabs(el.deltaEtaSC+el.eta)
            effarea = 0.

            if (abseta >= 0.0000 and abseta <= 1.0000): effarea = 0.1703
            if (abseta >  1.0000 and abseta <= 1.4790): effarea = 0.1715
            if (abseta >  1.4790 and abseta <= 2.0000): effarea = 0.1213
            if (abseta >  2.0000 and abseta <= 2.2000): effarea = 0.1230
            if (abseta >  2.2000 and abseta <= 2.3000): effarea = 0.1635
            if (abseta >  2.3000 and abseta <= 2.4000): effarea = 0.1937
            if (abseta >  2.4000 and abseta <= 5.0000): effarea = 0.2393

            effarea *= max(rho, 0.);
            pfiso2015 = el.pfRelIso03_all*el.pt + max( (el.pfRelIso03_all - el.pfRelIso03_chg)*el.pt - effarea, 0. )
            pfIsos.append(pfiso2015)

        return pfIsos

    def ID(self, idtyp):
        etas = [self.electrons[i].eta for i in range(self.nElectrons)]
        pass_eta = [ eta <= 2.5 for eta in etas ]

        idvals = []
        if idtyp == 'FAIL': idvals = [False for i in range(self.nElectrons)]
        elif idtyp == 'TIGHT_15': idvals = self.TightID25ns()
        elif idtyp == 'MEDIUM_15': idvals = self.MediumID25ns()
        elif idtyp == 'LOOSE_15': idvals = self.LooseID25ns()
        elif idtyp == 'VETO_15': idvals = self.VetoID25ns()
        elif idtyp == 'TIGHT_15_NoECAL_Gap': idvals = self.TightID25ns_noECAL_Gap()
        ## not sure how to deal with !Veto   elif idtyp == 'NOTVETO_15': idvals = !self.VetoID25ns()
        elif idtyp == 'FAKES': idvals = self.FakeID()
        else: idvals = [False for i in range(self.nElectrons)]

        return [eta and idval for eta, idval in zip(idvals, pass_eta) ]

