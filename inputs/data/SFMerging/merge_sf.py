'''
This script prepares the combined SF file to be used for lepton ID.
Given each time you have to do some different shit because streamlining is a concept
that does not permeates in this world, everything is hardcoded.

Enjoy
'''
from rootpy.io import root_open
from rootpy.plotting import Hist, Hist2D
from pdb import set_trace
from ROOT import TObjString
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('lepton', help='Choose between electron or muon.')
args = parser.parse_args()

lchoices = ['muon', 'electron']
if args.lepton not in lchoices:
    raise RuntimeError("Can only choose 'muon' or 'electron'!!")

def transpose(hist2d):
	xbins = set()
	ybins = set()
	for i in hist2d:
		xbins.add(i.x.low)
		xbins.add(i.x.high)
		ybins.add(i.y.low)
		ybins.add(i.y.high)

	xbins = sorted(list(xbins))[1:-1]
	ybins = sorted(list(ybins))[1:-1]
	ret = Hist2D(ybins, xbins)
	for i in hist2d:
		idx = ret.FindFixBin(i.y.center, i.x.center)
		ret[idx].error = i.error
		ret[idx].value = i.value
	return ret

def fill_oflow(hist2d):
	xbins = hist2d.nbins()
	ybins = hist2d.nbins(1)
	for i in range(xbins+2):
		hist2d[i, 0].value = hist2d[i, 1].value
		hist2d[i, 0].error = hist2d[i, 1].error
		hist2d[i, ybins+1].value = hist2d[i, ybins].value
		hist2d[i, ybins+1].error = hist2d[i, ybins].error
		
	for i in range(ybins+2):
		hist2d[0, i].value = hist2d[1, i].value
		hist2d[0, i].error = hist2d[1, i].error
		hist2d[xbins+1, i].value = hist2d[xbins, i].value
		hist2d[xbins+1, i].error = hist2d[xbins, i].error
	return hist2d

def graph2hist(graph):
	xs = [i for i, _ in graph]
	widths = [i for i in trk.ratio_eta.xerr()]
	lows = [x-w[0] for x, w in zip(xs, widths)]
	lows.append(widths[-1][1]+xs[-1])
	ret = Hist(lows)
	for i, xy in enumerate(graph):
		_, y = xy
		ret[i+1].value = y

#set_trace()
trg = 0
lepid = 0
lepiso = 0

info = Hist(3, 0, 3, type='I')
if args.lepton == 'muon':
    lid = root_open('Muons/Run2018ABCD_SF_ID.root')
    lepid = fill_oflow(lid.NUM_TightID_DEN_genTracks_pt_abseta.Clone())

    iso = root_open('Muons/Run2018ABCD_SF_ISO.root') #tracking eff as iso given is a 2D plot
    lepiso = fill_oflow(iso.NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta.Clone())

        ## weights are hard coded, need to check!
    trig_p1 = root_open('Muons/EfficienciesAndSF_2018Data_BeforeMuonHLTUpdate.root')
    trig_p2 = root_open('Muons/EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root')
    print '\n       Weights used for muon trigger effs are hardcoded. Need to check!!!\n'
    # Muon run splitting found here: https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffs2018
    p1_weight = 0.15105962910078496 # lumi ratio before 316361
    p2_weight = 0.8489403708992149 # lumi ratio after and including 316361
    trg = trig_p1.IsoMu24_PtEtaBins.pt_abseta_ratio*p1_weight + trig_p2.IsoMu24_PtEtaBins.pt_abseta_ratio*p2_weight

    info[0].value = 1 #0 pt as Y, 1 pt as X
    info[1].value = 1 #trig SF in |eta| (1) of full eta (0)
    info[2].value = 1 #ID SF in |eta| (1) of full eta (0)
    info[3].value = 1 #Iso SF in |eta| (1) of full eta (0)

if args.lepton == 'electron':
    lid = root_open('Electrons/2018_ElectronTight.root')
    lepid = fill_oflow(lid.EGamma_SF2D.Clone())

    trig = root_open('Electrons/electron_sf_Moriond17_TightCutID.root')
    trg = trig.trg.Clone()

    info[0].value = 0 #0 pt as Y, 1 pt as X
    info[1].value = 0 #trig SF in |eta| (1) of full eta (0)
    info[2].value = 0 #ID SF in |eta| (1) of full eta (0)
    info[3].value = 0 #Iso SF in |eta| (1) of full eta (0)

#set_trace()
if trg == 0 or lepid == 0:
    raise RuntimeError('Lepton trigger AND id info must be provided!')

#read itself and dump
code = TObjString(open('merge_sf.py').read())

with root_open('output.root', 'w') as out:
	#set_trace()
	out.WriteTObject(trg, 'trg')
	out.WriteTObject( lepid, 'id')
	if args.lepton == 'muon': out.WriteTObject( lepiso, 'iso')
	out.WriteTObject(info, 'info')
	out.WriteTObject(code, 'code')
