'''
This script prepares the combined SF file to be used for lepton ID.
Given each time you have to do some different shit because streamlining is a concept
that does not permeates in this world, everything is hardcoded.

Enjoy
'''
from rootpy.io import root_open, File
from rootpy.plotting import Hist, Hist2D
from pdb import set_trace
from ROOT import TObjString, TFile

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
trig = TFile.Open('EfficienciesAndSF_RunBtoF_Nov17Nov2017.root')
lepid = TFile.Open('Run2018ABCD_SF_ID.root')
iso = TFile.Open('Run2018ABCD_SF_ISO.root') #tracking eff as iso given is a 2D plot
#trk = root_open('ratios.root')

#trg = transpose(trig.Ele32_eta2p1_WPTight_Gsf__EffData)
#trg = transpose(trig.Get('IsoMu27_PtEtaBins/pt_abseta_ratio'))
# 
# htrk = graph2hist(trk.ratio_eta)

#read itself and dump
code = TObjString(open('merge_sf.py').read())

info = Hist(3, 0, 3, type='I')
info[0].value = 1 #0 pt as Y, 1 pt as X
info[1].value = 1 #trig SF in |eta| (1) of full eta (0)
info[2].value = 1 #ID SF in |eta| (1) of full eta (0)
info[3].value = 1 #Iso SF in |eta| (1) of full eta (0)
#info[4].value = 1 #tracking SF in |eta| (1) of full eta (0)
with File.open('output.root', 'w') as out:
#with TFile.Open('output.root', 'w') as out:
	set_trace()
	#out.WriteTObject(fill_oflow(trig.Get('IsoMu27_PtEtaBins/pt_abseta_ratio')), 'trg')
	#out.WriteTObject(fill_oflow(lepid.NUM_TightID_DEN_genTracks_pt_abseta.Clone()), 'id')
	#out.WriteTObject(fill_oflow(iso.NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta.Clone()), 'iso')
	out.WriteTObject(trig.Get('IsoMu27_PtEtaBins/pt_abseta_ratio').Clone(), 'trg')
	out.WriteTObject(lepid.Get('NUM_TightID_DEN_genTracks_pt_abseta').Clone(), 'id')
	out.WriteTObject(iso.Get('NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta').Clone(), 'iso')
	# out.WriteTObject(htrk, 'trk')
	out.WriteTObject(info, 'info')
	out.WriteTObject(code, 'code')
