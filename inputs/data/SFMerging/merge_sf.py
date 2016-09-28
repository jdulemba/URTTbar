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
		idx = ret.FindFixBin(i.x.center, i.y.center)
		ret[idx].error = i.error
		ret[idx].value = i.value
	return ret

def graph2hist(graph):
	xs = [i for i, _ in graph]
	widths = [i for i in trk.ratio_eta.xerr()]
	lows = [x-w[0] for x, w in zip(xs, widths)]
	lows.append(widths[-1][1]+xs[-1])
	ret = Hist(lows)
	for i, xy in enumerate(graph):
		_, y = xy
		ret[i+1].value = y

trig = root_open('TriggerSF_v1.root')
lepid = root_open('ideff.root')
iso = root_open('tracking_eff.root') #tracking eff as iso given is a 2D plot
#trk = root_open('ratios.root')

trg = transpose(trig.Ele32_eta2p1_WPTight_Gsf__EffData)
# trg1 = trig.IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093.efficienciesDATA.abseta_pt_DATA 
# trg2 = trig.IsoMu22_OR_IsoTkMu22_PtEtaBins_Run274094_to_276097.efficienciesDATA.abseta_pt_DATA
# trg = trg1*0.0482 + trg2*0.9517
# 
# htrk = graph2hist(trk.ratio_eta)

#read itself and dump
code = TObjString(open('merge_sf.py').read())

info = Hist(3, 0, 3, type='I')
info[0].value = 0 #0 pt as Y, 1 pt as X
info[1].value = 0 #trig SF in |eta| (1) of full eta (0)
info[2].value = 0 #ID SF in |eta| (1) of full eta (0)
info[3].value = 0 #Iso SF in |eta| (1) of full eta (0)
info[4].value = 0 #tracking SF in |eta| (1) of full eta (0)
with root_open('output.root', 'w') as out:
	out.WriteTObject(trg, 'trg')
	out.WriteTObject(lepid.EGamma_SF2D.Clone(), 'id')
	out.WriteTObject(iso.EGamma_SF2D.Clone(), 'iso')
	# out.WriteTObject(htrk, 'trk')
	out.WriteTObject(info, 'info')
	out.WriteTObject(code, 'code')
