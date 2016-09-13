'''
This script prepares the combined SF file to be used for lepton ID.
Given each time you have to do some different shit because streamlining is a concept
that does not permeates in this world, everything is hardcoded.

Enjoy
'''
from rootpy.io import root_open
from rootpy.plotting import Hist

trig = root_open()
lepid = root_open()
iso = root_open()
trk = root_open()

info = Hist(3, 0, 3, type='I')
info[0].value = 1 #0 pt as Y, 1 pt as X
info[1].value = 0 #trig SF in |eta| (1) of full eta (0)
info[2].value = 0 #ID SF in |eta| (1) of full eta (0)
info[3].value = 0 #Iso SF in |eta| (1) of full eta (0)
info[4].value = 0 #tracking SF in |eta| (1) of full eta (0)
with root_open('output.root', 'w') as out:
	out.WriteTObject(trig.some_histo.Clone(), 'trg')
	out.WriteTObject(lepid.some_histo.Clone(), 'id')
	out.WriteTObject(iso.some_histo.Clone(), 'iso')
	out.WriteTObject(trk.some_histo.Clone(), 'trk')
	out.WriteTObject(info, 'info')
