#! /bin/env python
from argparse import ArgumentParser
from rootpy.io import root_open

parser = ArgumentParser()
parser.add_argument('fst')
parser.add_argument('snd')
parser.add_argument('--labels')
parser.add_argument('--nolist', action='store_true')
args = parser.parse_args()

def load(fname):
	dd = {}
	with root_open(fname) as tf:
		for i in tf.sync:
			if i.RecoSuccess:
				dd[(i.Run, i.LumiSection, i.Event)] = i.MassTT
	return dd

fst = load(args.fst)
snd = load(args.snd)

l_fst, l_snd = ('first', 'second') if not args.labels else tuple(args.labels.split(','))

print 'Summary:'
print '   # events in file %s: %i' % (l_fst, len(fst))
print '   # events in file %s: %i' % (l_snd, len(snd))

sfst = set(fst.keys())
ssnd = set(snd.keys())
inter = sfst.intersection(ssnd)
print '   # events overlap: %i' % (len(inter))

massfail = []
for evt in inter:
	fmass = fst[evt]
	smass = snd[evt]	
	if abs(fmass / smass - 1) > 1e-3:
		massfail.append((evt, fmass, smass))
print '   # events mass overlap: %i' % (len(inter) - len(massfail))

if args.nolist:
	exit()

print '\n\n'
print 'Events in %s but not in %s:' % (l_fst, l_snd)
diff = sfst.difference(ssnd)
for evt in diff:
	print '  %i:%i:%i' % evt

print '\n\nEvents in %s but not in %s:' % (l_snd, l_fst)
diff = ssnd.difference(sfst)
for evt in diff:
	print '  %i:%i:%i' % evt

print '\n\nEvents with mass difference:'
for evt, i, j in massfail:
	lab = '%i:%i:%i' % evt
	print '  %s\t\t (%.2f vs %.2f)' % (lab, i, j)

