from glob import glob
import URAnalysis.Utilities.prettyjson as prettyjson
from uncertainties import ufloat

def convert(jm):
	ret = {}
	for 

jmap = prettyjson.loads(open('yields.json').read())
#convert to ufloats
for sam in jmap:
	for cat in jmap[sam]:
		vals = jmap[sam][cat]
		jmap[sam][cat] = ufloat(*vals)
