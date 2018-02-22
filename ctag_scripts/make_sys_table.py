import os, sys
from rootpy.io import root_open
from glob import glob
from argparse import ArgumentParser
from rootpy import asrootpy
from math import sqrt
from pdb import set_trace

parser = ArgumentParser()
parser.add_argument('wp')
parser.add_argument('--eras', default='All')
args = parser.parse_args()

if not (args.eras == 'All_Runs' or args.eras == 'Run_B' or args.eras == "Run_CtoE" or args.eras == "Run_EtoF"):
    logging.error('Not a valid era to choose from.')
    sys.exit()


def get_unc(fname):
	tf = root_open(fname)
	sf = asrootpy(tf.fit_s.floatParsFinal())['charmSF']
	unc = max(abs(i) for i in sf.error)
	return unc

jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']
indir = '%s/plots/%s/ctageff/%s/mass_discriminant/%s/' % (project, jobid, args.eras, args.wp)

total = get_unc('%s/MaxLikeFit.root' % indir)
stat = get_unc('%s/MaxLikeFitStatOnly.root' % indir)
shifts = glob('%s/sys_breakdown/*.root' % indir)
values = []

for shift in shifts:
	name = os.path.basename(shift).split('.')[0]
	if name == 'other': continue
	shifted = get_unc(shift)
	if shifted < stat:
		print name, 'has smaller uncertainties than expected'
		effect = 0
	else:
		effect = 100*sqrt(shifted**2 - stat**2)	
	values.append((name, effect))

values.sort(key=lambda x:x[1], reverse=True)

with open('%s/systematics.tex' % indir, 'w') as tab:
	tab.write(r'''\begin{tabular}{cc}
\hline
Source & effect \\
\hline
''')
	for info in values:
		tab.write('%s & %.1f\\%% \\\\ \n' % info)
	tab.write('\\hline \n')
	tab.write('Total & %.1f\\%% \\\\ \n' % (total*100))
	tab.write('\\hline \n')
	tab.write('\\end{tabular} \n')
