#! /bin/env python

import os
from glob import glob
import URAnalysis.Utilities.prettyjson as prettyjson

jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']
infiles = glob('%s/inputs/%s/*.txt' % (project, jobid))
names = [os.path.basename(i).replace('.txt','') for i in infiles]
data = [i for i in names if i.startswith('data')]
mcs  = [i for i in names if not i.startswith('data')]
samples = prettyjson.loads(
	open('%s/samples.json' % project).read()
)

table = open('%s/samples.tex' % project, 'w')
table.write(
	r'''
\begin{tabular}{cc}
\hline
\multicolumn{2}{c}{ Simulated } \\
\hline
Sample & $\sigma$ [pb] \\
\hline 
'''
)

for name in mcs:
	info = [i for i in samples if i['name'] == name][0]
	dbs = info['DBSName'].split('/')[1]
	dbs = dbs.replace('_',r'\_')
	xsec = info['xsection']
	table.write('%s & %.0f \\\\ \n' % (dbs, xsec))

table.write(r'''\hline
\multicolumn{2}{c}{ Real data } \\
\hline 
Sample & luminosity [pb$^{-1}$] \\ 
\hline 
''')

for name in data:
	info = [i for i in samples if i['name'] == name][0]
	dbs = info['DBSName'].split('/')
	dbs = '%s %s' % (dbs[1], dbs[2])
	dbs = dbs.replace('_',r'\_')
	lumi = float(
		open('%s/inputs/%s/%s.lumi' % (project, jobid, name)
				 ).read()
		)
	table.write('%s & %.0f \\\\ \n' % (dbs, lumi))
table.write(
	r'''\hline
\end{tabular}
''')

