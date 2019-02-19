import rootpy.io as io
import rootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
rootpy.log["/data_views"].setLevel(rootpy.log.WARNING)
log = rootpy.log["/make_postfit_plots"]
log.setLevel(rootpy.log.INFO)
import os
import re
import ROOT
from URAnalysis.Utilities.quad import quad
import URAnalysis.Utilities.latex as tex
from pdb import set_trace
from URAnalysis.Utilities.print_table import print_table
import URAnalysis.Utilities.prettyjson as prettyjson

def pipppo():
	pass

wps = [
	#('csvLoose'  , 'CSVv2 L'),
	#('csvMedium' , 'CSVv2 M'),
	#('csvTight'  , 'CSVv2 T'),
	('DeepctagLoose' , 'DeepCSV c-tagger L'),
	('DeepctagMedium', 'DeepCSV c-tagger M'),
	('DeepctagTight' , 'DeepCSV c-tagger T'),
	('DeepCSVLoose' , 'DeepCSV L'), 
	('DeepCSVMedium',	'DeepCSV M'),
	('DeepCSVTight' ,	'DeepCSV T'), 
	#('cmvaLoose' , 'cMVAv2 L' ),
	#('cmvaMedium', 'cMVAv2 M'),
	#('cmvaTight' , 'cMVAv2 T' ),
	]

categories = [
	'notag', 
	'leadtag', 
	'subtag',
	'ditag',
]

samples = [
	('right_whad', '\\ttbar, right W$_{\mathrm{had}}$'),
	('wrong_whad', '\\ttbar, wrong W$_{\mathrm{had}}$'),
	('nonsemi_tt', 'Other \\ttbar decay'),
	('single_top', 'single top'),
	('vjets', 'V+jets'),
	('qcd', 'QCD Multi-jet'),
]

jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']
table   = open('%s/plots/%s/ctageff/mass_discriminant/prefit_yields.tex' % (project, jobid), 'w')
results = open('%s/plots/%s/ctageff/mass_discriminant/results.tex'       % (project, jobid), 'w')
ncols = len(samples)+3
table.write(
	'''
\\begin{tabular}{%s}
\hline
''' % ('c'*ncols)
)
table.write(tex.tabline(['category', 'observed', 'total expected']+[i for _, i in samples], convert=False))
header = '\n \\begin{tabular}{cc} \n \n'
results.write('\n \\begin{tabular}{ccc} \n \\hline \n')
results.write(tex.tabline(['working point', 'charm SF', 'stat only unc.'], convert=False))
results.write('\\hline \n')

## to make txt file with results
rows = []
rows.append(("Working Point", "SF Mean", "SF Upper Error", "SF Lower Error", "Stat Upper Error", "Stat Lower Error"))
scale_factors = {}


for wp, wpname in wps:
	table.write('\hline \n')
	table.write('\multicolumn{%d}{c}{ %s } \\\\ \n' % (ncols, wp))
	table.write('\hline \n')
		
	wpdir = '%s/plots/%s/ctageff/mass_discriminant/%s' % (project, jobid, wp)
	data_f = io.root_open('%s/datacard.root'   % (wpdir))
	mc_f   = io.root_open('%s/MaxLikeFit.root' % (wpdir))
	stat_f = io.root_open('%s/MaxLikeFitStatOnly.root' % (wpdir))
	#set_trace()
	prefit = mc_f.norm_prefit
	pshapes = mc_f.shapes_prefit
	pars = rootpy.asrootpy(mc_f.fit_s.floatParsFinal())
	stat_pars = rootpy.asrootpy(stat_f.fit_s.floatParsFinal())
	results.write(
		tex.tabline(
			[
				wpname, 
				tex.format_roovar(pars['charmSF'], True), 
				tex.format_roovar(stat_pars['charmSF'], True, erronly=True)
				], convert=False
			)
		)

	if type(pars['charmSF'].error)==float:
		rows.append((wpname, format(pars['charmSF'].value, '.3f'),\
			format(pars['charmSF'].error, '.3f'),\
			format(pars['charmSF'].error, '.3f'),\
			format(stat_pars['charmSF'].error[0], '.3f'),\
			format(stat_pars['charmSF'].error[1], '.3f'))
		)
		scale_factor = {
			wpname : {
				'mean' : pars['charmSF'].value,
				'sf_upper' : pars['charmSF'].error,
				'sf_lower' : pars['charmSF'].error,
				'stat_upper' : stat_pars['charmSF'].error[0],
				'stat_lower' : stat_pars['charmSF'].error[1],
			}
		}
	else:
		rows.append((wpname, format(pars['charmSF'].value, '.3f'),\
			format(pars['charmSF'].error[0], '.3f'),\
			format(pars['charmSF'].error[1], '.3f'),\
			format(stat_pars['charmSF'].error[0], '.3f'),\
			format(stat_pars['charmSF'].error[1], '.3f'))
		)
		scale_factor = {
			wpname : {
				'mean' : pars['charmSF'].value,
				'sf_upper' : pars['charmSF'].error[0],
				'sf_lower' : pars['charmSF'].error[1],
				'stat_upper' : stat_pars['charmSF'].error[0],
				'stat_lower' : stat_pars['charmSF'].error[1],
			}
		}

	scale_factors[wp] = scale_factor

	
	for cat in categories:
		line = [cat, '%d' % data_f.Get('%s/data_obs' % cat).Integral()]
		total = pshapes.Get('%s/total' % cat)
		t_err = ROOT.Double()
		t_val = total.IntegralAndError(1, total.nbins(), t_err)
		line.append(tex.format_with_error(t_val, t_err))
		for sam, _ in samples:
			path = '%s/%s' % (cat, sam)
			#if wp == 'ctagMedium' and cat == 'notag' and sam == 'qcd':
			#   set_trace()
			if path not in prefit:
				print path, "in", wp, "does not exist!"
				line.append('$0 \pm 0$')
			else:
				var = prefit['%s/%s' % (cat, sam)]
				shape = pshapes.Get(path)
				line.append(
					tex.format_roovar(var)
					)
		table.write(tex.tabline(line, convert=False))
		

table.write('''\hline
\end{tabular}
''')
results.write('''\hline
\end{tabular}
''')

with open('%s/plots/%s/ctageff/mass_discriminant/results.json' % (project, jobid), 'w') as f:
    f.write(prettyjson.dumps(scale_factors))

print_table(rows, filename='%s/plots/%s/ctageff/mass_discriminant/results.txt' % (project, jobid))
