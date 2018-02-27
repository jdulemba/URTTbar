import rootpy.io as io
import rootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
rootpy.log["/data_views"].setLevel(rootpy.log.WARNING)
log = rootpy.log["/make_postfit_plots"]
log.setLevel(rootpy.log.INFO)
import os,sys
import re
import ROOT
from URAnalysis.Utilities.quad import quad
import URAnalysis.Utilities.latex as tex
from pdb import set_trace
import logging
import URAnalysis.Utilities.prettyjson as prettyjson
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--eras', default='All')

args = parser.parse_args()

if args.eras == 'All':
    dir_name = 'All_Runs'
elif args.eras == 'B' or args.eras == 'CtoE' or args.eras == 'EtoF':
    dir_name = 'Run_%s' % args.eras
else:
    logging.error('Not a valid era to choose from.')
    sys.exit()


def pipppo():
	pass

wps = [
	('csvLoose'  , 'CSVv2 L'),
	('csvMedium' , 'CSVv2 M'),
	('csvTight'  , 'CSVv2 T'),
	('ctagLoose' , 'c-tagger L'),
	('ctagMedium', 'c-tagger M'),
	('ctagTight' , 'c-tagger T'),
	('DeepCsvLoose' , 'DeepCSV L'), 
	('DeepCsvMedium',	'DeepCSV M'),
	('DeepCsvTight' ,	'DeepCSV T'), 
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
table   = open('%s/plots/%s/ctageff/%s/mass_discriminant/prefit_yields.tex' % (project, jobid, dir_name), 'w')
results = open('%s/plots/%s/ctageff/%s/mass_discriminant/results.tex'       % (project, jobid, dir_name), 'w')
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

rows = []
rows.append(("Working Point", "SF Mean", "SF Upper Error", "SF Lower Error"))#, "Stat Upper Error", "Stat Lower Error"))

scale_factors = {}

for wp, wpname in wps:
	table.write('\hline \n')
	table.write('\multicolumn{%d}{c}{ %s } \\\\ \n' % (ncols, wp))
	table.write('\hline \n')
		
	wpdir = '%s/plots/%s/ctageff/%s/mass_discriminant/%s' % (project, jobid, dir_name, wp)
	data_f = io.root_open('%s/datacard.root'   % (wpdir))
	mc_f   = io.root_open('%s/MaxLikeFit.root' % (wpdir))
	stat_f = io.root_open('%s/MaxLikeFitStatOnly.root' % (wpdir))
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

	rows.append((wpname, format(pars['charmSF'].value, '.3f'),\
                format(pars['charmSF'].error[0], '.3f'),\
                format(pars['charmSF'].error[1], '.3f')))#,\
                format(stat_pars['charmSF'].error[0], '.3f'),\
                format(stat_pars['charmSF'].error[1], '.3f')))

	scale_factor = {
		wpname : {
			'mean' : pars['charmSF'].value,
			'sf_upper' : pars['charmSF'].error[0],
			'sf_lower' : pars['charmSF'].error[1]
			'stat_upper' : stat_pars['charmSF'].error[0],
			'stat_lower' : stat_pars['charmSF'].error[1],
            }
        }
	scale_factors[wp] = scale_factor
	#set_trace()
	
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

with open('%s/plots/%s/ctageff/%s/mass_discriminant/results.json' % (project, jobid, dir_name), 'w') as f:
    f.write(prettyjson.dumps(scale_factors))
#set_trace()

def print_table(lines, filename, separate_head=True):
    """Prints a formatted table given a 2 dimensional array"""
    #Count the column width
    widths = []
    for line in lines:
        for i,size in enumerate([len(x) for x in line]):
            while i >= len(widths):
                widths.append(0)
            if size > widths[i]:
                widths[i] = size

    #Generate the format string to pad the columns
    print_string = ""
    for i,width in enumerate(widths):
        print_string += "{" + str(i) + ":" + str(width) + "} | "
    if (len(print_string) == 0):
        return
    print_string = print_string[:-3]

    with open(filename, 'w') as f:
        #Print the actual data
        for i,line in enumerate(lines):
            print >> f, print_string.format(*line)
            if (i == 0 and separate_head):
                print >> f, "-"*(sum(widths)+3*(len(widths)-1))


print_table(rows, filename='%s/plots/%s/ctageff/%s/mass_discriminant/results.txt' % (project, jobid, dir_name))
print '\nprefit_yields.tex\nresults.tex\nresults.txt\nresults.json\n\nhave been created and are located in\n     %s/plots/%s/ctageff/%s/mass_discriminant' % (project, jobid, dir_name)
