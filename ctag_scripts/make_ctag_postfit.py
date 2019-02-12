import rootpy.io as io
import rootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
rootpy.log["/data_views"].setLevel(rootpy.log.WARNING)
log = rootpy.log["/make_postfit_plots"]
log.setLevel(rootpy.log.INFO)
from rootpy.plotting import Hist
from argparse import ArgumentParser
import os, sys
import re
import ROOT
from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.Utilities.prettyjson as prettyjson
from URAnalysis.Utilities.quad import quad
from styles import styles
from pdb import set_trace
from labels import set_pretty_label
from URAnalysis.Utilities.tables import Table
from argparse import ArgumentParser
import numpy as np

def foo():
	pass

parser = ArgumentParser()
parser.add_argument('working_point')
parser.add_argument('plots', choices=['prefit', 'postfit', 'both'],
                    help='what to plot')
parser.add_argument('--eras', default='All')

args = parser.parse_args()
plots_to_do = []
if args.plots == 'prefit' or args.plots == 'both':
	plots_to_do.append(('prefit', 'shapes_prefit'))
if args.plots == 'postfit' or args.plots == 'both':
	plots_to_do.append(('postfit', 'shapes_fit_s'))

if args.eras == 'All':
	run_dir = 'All_Runs'

elif args.eras == 'B':
	run_dir = 'Run_B'

elif args.eras == 'CtoE':
	run_dir = 'Run_CtoE'

elif args.eras == 'EtoF':
	run_dir = 'Run_EtoF'

else:
    logging.error('Not a valid era to choose from.')
    sys.exit()

def fix_binning(postfit, correct_bin):
	ret = correct_bin.Clone()
	ret.Reset()
	ret.title = postfit.title
	ret.decorate(**postfit.decorators)
	assert(ret.nbins() <= postfit.nbins())
	for rbin, pbin in zip(ret, postfit):
		rbin.value = pbin.value
		rbin.error = pbin.error
	ret.xaxis.title = '#lambda_{M}'
	return ret

input_dir = 'plots/%s/ctageff/%s/mass_discriminant/%s' % (os.environ['jobid'], run_dir, args.working_point)

mlfit_file = io.root_open(
	'%s/MaxLikeFit.root' % input_dir
	)
datacard_file = io.root_open(
	'%s/datacard.root' % input_dir
) #contains data
#set_trace()
def ordering(histo):
	multiplier = 1.
	if histo.title.startswith('right_whad'):
		multiplier = 100000.
	elif histo.title.startswith('wrong_whad'):
		multiplier = 100.
	return multiplier*histo.Integral()

if args.eras == 'All':
    plotter = BasePlotter(
    	'%s/postfit' % input_dir,
    	defaults = {'save' : {'png' : True, 'pdf' : True}},
    	styles = {
            'right_whad *' : styles['right_whad *'],
    		'qcd *' : styles['QCD*'],
    		'single_top *' : styles['single*'],
            'wrong_whad *' : styles['wrong_whad *'],
            'nonsemi_tt *' : styles['nonsemi_tt *'],
    		'vjets *' : styles['[WZ]Jets*'],
    		'data_obs' : styles['data*'],
    		'Total signal+background *' : {
    			'legendstyle' : 'f',
    			'drawstyle' : 'PE2',
    			'linewidth' : 0,
    			'title' : "uncertainty",
    			'markersize' : 0,
    			'fillcolor' : 1,
    			'fillstyle' : 3013
    			}
    		}
    )
elif args.eras == 'B':
    plotter = BasePlotter(
    	'%s/postfit' % input_dir,
    	defaults = {'save' : {'png' : True, 'pdf' : True}},
    	styles = {
            'right_whad *' : styles['right_whad *'],
    		'qcd *' : styles['QCD*'],
    		'single_top *' : styles['single*'],
            'wrong_whad *' : styles['wrong_whad *'],
            'nonsemi_tt *' : styles['nonsemi_tt *'],
    		'vjets *' : styles['[WZ]Jets*'],
    		'data_obs' : styles['data*B*'],
    		'Total signal+background *' : {
    			'legendstyle' : 'f',
    			'drawstyle' : 'PE2',
    			'linewidth' : 0,
    			'title' : "uncertainty",
    			'markersize' : 0,
    			'fillcolor' : 1,
    			'fillstyle' : 3013
    			}
    		}
    )
elif args.eras == 'CtoE':
    plotter = BasePlotter(
    	'%s/postfit' % input_dir,
    	defaults = {'save' : {'png' : True, 'pdf' : True}},
    	styles = {
            'right_whad *' : styles['right_whad *'],
    		'qcd *' : styles['QCD*'],
    		'single_top *' : styles['single*'],
            'wrong_whad *' : styles['wrong_whad *'],
            'nonsemi_tt *' : styles['nonsemi_tt *'],
    		'vjets *' : styles['[WZ]Jets*'],
    		'data_obs' : styles['data*CtoE*'],
    		'Total signal+background *' : {
    			'legendstyle' : 'f',
    			'drawstyle' : 'PE2',
    			'linewidth' : 0,
    			'title' : "uncertainty",
    			'markersize' : 0,
    			'fillcolor' : 1,
    			'fillstyle' : 3013
    			}
    		}
    )
elif args.eras == 'EtoF':
    plotter = BasePlotter(
    	'%s/postfit' % input_dir,
    	defaults = {'save' : {'png' : True, 'pdf' : True}},
    	styles = {
            'right_whad *' : styles['right_whad *'],
    		'qcd *' : styles['QCD*'],
    		'single_top *' : styles['single*'],
            'wrong_whad *' : styles['wrong_whad *'],
            'nonsemi_tt *' : styles['nonsemi_tt *'],
    		'vjets *' : styles['[WZ]Jets*'],
    		'data_obs' : styles['data*EtoF*'],
    		'Total signal+background *' : {
    			'legendstyle' : 'f',
    			'drawstyle' : 'PE2',
    			'linewidth' : 0,
    			'title' : "uncertainty",
    			'markersize' : 0,
    			'fillcolor' : 1,
    			'fillstyle' : 3013
    			}
    		}
    )

def print_var_line(varname, result):
	var  = result[varname]
	merr = min(var.error) if hasattr(var.error, '__len__') else var.error
	Merr = max(var.error) if hasattr(var.error, '__len__') else var.error
	return 'best fit value %s: %.3f %.3f/%.3f\n' % (varname, var.value, merr, Merr)

for pfit, tdir in plots_to_do:
	out_dir = '%s/%s' % (input_dir,pfit)
	try:
		shapes = mlfit_file.Get(tdir)
	except io.file.DoesNotExist:
		with open('%s/yields.raw_txt' % out_dir,'w') as tab:
			tab.write('---------  FIT FAILED  ---------\n')
		print "FIT FAILED!"
		continue
	plotter.set_outdir(out_dir)
	categories = [i.name for i in shapes.keys()]
	sample_names = set()
	for c in categories:
		for i in shapes.Get(c).keys():
			if i.name.startswith('total'): continue
			sample_names.add(i.name)
	sample_names = list(sample_names)
	cols = ['category', 'observed']+sample_names
	cols = [i+':%'+str(len(i)+3)+'s' for i in cols]
	table = Table(*cols, title='%s %s yields' % (args.working_point, pfit))

	for cat_name in categories:
		log.info("making plots for %s" % cat_name)
		data_dir = datacard_file.Get(cat_name)
		cat_dir = shapes.Get(cat_name)
		data = data_dir.data_obs
		data.title = 'data_obs' 
		hsum = fix_binning(cat_dir.total, data)
		available_samples = [i.name for i in cat_dir.keys() if not i.name.startswith('total')]

		line = [cat_name, '%.1f' % data.Integral()]+['%.1f' % (cat_dir.Get(i).Integral() if i in available_samples else 0.0) for i in sample_names]
		table.add_line(*line)

		#set_trace()
		samples = [fix_binning(cat_dir.Get(i), data) for i in available_samples]
		samples.sort(key=ordering)

		stack = plotter.create_stack(*samples, sort=False)
		legend = LegendDefinition(position='NE')
		plotter.overlay_and_compare( 
			#set_trace(),
			[stack, hsum], data,
			legend_def = legend,
			lower_y_range=0.5,
            method='datamc',
		    xtitle='%s #lambda_{M}' % cat_name,
		    ytitle='Events',
			writeTo=cat_name,
			)
		#set_trace()
	table.add_separator()
	with open('%s/yields.raw_txt' % out_dir,'w') as tab:
		tab.write('%s\n' % table)
		if pfit == 'postfit':
			#set_trace()
			pars = rootpy.asrootpy(mlfit_file.fit_s.floatParsFinal())
			tab.write(print_var_line('charmSF', pars))
			if 'lightSF' in pars:
				tab.write(print_var_line('lightSF', pars))

