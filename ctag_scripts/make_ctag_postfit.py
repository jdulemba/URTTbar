import rootpy.io as io
import rootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)
rootpy.log["/data_views"].setLevel(rootpy.log.WARNING)
log = rootpy.log["/make_postfit_plots"]
log.setLevel(rootpy.log.INFO)
from rootpy.plotting import Hist
from rootpy import asrootpy
from argparse import ArgumentParser
import os
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

def foo():
	pass

parser = ArgumentParser()
parser.add_argument('working_point')
parser.add_argument('plots', choices=['prefit', 'postfit', 'both'],
                    help='what to plot')
args = parser.parse_args()
plots_to_do = []
if args.plots == 'prefit' or args.plots == 'both':
	plots_to_do.append(('prefit', 'shapes_prefit'))
if args.plots == 'postfit' or args.plots == 'both':
	plots_to_do.append(('postfit', 'shapes_fit_s'))

def fix_binning(postfit, correct_bin):
	#set_trace()
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

input_dir = 'plots/%s/ctageff/mass_discriminant/%s' % (os.environ['jobid'], args.working_point)

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

plotter = BasePlotter(
	'%s/postfit' % input_dir,
	defaults = {'save' : {'png' : True, 'pdf' : False}},
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

def print_var_line(varname, tree):
        ## Get branch entries
	tree.GetBranch(varname).GetEntry(0)
	tree.GetBranch('%sLoErr' % varname).GetEntry(0)
	tree.GetBranch('%sHiErr' % varname).GetEntry(0)

	var  = tree.GetLeaf(varname).GetValue()
	MinErr = tree.GetLeaf('%sLoErr' % varname).GetValue()
	MaxErr = tree.GetLeaf('%sHiErr' % varname).GetValue()
	return 'best fit value %s: %.3f -%.3f/%.3f\n' % (varname, var, MinErr, MaxErr)

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
			if i.name.startswith('total') or 'data' in i.name: continue
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
		available_samples = [i.name for i in cat_dir.keys() if not i.name.startswith('total') and 'data' not in i.name]

		line = [cat_name, '%.1f' % data.Integral()]+['%.1f' % (cat_dir.Get(i).Integral() if i in available_samples else 0.0) for i in sample_names]
		table.add_line(*line)

		#set_trace()
		samples = [fix_binning(cat_dir.Get(i), data) for i in available_samples]
		samples.sort(key=ordering)

		stack = plotter.create_stack(*samples, sort=False)
		legend = LegendDefinition(position='NE')
		plotter.overlay_and_compare( 
			[stack, hsum], data,
			writeTo=cat_name,
			legend_def = legend,
		    xtitle='%s #lambda_{M}' % cat_name,
		    ytitle='Events',
		    method='datamc'
		)
	table.add_separator()
	with open('%s/yields.raw_txt' % out_dir,'w') as tab:
		tab.write('%s\n' % table)
		if pfit == 'postfit':
			tree = mlfit_file.tree_fit_sb
			#set_trace()
			tab.write(print_var_line('charmSF', tree))
			lightSF = [ i.GetName() for i in tree.GetListOfLeaves() if 'lightSF' in i.GetName()]
			if lightSF:
				tab.write(print_var_line('lightSF', tree))

