#! /bin/env python
'''
Top Spin Gen Plotting macro
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
from rootpy.plotting import views, Graph, Hist
from argparse import ArgumentParser
from URAnalysis.Utilities.roottools import slice_hist
import ROOT

project = os.environ['URA_PROJECT']
parser = ArgumentParser()
parser.add_argument('todo', nargs='+', help='things to do', choices=['shapes', 'comparisons', 'effs'])
args = parser.parse_args()

jobid = os.environ['jobid']
plotter = BasePlotter(
   '%s/plots/%s/tsgen' % (project, jobid),
   defaults = {'save' : {'png' : True, 'pdf' : False}}
)

styles = {
   'lep' : {
      'title' : 'l',
      'linecolor' : 'blue',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'nu' : {
      'title' : '#nu',
      'linecolor' : 'cyan',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'utype' : {
      'title' : 'up-type',
      'linecolor' : 'red',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'dtype' : {
      'title' : 'down-type',
      'linecolor' : 'violet',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   '' : {
      'title' : 'right perm',
      'linecolor' : 'blue',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'wrong' : {
      'title' : 'wrong perm',
      'linecolor' : 'red',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'wrongdtype' : {
      'title' : 'wrong (d-type)',
      'linecolor' : 'black',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'wrongutype' : {
      'title' : 'wrong (u-type)',
      'linecolor' : '#00960a',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      }
}

tfile = root_open('%s/results/%s/topspin_gen/ttJets.root' % (project, jobid))

variables = [
   (16, 'NE', 'helframe_costheta_%s', 'cos #theta_{hel}', ['dtype', 'utype', 'lep', 'nu']),
   (8, 'SE', 'helframe_cosdelta_lep_%s'   , 'cos #delta(l, q)', ['dtype', 'utype']),
   (4, 'NE', 'helframe_prodcosth_lep_%s'  , 'cos #theta_{hel}(l) cos #theta_{hel}(q)', ['dtype', 'utype']),
   (4, 'NE', 'labframe_cosdeltaphi_lep_%s', 'cos #Delta#varphi(l, q)', ['dtype', 'utype']),
	 (4, 'NE', 'wframe_costheta_%s'  , 'cos #theta*, W frame', ['dtype', 'utype']), #, 'lep', 'nu']),
	 (4, 'NE', 'wframe_energy_%s'    , 'energy, W frame', ['dtype', 'utype']),
	 (4, 'NE', 'tframe_energy_%s'    , 'energy, t frame', ['dtype', 'utype']),
	 (4, 'NE', 'labframe_energy_%s'  , 'energy, lab frame', ['dtype', 'utype']),
	 (8, 'NE', 'tframe_wlep_costheta%s', '', ['']),
	 (8, 'NE', 'tframe_whad_costheta%s', '', ['']),
	 (8, 'NE', 'tframe_blep_costheta%s', '', ['']),
	 (8, 'NE', 'tframe_bhad_costheta%s', '', ['']),
]

if 'shapes' in args.todo:
	for tdir in ['matched', 'parton']:
		path_dir = tdir
		plotter.set_subdir(path_dir)
		for rbin, lpos, var, xaxis, types in variables:
			to_draw = []
			add = ['wrongutype', 'wrongdtype'] if tdir == 'matched' else []         
			if types == [''] and tdir == 'matched':
				add = ['wrong']
			base_path = 'semilep/%s/%s' % (path_dir, var)            
			for ptype in types+add:
				if ptype.startswith('wrong'):
					postfix = ptype.replace('wrong','')
					path = base_path % postfix
					path = path.replace(tdir, 'wrong')
				else:
					path = base_path % ptype
				histo = tfile.Get(path).Clone()
				if histo.Integral():
					histo.Scale(1./histo.Integral())
				else:
					print '%s has not integral!' % path
				histo.Rebin(rbin)
				plotter.set_histo_style(histo, **styles[ptype])
				to_draw.append(histo)
			legend = LegendDefinition(position=lpos)
			plotter.overlay(to_draw, legend_def=legend, xtitle=xaxis, ytitle='a.u.', y_range='shape')
			plotter.save(var.replace('_%s', '').replace('%s', ''))

def make_view(inview, title, color, line='solid'):
	return views.TitleView(
		views.StyleView( 
			inview,
			linecolor=color,
			drawstyle='hist',
			linewidth=3,
			legendstyle='l',
				 linestyle=line,
			),
		title
		)

shapes = [
	make_view(tfile, 'tt jets', 'blue'),
	make_view(
		root_open('results/%s/topspin_gen/HtoTTM400.root' % jobid),
		'H M400', '#0a6d3d'),
	make_view(
		root_open('results/%s/topspin_gen/AtoTTM400.root' % jobid),
		'A M400','#fc3936'),
	##make_view(
	##	root_open('results/%s/topspin_gen/HtoTTM400_Interf.root' % jobid),
	##	'H M400 (I)', '#0a6d3d', 'dashed'),
	##make_view(
	##	root_open('results/%s/topspin_gen/AtoTTM400_Interf.root' % jobid),
	##	'A M400 (I)', '#fc3936', 'dashed'),
	## make_view(
	##	 root_open('results/%s/topspin_gen/HtoTTM600.root' % jobid),
	##	 'A M600','#6d0a3a'),
	## make_view(
	##	 root_open('results/%s/topspin_gen/HtoTTM800.root' % jobid),
	##	 'A M800','#59cfdf'),
	]

masses = [i.Get('semilep/parton/ttm').Rebin(1) for i in shapes]
with root_open('mass_weights.root', 'w') as outf:
	for i in masses:
		outf.WriteTObject(i, i.GetTitle())
	weight = masses[1]
	weight.Scale(1./weight.Integral())
	masses[0].Scale(1./masses[0].Integral())
	weight.Divide(masses[0])
	outf.WriteTObject(weight, 'weight')

variables = [
	(16, 'NE', 'helframe_costheta_dtype', 'cos #theta_{hel}(d)'),
	(16, 'NE', 'helframe_costheta_utype', 'cos #theta_{hel}(u)'),
	(16, 'NE', 'helframe_costheta_lep'  , 'cos #theta_{hel}(l)'),
	(8, 'SE', 'helframe_cosdelta_lep_dtype', 'cos #delta(l, d)'),
	(8, 'SE', 'helframe_cosdelta_lep_utype', 'cos #delta(l, u)'),
	(4, 'NE', 'helframe_prodcosth_lep_dtype', 'cos #theta_{hel}(l) cos #theta_{hel}(d)'),
	(4, 'NE', 'helframe_prodcosth_lep_utype', 'cos #theta_{hel}(l) cos #theta_{hel}(u)'),
	(4, 'NE', 'labframe_cosdeltaphi_lep_dtype', 'cos #Delta#varphi(l, d)'),
	(4, 'NE', 'labframe_cosdeltaphi_lep_utype', 'cos #Delta#varphi(l, u)'),
	#Kinvars
	(range(0, 305, 15), 'NE',  "lpt",  'p_{T}(l)' ),
	(range(0, 305, 15), 'NE', "blpt",  'p_{T}(bl)'),
	(range(0, 305, 15), 'NE', "bhpt",  'p_{T}(bh)'),
	(range(0, 305, 15), 'NE', "wupt",  'p_{T}(wu)'),
	(range(0, 305, 15), 'NE', "wdpt",  'p_{T}(wd)'),
	(10, 'NE',  "leta", '|#eta(l)|' ),
	(10, 'NE', "bleta", '|#eta(bl)|'),
	(10, 'NE', "bheta", '|#eta(bh)|'),
	(10, 'NE', "wueta", '|#eta(wu)|'),
	(10, 'NE', "wdeta", '|#eta(wd)|'),      
	(range(300, 1010, 25), 'NE', "ttm",	'M(tt)' ),
	(range(0, 505, 30), 'NE', "ttpt",  'p_{T}(tt)'),
	(10, 'NE', "tteta", '|#eta(tt)|'),
	 (4, 'NE', 'wframe_costheta_dtype'  , 'cos #theta*, W frame'),
	 (4, 'NE', 'wframe_costheta_utype'  , 'cos #theta*, W frame'),
	 (4, 'NE', 'wframe_costheta_lep'	 , 'cos #theta*, W frame'),

	 (4, 'NE', 'wframe_energy_dtype', 	'Energy(d), W frame'),
	 (4, 'NE', 'wframe_energy_utype', 	'Energy(u), W frame'),
	 (4, 'NE', 'tframe_energy_dtype', 	'Energy(d), t frame'),
	 (4, 'NE', 'tframe_energy_utype', 	'Energy(u), t frame'),
	 (4, 'NE', 'labframe_energy_dtype', 'Energy(d), lab frame'),
	 (4, 'NE', 'labframe_energy_utype', 'Energy(u), lab frame'),

	 (16, 'NE', "tframe_wlep_costheta", 'cos #theta*(W_{lep})'),
	 (16, 'NE', "tframe_whad_costheta", 'cos #theta*(W_{had})'),
	 (16, 'NE', "tframe_blep_costheta", 'cos #theta*(b_{lep})'),
	 (16, 'NE', "tframe_bhad_costheta", 'cos #theta*(b_{had})'),
	 (8, 'NE', "ttframe_top_costheta", 'cos #theta*(t)'),
	 (8, 'NE', "ttframe_tba_costheta", 'cos #theta*(#bar{t})'),
	]

if 'comparisons' in args.todo:
	for tdir in ['parton', 'matched']:
		path_dir = tdir
		plotter.set_subdir('comparison_%s' % path_dir)
		for rbin, lpos, var, xaxis in variables:
			to_draw = []         
			path = ('semilep/%s/%s' % (path_dir, var))
			for vv in shapes:
				histo = vv.Get(path)
				if histo.get_dimension() == 2:
					prx = histo.ProjectionX()
					if prx is None:
						set_trace()
					prx = asrootpy(prx)
					prx.decorate(**histo.decorators)
					histo = prx.Clone()
				if histo.Integral():
					histo.Scale(1./histo.Integral())
				else:
					print '%s has not integral!' % path
				histo = urviews.RebinView.rebin(histo, rbin)
				to_draw.append(histo)
	
			legend = LegendDefinition(position=lpos)
			plotter.overlay(to_draw, legend_def=legend, xtitle=xaxis, ytitle='a.u.', y_range='shape', ignore_style=True)
			plotter.save(var)


def make_eff(num, den, xtitle='', ytitle='', **decorators):
	ret = asrootpy(ROOT.TGraphAsymmErrors(num, den, 'cp'))
	ret.xaxis.title = xtitle
	ret.yaxis.title = ytitle
	ret.decorate(**decorators)
	return ret

def make_slices(hist, nslices):
	nbins = hist.GetNbinsY()
	rebin = nbins/nslices
	if nbins % nslices != 0:
		raise RuntimeError("I cannot create %d slices" % nslices)
	newh = urviews.RebinView.rebin(hist, [[1], [rebin]])
	ret = []
	for ibin in range(1, nslices+1):
		m = newh.yaxis.GetBinLowEdge(ibin)
		M = newh.yaxis.GetBinUpEdge(ibin)
		hsl = slice_hist(newh, ibin)
		hsl.title = '%.0f < m(tt) < %.0f' % (m, M)
		ret.append(hsl)
	return ret

def plot_res(plotter, hist, nslices, nicevar, linestyle='solid', legend=True):
	colors = [
		'#0000ff',
		'#ff0000',
		'#008000',
		'#660066',
		'#ff9b02',
		]
	slices = make_slices(hist, nslices)
	for i in slices:
		if i.Integral():
			i.Scale(1./i.Integral())
	gen = nicevar % 'gen'
	reco = nicevar % 'reco'
	xtitle = '(%s - %s)/%s' % (reco, gen, gen)
	plotter.overlay(
		slices,
		LegendDefinition(position='NE') if legend else None,
		xtitle=xtitle,
		ytitle='A.U.',
		linestyle=linestyle,
		linecolor=colors[:len(slices)],
		fillstyle='hollow',
		drawstyle='hist' if legend else 'hist same',
		legendstyle='l',
		linewidth=2,
		markerstyle=0
		)

def single_bins(plotter, right, wrong, nicevar, nslices, fitpars, relative=True, skip=0):
	right_slices = make_slices(right, nslices)
	wrong_slices = make_slices(wrong, nslices)	
	ibin = 0
	resolution = []
	bias = []
	binning = set()
	for idx, info in enumerate(zip(right_slices, wrong_slices)):
		right, wrong = info
		if idx < skip: continue
		if right.Integral():
			right.Scale(1./right.Integral())
		else:
			continue
		if wrong.Integral():
			wrong.Scale(1./wrong.Integral())
		
		gen = nicevar % 'gen'
		reco = nicevar % 'reco'
		xtitle = '(%s - %s)/%s' % (reco, gen, gen)
		print '###############  Bin%d  #####################' % ibin
		x_range = [-2, 2] if relative else [-1000, 1000]
		tf1 = plotter.parse_formula(
			fitpars[0], fitpars[1],
			x_range
			)
		right.Fit(tf1, 'WL')		
		tf1.linewidth = 2
		tf1.linecolor = 'blue'

		rr = right.title.split('<')
		binning.add(float(rr[-1]))
		binning.add(float(rr[0]))
		resolution.append((
				tf1.GetParameter('width'), tf1.GetParError(tf1.GetParNumber('width'))
				))
		bias.append((
				tf1.GetParameter('mean'), tf1.GetParError(tf1.GetParNumber('mean'))
				))

		plotter.overlay(
			[right, wrong, tf1],
			None,
			xtitle=xtitle,
			ytitle='A.U.',
			linestyle=['solid', 'dashed'],
			fillstyle='hollow',
			drawstyle='hist',
			legendstyle='l',
			linewidth=2,
			markerstyle=0
			)
		plotter.save('Bin%d' % ibin)
		ibin += 1
	res = Hist(sorted(list(binning)))
	for idx, info in enumerate(resolution):
		res[idx+1].value, res[idx+1].error = info
		if not relative:
			res[idx+1].value /= res[idx+1].x.center
			res[idx+1].error /= res[idx+1].x.center
	
	plotter.plot(res, xtitle='m_{gen}(tt)', ytitle='resolution')
	plotter.save('resolution')

	res = Hist(sorted(list(binning)))
	for idx, info in enumerate(bias):
		res[idx+1].value, res[idx+1].error = info
		if not relative:
			res[idx+1].value /= res[idx+1].x.center
			res[idx+1].error /= res[idx+1].x.center
	
	plotter.plot(res, xtitle='m_{gen}(tt)', ytitle='bias')
	plotter.save('bias')


dgauss = 'amplitude*TMath::Gaus(x, mean, width)+a2*TMath::Gaus(x, mean, scale*width)', 'amplitude[0.1, 0, 1], mean[0, -1, 1], width[0.5, 0, 1], a2[0.1, 0, 1], scale[3, 1, 1000]'
gauss = 'amplitude*TMath::Gaus(x, mean, width)', 'amplitude[0.1, 0, 1], mean[0, -1, 1], width[0.5, 0, 1]'
bwinger = 'amplitude*TMath::BreitWigner(x, mean, width)', 'amplitude[0.1, 0, 1], mean[0, -1, 1], width[0.5, 0, 1]'
dgauss_abs = 'amplitude*TMath::Gaus(x, mean, width)+a2*TMath::Gaus(x, mean, scale*width)', 'amplitude[0.1, 0, 1], mean[0, -1000, 1000], width[50, 0, 800], a2[0.1, 0, 1], scale[3, 1, 1000]'
bwgauss = 'amplitude*TMath::BreitWigner(x, mean, width)+a2*TMath::Gaus(x, mean, scale*width)', 'amplitude[0.1, 0, 1], mean[0, -1, 1], width[0.5, 0, 1], a2[0.1, 0, 1], scale[3, 1, 1000]'
if 'effs' in args.todo:
	## plotter.set_subdir('efficiencies')
	## denominator = tfile.semilep.parton.eff_mtt.Clone()
	## for pic, name, dirname in [
	##	 ('obj_eff', 'Obj. selection | tt semi-lep', 'obj_selection'),
	##	 ('match_eff', 'MC matching | obj. selection', 'matched'),
	##	 ('mass_eff', 't_{had} mass | MC matching', 'mass_selection'),
	##	 ('perm_eff', 'Perm. selection | t_{had} mass', 'perm_selection'),
	##	 ('sel_eff', 'Correct perm. | Perm. selection', 'selected')
	##	 ]:
	##	 numerator = tfile.semilep.Get('%s/eff_mtt' % dirname)
	##	 obj_sel_eff = make_eff(numerator, denominator, xtitle='gen m(tt) [GeV]', ytitle='#varepsilon(%s)' % name)
	##	 obj_sel_eff.drawstyle = 'AP'
	##	 plotter.plot(obj_sel_eff)
	##	 plotter.save(pic)
	##	 denominator = numerator
	
	for method, rname, wname in [('norm', 'matched', 'wrong'), ('kinfit', 'matchfit', 'wrongfit')]:
		for nslices, fitpars, resname, xtit in [
			(16, dgauss, "res_mtt", "m^{tt}(%s)"), 
			#(16, dgauss, "res_mth", "m^{th}(%s)"),
			#(16, bwinger, 'res_deltatt', 'cos(#delta_{tt})(%s)'),
			#(16, dgauss, "res_ptl", 'p^{tlep}(%s)'), 
			#(16, dgauss, "res_pth", 'p^{thad}(%s)')
			#(10, dgauss, "res_jete", "E^{jet}(%s)"),
			#(10, gauss, 'res_deltaw_ptw', '%s'),
			#(10, gauss, 'res_mw_ptw', '%s'),
			#(8, gauss, "res_mwh", "m^{Wh}(%s)"),
			#(8, gauss, "res_pttt", "p_{T}^{tt}(%s)"), 
			]:
			plotter.set_subdir('efficiencies/%s' % method)
			right = tfile.Get('semilep/%s/%s' % (rname, resname))
			wrong = tfile.Get('semilep/%s/%s' % (wname, resname))
			plot_res(plotter, right, 2, xtit, 'solid', True)
			plot_res(plotter, wrong, 2, xtit, 'dashed', False)
			plotter.save(resname)
			plotter.set_subdir('efficiencies/%s/%s' % (method, resname))
			single_bins(
				plotter, right, wrong, xtit, 
				nslices, fitpars, not resname.endswith('_abs'),
				2
				)
		
