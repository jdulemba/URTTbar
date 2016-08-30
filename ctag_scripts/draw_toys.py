#! /bin/env python
# coding: utf-8
from rootpy.io import root_open
from URAnalysis.Utilities.roottools import Envelope
from pdb import set_trace
from argparse import ArgumentParser
from glob import glob
from rootpy.plotting import Hist, Hist2D
import re
import os 
from rootpy import asrootpy
from URAnalysis.PlotTools.BasePlotter import BasePlotter
from uncertainties import ufloat
from rootpy.plotting import Canvas, Graph
from rootpy import ROOT
import URAnalysis.Utilities.quad as quad
ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)

parser = ArgumentParser()
parser.add_argument('wp')
parser.add_argument('--asimov', '-a', action='store_true', help='run on asimov scan')
parser.add_argument('--conly', '-c', action='store_true', help='use only charmSF')
parser.add_argument('--usemean', '-m', action='store_true', help='use mean instead of gaussian fit')
parser.add_argument('--dir', '-d', default='toys_simple', help='directory name')
parser.add_argument('--precomputed', '-p', action='store_true', help='did I already got through the toys')
parser.add_argument('--systematic', '-s', type=float, default=0., help='Test the effect of an additional systematic')

args = parser.parse_args()

toydir= args.dir
basedir = '%s/plots/%s/ctageff/mass_discriminant/%s' % (os.environ['URA_PROJECT'], os.environ['jobid'], args.wp)
indir=os.path.join(basedir, toydir)
toy_regex=re.compile('toys_res\.\d+\.(?P<csf>[01]\.\d)\.(?P<lsf>[01]\.\d)\.root$') if not args.conly else \
	 re.compile('toys_res\.\d+\.(?P<csf>[01]\.\d)\.root$')

def get_toy_point(fname):
	base=os.path.basename(fname)
	m = toy_regex.match(base)
	if not m:
		raise RuntimeError('problem matching %s' % base)
	return float(m.group('csf')), float(m.group('lsf')) if not args.conly else None

def covers(realvar, val, sys=None):
	delta = realvar.value - val   
	err_idx = 1 if delta >= 0 else 0
	delta = abs(delta)
	if hasattr(realvar.error, '__getitem__'): #two sided error        
		err = abs(realvar.error[err_idx])
	else: #one sided error
		err = realvar.error
	if sys is not None:
		err = quad.quad(err, realvar.value*sys)
	return delta < err
	
if not args.precomputed:
	toys_files = glob('%s/toys_res*.root' % indir)
	csf_points = set()
	lsf_points = set()
	points = {}
	for toy in toys_files:
		key = get_toy_point(toy)
		if key in points: points[key].append(toy)
		else: points[key] = [toy]
		csf, lsf = key
		csf_points.add(csf)
		if not args.conly:
			lsf_points.add(lsf)
	
	csf_points = sorted(list(csf_points))
	lsf_points = sorted(list(lsf_points))

	def make_binning(lst):
		lst = [(2*lst[0] - lst[1])] + lst
		lst.append(
			(2*lst[-1] - lst[-2])
			)
		ll = len(lst)
		ret = []
		for idx in range(ll-1):
			point = lst[idx]
			next  = lst[idx+1]
			ret.append((next+point)/2)
		return ret
		
	csf_points = make_binning(csf_points)
	lsf_points = make_binning(lsf_points) if not args.conly else []

	cbias   = Hist(csf_points) if args.conly else Hist2D(csf_points, lsf_points)
	lbias   = None if args.conly else Hist2D(csf_points, lsf_points)
	ccover  = Hist(csf_points) if args.conly else Hist2D(csf_points, lsf_points)
	lcover  = None if args.conly else Hist2D(csf_points, lsf_points)
	failing = Hist(csf_points) if args.conly else Hist2D(csf_points, lsf_points)
	if args.systematic:
		ccover_p = Hist(csf_points) if args.conly else Hist2D(csf_points, lsf_points)
		lcover_p = None if args.conly else Hist2D(csf_points, lsf_points)
		

	out = root_open('%s/toys_summary.root' % indir, 'w')
	hdir = out.mkdir('histograms')

	it=0
	mtot =0;
	for key, files in points.iteritems():
		#print key
		it +=1
		cval, lval = key
		csf = Hist(200, 0, 2) if not args.asimov else None
		lsf = Hist(200, 0, 2) if not args.asimov and not args.conly else None
		csf_prefit = Hist(200, 0, 2) if not args.asimov else None
		lsf_prefit = Hist(200, 0, 2) if not args.asimov and not args.conly else None
		csf2d = Hist2D(100, 0, 2, 100, 0, 2) if not args.asimov else None
		lsf2d = Hist2D(100, 0, 2, 100, 0, 2) if not args.asimov and not args.conly else None
		failed = 0
		tot = 0
		cover_tot = 0
		lcover_ok = 0
		ccover_ok = 0
		lcoverp_ok = 0
		ccoverp_ok = 0
		for toy_file in files:
			with root_open(toy_file) as mlfit:
				toys = [''] if args.asimov else [i.name for i in mlfit.keys() if i.name.startswith('toy_')]         
				tot += len(toys)
				for i in toys:            
					path = os.path.join(i, 'fit_s')
					toy = mlfit.Get(path)
					pars = asrootpy(toy.floatParsFinal())
					prefit = asrootpy(toy.floatParsInit())
					cover_tot += 1
					if args.asimov:
						csf = float(pars['charmSF'].value)
						if not args.conly:
							lsf = float(pars['lightSF'].value)
					else:
						csf.Fill(float(pars['charmSF'].value))
						if not args.conly:
							lsf.Fill(float(pars['lightSF'].value))
					if not args.conly and covers(pars['lightSF'], lval):
						lcover_ok += 1
					if covers(pars['charmSF'], cval):
						ccover_ok += 1
					if args.systematic:
						if not args.conly and covers(pars['lightSF'], lval, args.systematic):
							lcoverp_ok += 1
						if covers(pars['charmSF'], cval, args.systematic):
							ccoverp_ok += 1						
					if not args.asimov:
						csf_prefit.Fill(float(prefit['charmSF'].value))
						csf2d.Fill(prefit['charmSF'].value, pars['charmSF'].value)
						if not args.conly:
							lsf_prefit.Fill(float(prefit['lightSF'].value))
							lsf2d.Fill(prefit['lightSF'].value, pars['lightSF'].value)
					toy.IsA().Destructor(toy)


		if args.conly:
			ccover.Fill(cval, float(ccover_ok)/cover_tot)
			if args.systematic:
				ccover_p.Fill(cval, float(ccoverp_ok)/cover_tot)			
		else:
			ccover.Fill(cval, lval, float(ccover_ok)/cover_tot)
			lcover.Fill(cval, lval, float(lcover_ok)/cover_tot)
			if args.systematic:
				ccover_p.Fill(cval, lval, float(ccoverp_ok)/cover_tot)
				lcover_p.Fill(cval, lval, float(lcoverp_ok)/cover_tot)

		mtot = max(tot, mtot) if not args.asimov else 1
		fcn = BasePlotter.parse_formula(
			'amplitude*TMath::Gaus(x, mean, width)', 
			'amplitude[100, 0, 1000], mean[1, 0, 2], width[0.5, 0, 2]',
			[0,2]
			)
		if args.conly:
			failing.Fill(cval, failed)
		else:
			failing.Fill(cval, lval, failed)
		xb = cbias.xaxis.find_bin(cval)
		if not args.conly:
			yb = cbias.yaxis.find_bin(lval)
		
		if args.usemean:
			bias = ufloat(csf.GetMean(), 0)
		elif args.asimov:
			bias = ufloat(csf, 0)
		else:
			csf.Fit(fcn, 'LQN0')
			bias = ufloat(fcn['mean'].value, fcn['mean'].error)
		bias -= cval
		bias /= cval
		bin = cbias[xb] if args.conly else cbias[xb,yb]
		bin.value = bias.nominal_value
		bin.error = bias.std_dev
		
		if not args.conly:
			if args.usemean:
				bias = ufloat(lsf.GetMean(), 0)
			elif args.asimov:
				bias = ufloat(lsf, 0)
			else:
				lsf.Fit(fcn, 'LQN0')
				bias = ufloat(fcn['mean'].value, fcn['mean'].error)
			bias -= lval
			bias /= lval
			lbias[xb,yb].value = bias.nominal_value
			lbias[xb,yb].error = bias.std_dev		
				
		if not args.asimov:
			tag = '%.1f' % cval if args.conly else '%.1f_%.1f' % (cval, lval)
			hdir.WriteTObject(csf, 'csf_%s' % tag)
			hdir.WriteTObject(csf_prefit, 'csf_prefit_%s' % tag)
			hdir.WriteTObject(csf2d, 'csf2D_%s' % tag)
			if not args.conly:
				hdir.WriteTObject(lsf2d, 'lsf2D_%s' % tag)
				hdir.WriteTObject(lsf, 'lsf_%s' % tag)
				hdir.WriteTObject(lsf_prefit, 'lsf_prefit_%s' % tag)
					 
	for i in failing:
		new_val = i.value/mtot
		i.value = new_val
	out.WriteTObject(cbias.Clone(), 'cbias')
	out.WriteTObject(failing.Clone(), 'failing')
	out.WriteTObject(ccover.Clone(), 'ccover')
	if args.systematic:
		out.WriteTObject(ccover_p.Clone(), 'ccover_p')		
	if not args.conly:
		out.WriteTObject(lbias.Clone(), 'lbias')
		out.WriteTObject(lcover.Clone(), 'lcover')
		if args.systematic:
			out.WriteTObject(lcoverp.Clone(), 'lcoverp')		
else:
	infile = root_open('%s/toys_summary.root' % indir)
	cbias = infile.cbias
	failing = infile.failing
	ccover = infile.ccover
	if not args.conly:
		lbias = infile.lbias
		lcover = infile.lcover

#
# Draw plots
#
ROOT.gStyle.SetOptStat(0)

fitfile = root_open('%s/MaxLikeFit.root' % basedir)
fit = fitfile.fit_s
pars = asrootpy(fit.floatParsFinal())
fitted = None
if args.conly:
	fitted = ROOT.TArrow(pars['charmSF'].value, 2, pars['charmSF'].value, 1, 0.025, '|>')
	fitted.SetLineWidth(3)
	fitted.SetLineColor(ROOT.kBlue)
	fitted.SetFillColor(ROOT.kBlue)
else:
	fitted = Graph(1)
	fitted.SetPoint(0, pars['charmSF'].value, pars['lightSF'].value)
	fitted.markerstyle = 20
	fitted.markersize  = 3
	fitted.markercolor = '#009600'

def addlabels(h):
	 h.xaxis.title = 'charm SF'
	 h.yaxis.title = 'light SF'
	 h.xaxis.SetTitleOffset(0.9)
	 h.yaxis.SetTitleOffset(1.2 if args.conly else 1.)

addlabels(cbias)
addlabels(failing)
addlabels(ccover)
if args.conly:
	cbias.yaxis.title = 'bias'
	failing.yaxis.title = 'filure rate'
	ccover.yaxis.title = 'coverage'
	if args.systematic:
		addlabels(ccover_p)
		ccover_p.yaxis.title = 'coverage'
else:
	addlabels(lbias)
	addlabels(lcover)
	if args.systematic:
		addlabels(ccover_p)
		addlabels(lcover_p)

mask = failing.Clone()
mask.Reset()

#for lb, cb, mb, fb in zip(lbias, cbias, mask, failing):
#	if fb.value > 0.2:
#		lb.value = -100.
#		cb.value = -100.
#		mb.value = 1.
#	else:
#		mb.value = -1.

ROOT.gStyle.SetPaintTextFormat('.2f')
axis = 'yaxis' if args.conly else 'zaxis'
yrange =  (-0.2, 0.2) if args.conly else (-1, 1)
getattr(cbias, axis).range_user = yrange
if not args.conly:
	lbias.zaxis.range_user = (-1, 1)
canvas = Canvas(800, 800)
canvas.margin = (0.15, 0.15, 0.15, 0.1)
from array import array
ROOT.TColor.CreateGradientColorTable(
	 5, array('d', [0., 0.3, .5, 0.7, 1.]), 
	 array('d', [0., 0., .95, 1., 1.]), 
	 array('d', [0., 0., .95, 0., 0.]), 
	 array('d', [1., 1., .95, 0., 0.]), 255
	 )
ROOT.gStyle.SetNumberContours(255);

#mask.fillstyle = 3354
#mask.fillcolor = 1
#mask.zaxis.range_user = (-0.9, 0.9)

hdrawopts = 'HIST text' if args.conly else 'colz text'
fitdrawopts = '' if args.conly else 'P SAME'

if not args.conly:
	lbias.Draw(hdrawopts)
	fitted.Draw(fitdrawopts)
	canvas.SaveAs('%s/lbias.png' % indir)
	canvas.SaveAs('%s/lbias.pdf' % indir)

cbias.Draw(hdrawopts)
if args.conly:
	fitted.SetY1(-0.13)
	fitted.SetY2(-0.19)	
fitted.Draw(fitdrawopts)
canvas.SaveAs('%s/cbias.png' % indir)
canvas.SaveAs('%s/cbias.pdf' % indir)

ROOT.TColor.CreateGradientColorTable(
	2, array('d', [0., 1.]), 
	array('d', [0., 1.]), 
	array('d', [0., 0.]), 
	array('d', [1., 0.]), 256
	)
failing.Draw(hdrawopts) #'HIST' if args.conly else 'colz')
getattr(failing, axis).range_user = (-0.1, 1)
if args.conly:
	fitted.SetY1(0.075)
	fitted.SetY2(-0.075)
fitted.Draw(fitdrawopts) #'' if args.conly else 'P SAME')
canvas.SaveAs('%s/fail.png' % indir)
canvas.SaveAs('%s/fail.pdf' % indir)

ROOT.TColor.CreateGradientColorTable(
	3, array('d', [0., 0.682, 1.]), 
	array('d', [0., .95, 1.]), 
	array('d', [0., .95, 0.]), 
	array('d', [1., .95, 0.]), 255
	)
ccover.Draw(hdrawopts)#'HIST' if args.conly else 'colz text')
getattr(ccover, axis).range_user = (0., 1)
add = None
if args.conly:
	fitted.SetY1(0.175)
	fitted.SetY2(0.025)
	add = ROOT.TLine(ccover[1].x.low, 0.682, ccover[-1].x.low, 0.682)
	add.SetLineWidth(2)
	add.SetLineColor(ROOT.kRed)
	add.SetLineStyle(2)
	add.Draw()
fitted.Draw(fitdrawopts)#'' if args.conly else 'P SAME')
canvas.SaveAs('%s/ccover.png' % indir)
canvas.SaveAs('%s/ccover.pdf' % indir)

if args.systematic:
	ccover_p.Draw(hdrawopts)#'HIST' if args.conly else 'colz text')
	getattr(ccover_p, axis).range_user = (0., 1)
	if add is not None:
		add.Draw()
	fitted.Draw(fitdrawopts)
	canvas.SaveAs('%s/ccover_p.png' % indir)	
	canvas.SaveAs('%s/ccover_p.pdf' % indir)	

if not args.conly:
	lcover.Draw(hdrawopts)
	lcover.zaxis.range_user = (0., 1)
	fitted.Draw(fitdrawopts)
	canvas.SaveAs('%s/lcover.png' % indir)
	canvas.SaveAs('%s/lcover.pdf' % indir)

	if args.systematic:
		lcover_p.Draw(hdrawopts)#'HIST' if args.conly else 'colz text')
		getattr(lcover_p, axis).range_user = (0., 1)
		if add is not None:
			add.Draw()
		fitted.Draw(fitdrawopts)
		canvas.SaveAs('%s/lcover_p.png' % indir)	
		canvas.SaveAs('%s/lcover_p.pdf' % indir)	
