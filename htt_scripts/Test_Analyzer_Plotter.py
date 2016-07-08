'''
Test Analyzer Plotter macro
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
from ROOT import TFile

plotter = BasePlotter(
   'plots/test_analyzer',
   defaults = {'save' : {'png' : True, 'pdf' : False}}
)

myfile = root_open('../results/2016May30/Test_Analyzer/ttJets.root', 'read') #ttJets full analysis file
#myfile = root_open('../ttJets.Test_Analyzer.test.root', 'read') #ttJets sample file

variables = [
	('helframe_costheta_%s', 'cos #theta*', 'helframe_costheta'),
	('labframe_cosdeltaphi_lep_%s', 'cos #Delta#varphi(l, q)', 'labframe_cosdeltaphi_lep'),
	('labframe_deltaphi_lep_%s', '#Delta#varphi(l, q)', 'labframe_deltaphi_lep'),
	('helframe_cosdelta_lep_%s', 'cos #delta(l,q)', 'helframe_cosdelta_lep'),
	('helframe_prodcosth_lep_%s', 'cos #theta*_{l} cos #theta*_{q}', 'helframe_prodcosth_lep'),
]

defcol = 'black' # default color
yaxis = 'a.u.' # yaxis title
lpos = 'NE' #legend position
types = ['dtype', 'utype'] #variable is down- or up-type jet
colors = ['red', 'blue'] #plot legend colors (down is red, up is blue)
labels = ['down', 'up'] #legend variable labels

for var, xaxis, title in variables: #variable name in histogram, xaxis title, saved title of png
	j = 0
	to_draw = []
	for i in types:
		hist = asrootpy(myfile.Get(var % i)).Clone()
		plotter.set_histo_style(hist, title=labels[j], color=colors[j])
		to_draw.append(hist)
		j = j+1
	plotter.overlay(to_draw, legend_def=LegendDefinition(position=lpos), xtitle=xaxis, ytitle=yaxis)
	plotter.save(title)

#delta plots
deltvars = [
	('delta_helframe_costheta_hist', 'Helframe cos #theta*', 'blue'),
        ('delta_labframe_cosdeltaphi_hist', 'Labframe cos #Delta#varphi(l, q)', 'cyan'),
        ('delta_labframe_deltaphi_hist', 'Labframe #Delta#varphi(l, q)', 'red'),
        ('delta_helframe_cosdelta_hist', 'Helframe cos #delta(l,q)', 'violet'),
        ('delta_helframe_prodcosth_hist', 'Helframe cos #theta*_{l} cos #theta*_{q}', 'black'),
	('frac_delta_E_hist', 'Frac. #Delta E', '#00960a')
]

to_draw = []

for var, labels, colors in deltvars: #variable name in histogram, legend label, variable color
	hist = asrootpy(myfile.Get(var)).Clone()
	plotter.set_histo_style(hist, title=labels, color=colors)
	to_draw.append(hist)

plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle='#Delta', ytitle=yaxis)
plotter.save('Delta')

#fractional delta plots
frac_deltvars = [
        ('frac_delta_helframe_costheta_hist', 'Helframe cos #theta*', 'blue'),
        ('frac_delta_labframe_cosdeltaphi_hist', 'Labframe cos #Delta#varphi(l, q)', 'cyan'),
        ('frac_delta_labframe_deltaphi_hist', 'Labframe #Delta#varphi(l, q)', 'red'),
        ('frac_delta_helframe_cosdelta_hist', 'Helframe cos #delta(l,q)', 'violet'),
        ('frac_delta_helframe_prodcosth_hist', 'Helframe cos #theta*_{l} cos #theta*_{q}', 'black'),
]

to_draw = []

for var, labels, colors in frac_deltvars: #variable name in histogram, legend label, variable color
        hist = asrootpy(myfile.Get(var)).Clone()
        plotter.set_histo_style(hist, title=labels, color=colors)
        to_draw.append(hist)

plotter.overlay(to_draw, legend_def=LegendDefinition(position='NE'), xtitle='Fractional #Delta', ytitle=yaxis)
plotter.save('Frac_Delta')


#Other variables
lep_vars = [
	('lep_%s', ['charge', 'pt'], ['lepton charge', 'lepton pt'])
]
for var, types, xaxis in lep_vars: #variable name in histogram, type of lepton variable, xaxis title
	j = 0
	for i in types:
		hist = asrootpy(myfile.Get(var % i)).Clone()
		plotter.set_histo_style(hist, color=defcol)
		plotter.plot(hist, xtitle=xaxis[j], ytitle=yaxis)
		plotter.save(var % i)
		j = j+1

ang_vars = [ #other angular variables
	('helframe_costheta_%s', ['lep', 'nu'], 'cos #theta*')
]

for var, types, xaxis in ang_vars: #variable name, type of ang. variable, xaxis title
	for i in types:
		hist = asrootpy(myfile.Get(var % i)).Clone()
		plotter.set_histo_style(hist, color=defcol)
                plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
                plotter.save(var % i)

#variables that have no pairs
SA_vars = [
	('delta_E_hist', '#Delta E', 'Delta_E'),
	('frac_delta_E_hist', 'Frac. #Delta E', 'Frac_Delta_E'),
	('nJets', 'nJets', 'nJets'),
	('cut_flow', 'Cut Flow', 'Cut_Flow')
]

for var, xaxis, name in SA_vars: #variable name, xaxis title, png saved name
	hist = asrootpy(myfile.Get(var)).Clone()
	plotter.set_histo_style(hist, color=defcol)
	plotter.plot(hist, xtitle=xaxis, ytitle=yaxis)
	plotter.save(name)
