'''
Top Spin Gen Plotting macro
'''

from URAnalysis.PlotTools.BasePlotter import BasePlotter, LegendDefinition
import URAnalysis.PlotTools.views as urviews
import os
from rootpy.io import root_open
from rootpy import asrootpy
from pdb import set_trace
from rootpy.plotting import views
from argparse import ArgumentParser
from URAnalysis.Utilities.roottools import slice_hist
import ROOT

parser = ArgumentParser()
parser.add_argument('todo', nargs='+', help='things to do')
args = parser.parse_args()

allowed_to_do = set(['shapes', 'comparisons', 'effs'])
for i in args.todo:
   if i not in allowed_to_do:
      raise RuntimeError('%s is not an allowed command! You have: %s' % (i, ', '.join(allowed_to_do)))

jobid = os.environ['jobid']
plotter = BasePlotter(
   'plots/%s/tsgen' % jobid,
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
   'wrongd' : {
      'title' : 'wrong (d-type)',
      'linecolor' : 'black',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      },
   'wrongu' : {
      'title' : 'wrong (u-type)',
      'linecolor' : '#00960a',
      'drawstyle' : 'hist',
      'legendstyle' : 'l',
      }
}

tfile = root_open('results/%s/topspin_gen/ttJets.root' % jobid)

variables = [
   (16, 'NE', 'helframe_costheta_%s', 'cos #theta*', ['dtype', 'utype', 'lep', 'nu']),
   (8, 'SE', 'helframe_cosdelta_lep_%s'   , 'cos #delta(l, q)', ['dtype', 'utype']),
   (4, 'NE', 'helframe_prodcosth_lep_%s'  , 'cos #theta*_{l} cos #theta*_{q}', ['dtype', 'utype']),
   (4, 'NE', 'labframe_cosdeltaphi_lep_%s', 'cos #Delta#varphi(l, q)', ['dtype', 'utype']),
]

if 'shapes' in args.todo:
   for tdir in ['matched', 'part_acceptance', 'parton']:
      path_dir = tdir
      plotter.set_subdir(path_dir)
      for rbin, lpos, var, xaxis, types in variables:
         to_draw = []
         add = ['wrongu', 'wrongd'] if tdir == 'matched' else []         
         base_path = '%s/%s' % (path_dir, var)            
         for ptype in types+add:
            if ptype.startswith('wrong'):
               path = base_path % ('dtype' if ptype.endswith('d') else 'utype')
               path = path.replace(tdir, 'wrong')
            else:
               path = base_path % ptype
            histo = asrootpy(tfile.Get(path).ProjectionX()).Clone()
            if histo.Integral():
               histo.Scale(1./histo.Integral())
            else:
               print '%s has not integral!' % path
            histo.Rebin(rbin)
            plotter.set_histo_style(histo, **styles[ptype])
            to_draw.append(histo)
         legend = LegendDefinition(position=lpos)
         plotter.overlay(to_draw, legend_def=legend, xtitle=xaxis, ytitle='a.u.', y_range='shape')
         plotter.save(var.replace('_%s', ''))

def make_view(inview, title, color):
   return views.TitleView(
      views.StyleView( 
         inview,
         linecolor=color,
         drawstyle='hist',
         linewidth=3,
         legendstyle='l'
         ),
      title
      )

shapes = [
   make_view(tfile, 'tt jets', 'blue'),
   make_view(
      root_open('results/%s/topspin_gen/HtoTTM500_76X.root' % jobid), 
      'H M500', '#0a6d3d'),
   ## make_view(
   ##    root_open('results/%s/topspin_gen/HtoTTM400.root' % jobid),
   ##    'H M400', '#0a6d3d'),
   ## make_view(
   ##    root_open('results/%s/topspin_gen/AtoTTM400.root' % jobid),
   ##    'A M400','#fc3936'),
   ## make_view(
   ##    root_open('results/%s/topspin_gen/HtoTTM600.root' % jobid),
   ##    'A M600','#6d0a3a'),
   ## make_view(
   ##    root_open('results/%s/topspin_gen/HtoTTM800.root' % jobid),
   ##    'A M800','#59cfdf'),
   ]

variables = [
   (16, 'NE', 'helframe_costheta_dtype', 'cos #theta*(d)'),
   (16, 'NE', 'helframe_costheta_utype', 'cos #theta*(u)'),
   (16, 'NE', 'helframe_costheta_lep'  , 'cos #theta*(l)'),
   (8, 'SE', 'helframe_cosdelta_lep_dtype', 'cos #delta(l, d)'),
   (8, 'SE', 'helframe_cosdelta_lep_utype', 'cos #delta(l, u)'),
   (4, 'NE', 'helframe_prodcosth_lep_dtype', 'cos #theta*_{l} cos #theta*_{d}'),
   (4, 'NE', 'helframe_prodcosth_lep_utype', 'cos #theta*_{l} cos #theta*_{u}'),
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
   (range(300, 1010, 50), 'NE', "ttm",   'M(tt)' ),
   (range(0, 505, 30), 'NE', "ttpt",  'p_{T}(tt)'),
   (10, 'NE', "tteta", '|#eta(tt)|'),
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

if 'effs' in args.todo:
   plotter.set_subdir('efficiencies')
   denominator = tfile.semilep.parton.eff_mtt.Clone()
   for pic, name, dirname in [
      ('obj_eff', 'Obj. selection | tt semi-lep', 'obj_selection'),
      ('match_eff', 'MC matching | obj. selection', 'matched'),
      ('perm_eff', 'Perm. selection | MC matching', 'perm_selection'),
      ('sel_eff', 'Correct perm. | Perm. selection', 'selected')
      ]:
      numerator = tfile.semilep.Get('%s/eff_mtt' % dirname)
      obj_sel_eff = make_eff(numerator, denominator, xtitle='gen m(tt) [GeV]', ytitle='#varepsilon(%s)' % name)
      obj_sel_eff.drawstyle = 'AP'
      plotter.plot(obj_sel_eff)
      plotter.save(pic)
      denominator = numerator
   
   for resname, xtit in [
      ("res_mtt", "m^{tt}(%s)"), 
      ("res_pttt", "p_{T}^{tt}(%s)"), 
      ("res_ptl", 'p^{tlep}(%s)'), 
      ("res_pth", 'p^{thad}(%s)')]:
      right = tfile.Get('semilep/selected/%s' % resname)
      wrong = tfile.Get('semilep/wrong/%s' % resname)
      plot_res(plotter, right, 2, xtit, 'solid', True)
      plot_res(plotter, wrong, 2, xtit, 'dashed', False)
      plotter.save(resname)
