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

wps = [
   'csvLoose',
   'csvMedium',
   'csvTight',
   'ctagLoose',
   'ctagMedium',
   'ctagTight',
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
table = open('plots/%s/ctageff/mass_discriminant/prefit_yields.tex' % jobid, 'w')
results = open('plots/%s/ctageff/mass_discriminant/results.tex' % jobid, 'w')
ncols = len(samples)+3
table.write(
   '''
\\begin{tabular}{%s}
\hline
''' % ('c'*ncols)
)
table.write(tex.tabline(['category', 'observed', 'total expected']+[i for _, i in samples], convert=False))
results.write(
   '''
\\begin{tabular}{cc}
'''
)


for wp in wps:
   table.write('\hline \n')
   table.write('\multicolumn{%d}{c}{ %s } \\\\ \n' % (ncols, wp))
   table.write('\hline \n')

   results.write('''\hline
\multicolumn{2}{c}{ %s } \\\\
\hline
''' % wp)
      
   data_f = io.root_open('plots/%s/ctageff/mass_discriminant/%s/datacard.root' % (jobid, wp))
   mc_f = io.root_open('plots/%s/ctageff/mass_discriminant/%s/MaxLikeFit.root' % (jobid, wp))
   prefit = mc_f.norm_prefit
   pshapes = mc_f.shapes_prefit
   pars = rootpy.asrootpy(mc_f.fit_s.floatParsFinal())
   results.write(
      tex.tabline(
         ['c--jet scale factor', tex.format_roovar(pars['charmSF'], True)], convert=False
         )
      )
   if 'lightSF' in pars:
      results.write(
         tex.tabline(
            ['light jet scale factor', tex.format_roovar(pars['lightSF'], True)], convert=False
            )
         )

   for cat in categories:
      line = [cat, '%d' % data_f.Get('%s/data_obs' % cat).Integral()]
      total = pshapes.Get('%s/total' % cat)
      t_err = ROOT.Double()
      t_val = total.IntegralAndError(1, total.nbins(), t_err)
      line.append(tex.format_with_error(t_val, t_err))
      for sam, _ in samples:
         path = '%s/%s' % (cat, sam)
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
