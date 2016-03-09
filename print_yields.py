import URAnalysis.Utilities.prettyjson as prettyjson
import os
from uncertainties import ufloat
from pdb import set_trace
import math
from URAnalysis.Utilities.quad import quad
import URAnalysis.Utilities.tables as tables

json_ufloat = lambda d, n : ufloat(d[n], d['%s_err' % n])

jobid = os.environ['jobid']
summary = prettyjson.loads(
   open('plots/%s/ctageff/summary.json' % jobid).read()
)
summary = summary["mass_discriminant"]
wpoints = summary['working_points']

mc_effs = prettyjson.loads(
   open('plots/%s/ctageff/mc_effs.json' % jobid).read()
)

#set_trace()
nevts  = json_ufloat(summary, 'nevts'    )  #ufloat(summary['nevts'    ], summary['nevts_err'])
lcfrac = json_ufloat(summary, 'leadCFrac')  #ufloat(summary['leadCFrac'], summary['leadCFrac_err'])
scfrac = json_ufloat(summary, 'subCFrac' )  #ufloat(summary['subCFrac' ], summary['subCFrac_err'])
llfrac = 1 - lcfrac - scfrac #all-light fraction

factor = lambda leff, ceff, cfrac: ceff*cfrac + leff*(1-cfrac)
distance = lambda x, y: (x.n - y.n)/quad(x.s)
#line = lambda name, val, est: '%20s%20s%20s%20s' % (name, val, est, distance(val, est))
line = lambda name, val, est, est2, est3: ('%20s'*2+'%20.2f'*3) % (name, val, est.n, est2.n, est3.n) #distance(val, est), distance(val, est2))
header = ('%20s'*5) % ('region', 'mc evts', 'from closure', 'lead/sublead', 'avg eff') #, 'delta 2 (sigma)')
header_size = 6*20

for wp, info in wpoints.iteritems():
   table = tables.Table(
      'region:%20s', 
      'mc evts:%20s', 
      ## 'from closure:%20.2f', 
      'lead/sublead:%20.2f', 
      ## 'avg eff:%20.2f',
      'pair stat:%20.2f',
      'pair stat avg:%20.2f',
      title=wp
      )
   
   ceff = 0.5 ## json_ufloat(info, 'charmEff')  #ufloat(info['charmEff'], info['charmEff_err'])
   leff = 0.5 ## json_ufloat(info, 'lightEff')  #ufloat(info['lightEff'], info['lightEff_err'])
   F1   = 0.5 ## factor(leff, ceff, lcfrac)
   F2   = 0.5 ## factor(leff, ceff, scfrac)
   
   F1_mc = factor(mc_effs['light'][wp]['leading'], mc_effs['charm'][wp]['leading'], lcfrac)
   F2_mc = factor(mc_effs['light'][wp]['subleading'], mc_effs['charm'][wp]['subleading'], scfrac)

   F1_avg = factor(mc_effs['light'][wp]['average'], mc_effs['charm'][wp]['average'], lcfrac)
   F2_avg = factor(mc_effs['light'][wp]['average'], mc_effs['charm'][wp]['average'], scfrac)
   ## print 'closure', leff, ceff
   ## print 'lead', mc_effs['light'][wp]['leading'], mc_effs['charm'][wp]['leading']
   ## print 'sublead', mc_effs['light'][wp]['subleading'], mc_effs['charm'][wp]['subleading']
   ## print 'avg', mc_effs['light'][wp]['average'], mc_effs['charm'][wp]['average']

   A = json_ufloat(info, 'both_untagged')
   A_est = nevts*(1 -F1 -F2 +F1*F2)
   A_mc  = nevts*(1 -F1_mc -F2_mc +F1_mc*F2_mc)
   A_avg = nevts*(1 -F1_avg -F2_avg +F1_avg*F2_avg)

   ceff_l = mc_effs['charm'][wp]['leading']
   ceff_s = mc_effs['charm'][wp]['subleading']
   leff_l = mc_effs['light'][wp]['leading'] 
   leff_s = mc_effs['light'][wp]['subleading']
   ceff_a = mc_effs['charm'][wp]['average']
   leff_a = mc_effs['light'][wp]['average'] 
   
   A_pstat = nevts*(
      (1-ceff_l)*(1-leff_s)*lcfrac +
      (1-leff_l)*(1-ceff_s)*scfrac +
      (1-leff_l)*(1-leff_s)*llfrac
      )
   A_pstat_avg = nevts*(
      (1-ceff_a)*(1-leff_a)*lcfrac +
      (1-leff_a)*(1-ceff_a)*scfrac +
      (1-leff_a)*(1-leff_a)*llfrac
      )

   table.add_line(
      'both_untagged', 
      A, 
      ## A_est.n, 
      A_mc.n, 
      ## A_avg.n,
      A_pstat.n,
      A_pstat_avg.n
      )
   
   B = json_ufloat(info, 'lead_tagged')
   B_est = nevts*F1*(1-F2)
   B_mc  = nevts*F1_mc*(1-F2_mc)
   B_avg = nevts*F1_avg*(1-F2_avg)

   B_pstat = nevts*(
      ceff_a*(1-leff_s)*lcfrac +
      leff_l*(1-ceff_s)*scfrac +
      leff_l*(1-leff_s)*llfrac
      )
   B_pstat_avg = nevts*(
      ceff_a*(1-leff_a)*lcfrac +
      leff_a*(1-ceff_a)*scfrac +
      leff_a*(1-leff_a)*llfrac
      )
   table.add_line(
      'lead_tagged', 
      B, 
      ## B_est.n, 
      ## B_mc.n, 
      ## B_avg.n,
      B_pstat.n,
      B_pstat_avg.n
      )

   C = json_ufloat(info, 'sublead_tagged')
   C_est = nevts*F2*(1-F1)
   C_mc  = nevts*F2_mc*(1-F1_mc)
   C_avg = nevts*F2_avg*(1-F1_avg)

   C_pstat = nevts*(
      (1-ceff_l)*leff_s*lcfrac +
      (1-leff_l)*ceff_s*scfrac +
      (1-leff_l)*leff_s*llfrac
      )
   C_pstat_avg = nevts*(
      (1-ceff_a)*leff_a*lcfrac +
      (1-leff_a)*ceff_a*scfrac +
      (1-leff_a)*leff_a*llfrac
      )
   table.add_line(
      'sublead_tagged', 
      C, 
      ## C_est.n, 
      C_mc.n, 
      ## C_avg.n,
      C_pstat.n,
      C_pstat_avg.n
      )

   D = json_ufloat(info, 'both_tagged')
   D_est = nevts*F2*F1
   D_mc  = nevts*F2_mc*F1_mc
   D_avg = nevts*F2_avg*F1_avg

   D_pstat = nevts*(
      ceff_l*leff_s*lcfrac +
      leff_l*ceff_s*scfrac +
      leff_l*leff_s*llfrac
      )
   D_pstat_avg = nevts*(
      ceff_a*leff_a*lcfrac +
      leff_a*ceff_a*scfrac +
      leff_a*leff_a*llfrac
      )
   table.add_line(
      'both_tagged', 
      D, 
      ## D_est.n, 
      D_mc.n, 
      ## D_avg.n,
      D_pstat.n,
      D_pstat_avg.n
      )

   table.add_line(
      'lead notag', 
      A+C, 
      ## (A_est+C_est).n, 
      (A_mc+C_mc).n, 
      ## (A_avg+C_avg).n,
      (A_pstat+C_pstat).n,
      (A_pstat_avg+C_pstat_avg).n
      )
   table.add_line(
      'lead anytag', 
      B+D, 
      ## (B_est+D_est).n, 
      (B_mc+D_mc).n, 
      ## (B_avg+D_avg).n,
      (B_pstat+D_pstat).n,
      (B_pstat_avg+D_pstat_avg).n
      )
   table.add_line(
      'sublead notag', 
      A+B, 
      ## (A_est+B_est).n, 
      (A_mc+B_mc).n, 
      ## (A_avg+B_avg).n,
      (A_pstat+B_pstat).n,
      (A_pstat_avg+B_pstat_avg).n
      )
   table.add_line(
      'sublead anytag', 
      C+D, 
      ## (C_est+D_est).n, 
      (C_mc+D_mc).n, 
      ## (C_avg+D_avg).n,
      (C_pstat+D_pstat).n,
      (C_pstat_avg+D_pstat_avg).n
      )
   table.add_line(
      'total', 
      A+B+C+D, 
      ## (A_est+B_est+C_est+D_est).n, 
      (A_mc+B_mc+C_mc+D_mc).n, 
      ## (A_avg+B_avg+C_avg+D_avg).n,
      (A_pstat+B_pstat+C_pstat+D_pstat).n,
      (A_pstat_avg+B_pstat_avg+C_pstat_avg+D_pstat_avg).n
      )
   print table

