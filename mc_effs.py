import URAnalysis.Utilities.prettyjson as prettyjson
from glob import glob
import os
from pdb import set_trace
from uncertainties import ufloat
import URAnalysis.Utilities.tables as tables

def compute_eff(jmap, jrank, qtype):
   total = ['both_tagged', 'both_untagged', 'lead_tagged', 'sublead_tagged']
   passing = ['both_tagged']
   if jrank == 'leading': passing.append('lead_tagged')
   else: passing.append('sublead_tagged')
   
   n_total = sum(jmap[i][jrank][qtype] for i in total)
   n_passing = sum(jmap[i][jrank][qtype] for i in passing)
   return n_passing / n_total

jobid = os.environ['jobid']
jsons = glob('plots/%s/btageff/mass_discriminant/*/*/*/flavors_rightw.json' % jobid)

json_ufloat = lambda d, n : ufloat(d[n], d['%s_err' % n])

jobid = os.environ['jobid']
summary = prettyjson.loads(
   open('plots/%s/btageff/summary.json' % jobid).read()
)
wpoints = summary["mass_discriminant"]['working_points']

jmaps = {}
for jfile in jsons:
   wpoint, category, jetrank = tuple(jfile.split('/')[-4:-1])
   if not wpoint in jmaps: jmaps[wpoint] = {}
   if not category in jmaps[wpoint]: jmaps[wpoint][category] = {}
   if not jetrank in jmaps[wpoint][category]: jmaps[wpoint][category][jetrank] = {}

   jmaps[wpoint][category][jetrank] = prettyjson.loads(
      open(jfile).read()
      )

#filter out incomplete (service) working points
categories = ['sublead_tagged', 'lead_tagged', 'both_untagged', 'both_tagged']
to_rm = [i for i in jmaps if jmaps[i].keys() != categories]
for rm in to_rm:
   del jmaps[rm]

effs = {}
for qtype in ['light', 'charm', 'strange']:
   table = tables.Table(
      'WP:%10s', 
      'clead:cards lead:%10.1f',
      'csub:cards sub:%10.1f',
      'lead:lead (%):%10.1f', 
      'sub:sub (%):%10.1f', 
      'avg:avg (%):%10.1f', 
      'diff:rel. diff (%):%15.1f',
      title = '%s efficiencies' % qtype
      )
   effs[qtype] = {}
   for wpoint, jmap in jmaps.iteritems():
      line = table.new_line()
      line.WP = wpoint
      isThere = ("lead_%sEff" % qtype) in wpoints[wpoint]
      line.clead = 100*wpoints[wpoint]["lead_%sEff" % qtype] if isThere else 0.
      line.csub  = 100*wpoints[wpoint]["sub_%sEff"  % qtype]   if isThere else 0.

      leff = compute_eff(jmap, 'leading', qtype)
      seff = compute_eff(jmap, 'subleading', qtype)
      teff = (leff+seff)/2.

      line.lead = leff*100
      line.sub = seff*100
      line.avg = teff*100
      line.diff = 100*(leff-seff)/teff
      effs[qtype][wpoint] = {
         'leading' : leff,
         'subleading' : seff,
         'average' : teff
         }
   del line
   print '\n\n'
   print table

with open('plots/%s/btageff/mc_effs.json' % jobid, 'w') as out:
   out.write(prettyjson.dumps(effs))
   
