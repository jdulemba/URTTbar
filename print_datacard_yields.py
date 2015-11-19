from URAnalysis.Utilities.datacard import load
from URAnalysis.Utilities.tables import Table
from argparse import ArgumentParser
import os

parser = ArgumentParser()
parser.add_argument('varname')

args = parser.parse_args()
input_file = 'plots/{0}/ttxsec/{1}/{1}.txt'.format(os.environ['jobid'], args.varname)
card = load(input_file)

cols = ['category']+card.processes
cols = [i+':%'+str(len(i)+3)+'s' for i in cols]
yields = Table(*cols, title=args.varname)

categories = card.bins
categories.sort()

base_cat = ''
for cat in categories:
   base = cat.split('_')[0]
   if base != base_cat:
      yields.add_separator()
      base_cat = base
   line = [cat]+['%.1f' % card.rate(cat, proc) for proc in card.processes]
   yields.add_line(*line)
yields.add_separator()

print yields
