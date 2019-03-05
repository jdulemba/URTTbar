import URAnalysis.Utilities.prettyjson as prettyjson

import os, glob, sys
from pdb import set_trace
import logging as log
log.basicConfig(format='%(message)s', level=log.INFO)

jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']

tt_file = 'ttJetsSL'

input_file = '%s/inputs/%s/%s.meta.json' % (project, jobid, tt_file) 
lumi_file = '%s/inputs/%s/%s.lumi' % (project, jobid, tt_file)

## read contents of input file + create dictionary
tt_dict = prettyjson.loads(open(input_file).read() )
if 'sum_weights' not in tt_dict.keys():
    raise KeyError( "Can't find the correct key to find lhe weights - these are present: " + repr(tt_dict.keys()) )
sum_weights = tt_dict['sum_weights']

## create dictionary to be written into the weights.json file
lhe_weights_dict = {}
for idx in range(len(sum_weights)):
    lhe_weights_dict[idx] = sum_weights[idx]/sum_weights[0]

## write dict into file
#set_trace()
out_fname = '%s.weights.json' % os.path.basename(input_file).split('.')[0]
with open('%s/inputs/%s/%s' % (project, jobid, out_fname), 'w') as out:
    #set_trace()
    out.write(prettyjson.dumps(lhe_weights_dict))

log.info("\n----- %s created in inputs/%s/   -----\n" % (out_fname, jobid))
