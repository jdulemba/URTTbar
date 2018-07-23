import URAnalysis.Utilities.prettyjson as prettyjson

import os, glob, sys
from pdb import set_trace
import logging as log
log.basicConfig(format='%(message)s', level=log.INFO)

jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']

tt_files = ['ttJetsDiLep', 'ttJetsHad', 'ttJetsSL']

input_files = ['%s/inputs/%s/%s.meta.json' % (project, jobid, x) for x in tt_files ]
lumi_files = ['%s/inputs/%s/%s.lumi' % (project, jobid, x) for x in tt_files ]

for input_file in input_files:
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
    out_fname = '%s.test_weights.json' % os.path.basename(input_file).split('.')[0]
    with open('%s/inputs/%s/%s' % (project, jobid, out_fname), 'w') as out:
        #set_trace()
        out.write(prettyjson.dumps(lhe_weights_dict))
    
    log.info("\n----- %s created in inputs/%s/   -----\n" % (out_fname, jobid))



### create weighted json file from 3 ttJets files

    ## semilep file
SL_file = [x for x in input_files if 'SL' in x][0]
SL_dict = prettyjson.loads(open(SL_file).read() )
SL_lumi_file = [x for x in lumi_files if 'SL' in x][0]
SL_lumi = float(open(SL_lumi_file, 'read').readline())

    ## dilep file
DiLep_file = [x for x in input_files if 'DiLep' in x][0]
DiLep_dict = prettyjson.loads(open(DiLep_file).read() )
DiLep_lumi_file = [x for x in lumi_files if 'DiLep' in x][0]
DiLep_lumi = float(open(DiLep_lumi_file, 'read').readline())

    ## had file
Had_file = [x for x in input_files if 'Had' in x][0]
Had_dict = prettyjson.loads(open(Had_file).read() )
Had_lumi_file = [x for x in lumi_files if 'Had' in x][0]
Had_lumi = float(open(Had_lumi_file, 'read').readline())

    ## create "weights" based on lumi
SL_weight = SL_lumi/(SL_lumi+DiLep_lumi+Had_lumi)
DiLep_weight = DiLep_lumi/(SL_lumi+DiLep_lumi+Had_lumi)
Had_weight = Had_lumi/(SL_lumi+DiLep_lumi+Had_lumi)

SL_sum_weights = SL_dict['sum_weights']
DiLep_sum_weights = DiLep_dict['sum_weights']
Had_sum_weights = Had_dict['sum_weights']

sum_weights = [SL_sum_weights[x]*SL_weight+DiLep_sum_weights[x]*DiLep_weight+Had_sum_weights[x]*Had_weight for x in range(len(SL_sum_weights)) ]

## create dictionary to be written into the weights.json file
lhe_weights_dict = {}
for idx in range(len(sum_weights)):
    lhe_weights_dict[idx] = sum_weights[idx]/sum_weights[0]

## write dict into file
#set_trace()
out_fname = 'ttJets.test_weights.json'
with open('%s/inputs/%s/%s' % (project, jobid, out_fname), 'w') as out:
    #set_trace()
    out.write(prettyjson.dumps(lhe_weights_dict))

log.info("\n----- %s created in inputs/%s/   -----\n" % (out_fname, jobid))


#set_trace()

