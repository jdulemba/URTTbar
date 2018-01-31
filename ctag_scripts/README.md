# Instructions to run the ctagging efficiency from scratch

## Getting meta information for new NTuples

First of all you need to gather the info about all the files and samples you have. Each Ntuple production is clustered in one jobid. You can get the current jobid by issuing:
```
echo $jobid
```

To get the info:

```
echo export jobid=SOME_JOBID_TAG > jobid.sh
rake getfiles[USER]  ** user argument is only needed if you're getting files from someone else**
#depending on the number of files, might take a couple of hours
rake meta_batch
#depending on the number of files, could take a very long time
rake getlumi  **make sure symlink to json is correct**
```

Plan on this taking a day or more depending on the number of files.

Get lumi computes lumi for MC only, but creates .run.json files for data needed to compute lumi. **NOT FOR DATA**, for data you need to check every time (and works only on lxplus using BRIL).
Follow the commands from this webpage: https://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html
    or from this TWIKI: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchoolPreExerciseThirdSet#Exercise_16_Combining_the_data_a

Alternatively, download https://github.com/urcms/URBril on lxplus and follow those commands.

```
bash data_pileup.sh
# creates data.meta.pu(_up/_down).root files for all data
# replaces compute_lumi_...sh lines in getting lumi
# 1. check pileup_latest.txt file in data_pileup.sh for correctness
```

Another option is to directly copy the inputs/$jobid directory from someone who already did it.

## Running the analyzers

The analysis proceeds in three steps:
   * Compute shapes for the likelihood discriminant
   * Compute the flavour probabilities
   * Run the proper analyzer

Every new iteration check:
   * Update the lepton SF for Trigger, ID, Isolation
	 * Update Lepton ID/Isolation working points
	 * Update CSV file with the BTagging SF
   * Check that the CSV/cMVA/CTag/DeepCSV/Whatever btagging working points did not change defeinition

These values are found in ctag_eff.cfg [general] section.

### Compute shapes for the likelihood discriminant

This makes the distribution for TTBarSolver

```
rake 'analyze_batch[bin/permProbComputer.cc, ttJets$, ctag_scripts/ctag_eff.cfg]'
python make_permutation_distros.py  'what you want ttsolver filename in cfg file to be'
# update ttsolver file in ctag_eff.cfg
```

### Compute the flavour probabilities

```
rake 'analyze_batch[btag_topology_effs.cc, ttJets$, ctag_scripts/ctag_eff.cfg]'
python ctag_scripts/make_btag_efficiencies.py
```

### Run the proper analyzer

Takes a **long** time
```
rake 'analyze_batch[bin/ctag_eff.cc, *, ctag_scripts/ctag_eff.cfg]'
```

## Making plots

General control plots

```
python ctag_scripts/CTagEffPlotter.py --plots  --shapes --wps="notag" --noPOIpropagation
```

Makes plots and input stuff for the fit for every WP
```
rake ctag_shapes
```

## FIT

Works only in CMSSW 7.1.5 and with the master branch.

You can copy the plots directory from any CMSSW into the one you use for fitting.

Fitting one WP
```
rake ctag_fit[WPName]
```

Fitting all WP
```
rake ctag_fitall
```

Getting the results
```
rake breakdown_all #breaks down the systematics effects
python ctag_scripts/make_ctag_tables.py #this makes a latex file with SF values etc...
python ctag_scripts/make_ctag_postfit.py #postfit plots
python ctag_scripts/write_csv.py ctag|CSV|cMVA #this makes the CSV file to be given to POG
```
