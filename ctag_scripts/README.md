# Instructions to run the ctagging efficiency from scratch

## Getting meta information for new NTuples

First of all you need to gather the info about all the files and samples you have. Each Ntuple production is clustered in one jobid. You can get the current jobid by issuing:
```
echo $jobid
```

To get the info:

```
echo export jobid=SOME_JOBID_TAG > jobid.sh
rake getfiles[USER, sample]  ** if user==group then lpcbtag will be checked else a personal eos is used **
#depending on the number of files, might take a couple of hours
rake meta_batch   OR  meta_batch[sample]    for individual samples (can only be used for one at a time)
#depending on the number of files, could take a very long time
rake getlumi
```

Plan on this taking a day or more depending on the number of files.

Get lumi computes lumi for MC only. **NOT FOR DATA**
Lumi ** FOR DATA **
Follow the commands from this webpage: https://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html
    or from this TWIKI: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchoolPreExerciseThirdSet#Exercise_16_Combining_the_data_a

Alternatively, download https://github.com/urcms/URBril on lxplus and follow those commands.

Another option is to directly copy the inputs/$jobid directory from someone who already did it.

```
    ** FOR DATA **
To create final_run.json files,
check golden json file for correctness (/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt) and run

bash data_run_json_creator.sh

Pileup ** FOR DATA **
Check input files (final_run.json and pileup_latest.txt) for correctness.

bash data_pileup.sh
# creates data.meta.pu(_up/_down).root files for all data
# replaces compute_lumi_...sh lines in getting lumi
# 1. check pileup_latest.txt file in data_pileup.sh for correctness
#    and make sure it's actually the latest version (from /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PileUp/)
```

** If new ntuples have different branches, must run rake 'proxy[ttJets]' **


## Running the analyzers

The analysis proceeds in three steps:
   * Compute shapes for the likelihood discriminant
   * Compute the flavour probabilities
   * Run the proper analyzer

Every new iteration check:
   * Update the lepton SF for Trigger, ID, Isolation
	 * Update Lepton ID/Isolation working points
	    * Muons: https://twiki.cern.ch/twiki/bin/view/CMS/MuonPOG
	    * Electrons: 
	 * Update CSV file with the BTagging SF: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation
	    * MUST CHANGE 'comb' and 'incl' operating points to 'used' 
   * Check that the CSV/cMVA/CTag/Whatever btagging working points did not change defeinition
   * Check MET Filter recommendations: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
   * Check Triggers: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016 (2017, etc...)


These values are found in ctag_eff.cfg [general] section.

### Compute shapes for the likelihood discriminant

This makes the distribution for TTBarSolver

```
rake 'test[bin/permProbComputer.cc, ttJetsSL$, ctag_scripts/ctag_eff.cfg]'

rake 'analyze_batch[bin/permProbComputer.cc, ttJetsSL$, ctag_scripts/ctag_eff.cfg]'
python make_permutation_distros.py 'outname'
```

### Compute the flavour probabilities

```
rake 'test[bin/btag_topology_effs.cc, ttJetsSL$, ctag_scripts/ctag_eff.cfg]'

rake 'analyze_batch[bin/btag_topology_effs.cc, ttJets[SL, Had, DiLep]$, ctag_scripts/ctag_eff.cfg]'
python ctag_scripts/make_btag_efficiencies.py
```

### Run the proper analyzer

Takes a **long** time
```
rake 'test[bin/ctag_eff.cc, ttJetsSL$, ctag_scripts/ctag_eff.cfg]'

rake 'analyze_batch[bin/ctag_eff.cc, *, ctag_scripts/ctag_eff.cfg]'
```

## Making plots

General control plots

```
python CTagEffPlotter.py --plots  --shapes --wps="notag" --noPOIpropagation
```

Makes plots and input stuff for the fit for every WP
```
rake ctag_plotfit[DeepCSV or DeepJet]
```

## FIT

MAKE SURE jobid and branches ARE CORRECT!!!

Must be run in CMSSW 8.1.0 with 810_miniaod_2018ctag_combine branch (URTTbar) and ura_810_2018 in URAnalysis

You can copy the plots directory from any CMSSW into the one you use for fitting or use a symbolic link.

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
