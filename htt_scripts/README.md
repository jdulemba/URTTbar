## Instructions for running ttbar reconstruction from scratch

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
rake getlumi
```

Plan on this taking a day or more depending on the number of files.

Get lumi computes lumi for MC only. **NOT FOR DATA**, for data you need to check every time (and works only on lxplus using BRIL).
Follow the commands from this webpage: https://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html
    or from this TWIKI: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideCMSDataAnalysisSchoolPreExerciseThirdSet#Exercise_16_Combining_the_data_a

Alternatively, download https://github.com/urcms/URBril on lxplus and follow those commands.

Another option is to directly copy the inputs/$jobid directory from someone who already did it.

## Running the analyzers

The analysis has at least four steps:
    * Compute shapes for likelihood discriminant
    * Run linear discriminant analysis
    * Run proper analyzer

### Compute shapes for the likelihood discriminant

This makes the distribution for TTBarSolver

```
rake 'analyze_batch[bin/permProbComputer.cc, ttJets$, cfg_files/ttbar_reco_3J.cfg]'   # ttJets file could be ttJetsM0, ttJetsM700, or ttJetsM1000
python make_3J_permutation_distros.py [output file name]  # makes root file with likelihood shapes
python htt_scripts/permProbComputer_plotter.py  # makes formatted hists for 3 jet lost/merged events
```

### Compute distributions for merged/lost jet events with matched permutations to be used for linear discriminant analysis

```
rake 'test[bin/ttbar_matched_perms.cc, ttJetsM0$, cfg_files/ttbar_reco_3J.cfg]'  -> use    'time ... ' command to run over multiple files to get enough statistics

# copy file locally to run LDA
python LinearDiscr_Split.py --file=[ttJets sample] --disc_type=Total --comp_type=event
# use results from json in both_bp_cats() function in ttbar_bestperm_solutions.cc
rake 'analyze_batch[bin/ttbar_bestperm_solutions.cc, ttJetsM0$, cfg_files/ttbar_reco_3J.cfg]'
python htt_scripts/ttbar_reco_3J_plotter.py --analysis=Full --sample=ttJetsM0 --perm=Best_Perm --plot=All


# find thad mass cut to use in opposite classification events in both_bp_cats function in ttbar_final_reco_3J.cc
rake 'analyze_batch[bin/ttbar_final_reco_3J.cc, ttJetsM0$, cfg_files/ttbar_reco_3J.cfg]'

# use alpha_E and alpha_P linear fit results from alpha_correction_plots
# rescale reco thad in post_alpha_correction() by whichever alpha is consistent across mass ranges
rake 'analyze_batch[bin/ttbar_bestperm_solutions.cc, ttJetsM0$, cfg_files/ttbar_reco_3J.cfg]'



