#!/bin/bash
cd ../

echo -e "\nAll ttbar_reco_3J plots will be made from ttJetsM0, ttJetsM700, and ttJetsM1000 samples.\n"

for fname in ../results/$jobid/jet_perm_disc/ttJetsM*;
do
    filename=$(basename $fname .root);
    python jet_perm_disc_plotter.py Full $filename Everything;
done

echo -e "\n     FINISHED\n"

exec bash
