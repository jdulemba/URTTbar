#!/bin/bash
cd ../

for fname in ../results/$jobid/ttbar_reco_3J/AtoTT*;
do
    filename=$(basename $fname .root);
    python ttbar_reco_3J_plotter.py Full $filename Everything;
done

exec bash
