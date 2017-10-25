#!/bin/bash
cd ../ 

echo -e "All gen_partons plots will be made from ttJetsM0, ttJetsM700, and ttJetsM1000 samples.\n"

for fname in ../results/$jobid/gen_partons/ttJetsM*;
do
    filename=$(basename $fname .root);
    python gen_partons_plotter.py Full $filename Everything;
done

echo -e "\n3J, 4J, 5PJ, and AllJets plots for the combined ttJetsM samples.\n"

python gen_partons_plotter.py Full Combined 3J;
python gen_partons_plotter.py Full Combined 4J;
python gen_partons_plotter.py Full Combined 5PJ;
python gen_partons_plotter.py Full Combined AllJets;


echo -e "\n     FINISHED\n"

exec bash
