echo -e "This will create plots from each part of the analysis (gen_partons, jet_perm_dsic, ttbar_reco_3J).\n\nAll gen_partons plots will be made from ttJetsM0, ttJetsM700, and ttJetsM1000 samples.\n"

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

echo -e "\nAll jet_perm_disc plots will be made from ttJetsM0, ttJetsM700, and ttJetsM1000 samples.\n"

for fname in ../results/$jobid/jet_perm_disc/ttJetsM*;
do
    filename=$(basename $fname .root);
    python jet_perm_disc_plotter.py Full $filename Everything;
done

echo -e "\nAll ttbar_reco_3J plots will be made from ttJetsM0, ttJetsM700, and ttJetsM1000 samples.\n"

for fname in ../results/$jobid/ttbar_reco_3J/ttJetsM*;
do
    filename=$(basename $fname .root);
    python ttbar_reco_3J_plotter.py Full $filename Everything;
done

echo -e "\nAll ttbar_reco_3J plots will be made from AtoTT samples.\n"

for fname in ../results/$jobid/ttbar_reco_3J/AtoTT*;
do
    filename=$(basename $fname .root);
    python ttbar_reco_3J_plotter.py Full $filename Everything;
done

echo -e "\nAll ttbar_reco_3J plots will be made from HtoTT samples.\n"

for fname in ../results/$jobid/ttbar_reco_3J/HtoTT*;
do
    filename=$(basename $fname .root);
    python ttbar_reco_3J_plotter.py Full $filename Everything;
done

echo -e "\nThis has finished."
