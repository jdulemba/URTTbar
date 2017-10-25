echo -e "This script creates all of the plots from gen_partons, jet_perm_disc, and ttbar_reco_3J for the ttJetsM0, ttJetsM700, and ttJetsM1000 samples.\n"

echo -e "\nAll gen_partons plots will be made from ttJetsM0, ttJetsM700, and ttJetsM1000 samples.\n"

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
