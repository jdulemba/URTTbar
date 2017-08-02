for fname in ../results/2017Ap19/ttbar_reco_3J/AtoTT_*;
do
    filename=$(basename $fname .root);
#    python ttbar_reco_3J_plotter.py Full $filename Gen_Plots;
    python ttbar_reco_3J_plotter.py Full $filename Reco_Plots;
#    python ttbar_reco_3J_plotter.py Full $filename Resolution_Plots;
done
