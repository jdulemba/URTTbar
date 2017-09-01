for fname in ../results/2017Ap19/gen_partons/ttJetsM*;
do
    filename=$(basename $fname .root);
    python gen_partons_plotter.py Full $filename AllJets;
    python gen_partons_plotter.py Full $filename 3J;
    python gen_partons_plotter.py Full $filename 4J;
    python gen_partons_plotter.py Full $filename 5PJ;
done

