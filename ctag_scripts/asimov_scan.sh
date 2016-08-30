#! /bin/env bash  

set -o nounset

wp=$1
wpd=$URA_PROJECT/plots/$jobid/ctageff/mass_discriminant/$wp
toyd=$wpd/asimovs
mkdir -p $toyd

model=$wpd/fitModel.root
rake $model
cp $model $toyd/.

pushd $toyd

for csf in 0.5 0.7 0.9 1.0 1.1 1.3 1.5; do
    for lsf in 0.5 0.7 0.9 1.0 1.1 1.3 1.5; do        
        echo scanning $csf $lsf
        combine ../fitModel.root -M MaxLikelihoodFit --skipBOnlyFit --saveToys --expectSignal 1 -t -1 --setPhysicsModelParameters charmSF=$csf,lightSF=$lsf &> toy_$csf.$lsf.log
        ecode=$?
        if [ $ecode -ne 0 ]; then
            echo '  --> failed'
        else
            mv higgsCombineTest.MaxLikelihoodFit.mH120.123456.root toys_src.$csf.$lsf.root 
            combine ../fitModel.root -M MaxLikelihoodFit --saveToys --skipBOnlyFit -t -1 --minos=all --toysFile toys_src.$csf.$lsf.root &> fit_$csf.$lsf.log
            mv mlfit.root toys_res.123456.$csf.$lsf.root
        fi
    done
done

popd