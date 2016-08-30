#! /bin/env bash

set -o nounset
set -o errexit

wp=$1
wpd=$URA_PROJECT/plots/$jobid/ctageff/mass_discriminant/$wp
toyd=$wpd/toys_simple_initfuzzy2d
mkdir -p $toyd

model=$wpd/fitModel.root
rake $model
#cp $model toyd/.

#
# Build batch script
#
bfile=$toyd/batch_job.sh
echo "#! /bin/bash

pushd $URA_PROJECT" > $bfile  

echo 'source environment.sh
popd

seed=$1
csf=$2
lsf=$3
ntoys=200
combine fitModel.root -M MaxLikelihoodFit --skipBOnlyFit --saveToys --expectSignal 1 -t $ntoys --seed $seed --setPhysicsModelParameters charmSF=$csf,lightSF=$lsf -S 0

mv higgsCombineTest.MaxLikelihoodFit.mH120.$seed.root toys_src.$seed.$csf.$lsf.root
combine fitModel.root -M MaxLikelihoodFit --saveToys --skipBOnlyFit -t $ntoys --minos=all --toysFile toys_src.$seed.$csf.$lsf.root --randomInitPOI lightSF=$lsf,0.1:charmSF=$csf,0.1
mv mlfit.root toys_res.$seed.$csf.$lsf.root
rm toys_src.$seed.$csf.$lsf.root
' >> $bfile

#
# Build condor jdl
#
jdl=$toyd/condor.jdl
echo "universe = vanilla
Executable = batch_job.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Transfer_Input_Files = $model
" > $jdl


val=0
for csf in 0.5 0.7 0.9 1.0 1.1 1.3 1.5; do
    for lsf in 0.5 0.7 0.9 1.0 1.1 1.3 1.5; do
        for seed in 2345678 2445678 2555678; do
            echo output = con_$val.stdout >> $jdl
            echo error = con_$val.stderr >> $jdl
            echo log = con_$val.log >> $jdl
            echo arguments = $seed $csf $lsf >> $jdl
            echo queue >> $jdl
            echo  >> $jdl
            val=`echo "$(($val+1))"`
        done
    done
done

echo $toyd
## pushd $toyd
## condor_submit condor.jdl
## popd
## 
## hold.py
