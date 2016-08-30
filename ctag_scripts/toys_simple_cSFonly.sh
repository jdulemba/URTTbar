#! /bin/env bash

set -o nounset
set -o errexit

wp=$1
wpd=$URA_PROJECT/plots/$jobid/ctageff/mass_discriminant/$wp
toyd=$wpd/toys_wsys
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
ntoys=200
combine fitModel.root -M MaxLikelihoodFit --skipBOnlyFit --saveToys --expectSignal 1 -t $ntoys --seed -1 --setPhysicsModelParameters charmSF=$csf

mv higgsCombineTest.MaxLikelihoodFit.mH120.*.root toys_src.$seed.$csf.root
combine fitModel.root -M MaxLikelihoodFit --saveToys --skipBOnlyFit -t $ntoys --minos=all --toysFile toys_src.$seed.$csf.root
 #--randomInitPOI charmSF=$csf,0.1
mv mlfit.root toys_res.$seed.$csf.root
rm toys_src.$seed.$csf.root
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
for csf in 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5; do
    for seed in 1 2 3 4 5; do
        echo output = con_$val.stdout >> $jdl
        echo error = con_$val.stderr >> $jdl
        echo log = con_$val.log >> $jdl
        echo arguments = $seed $csf >> $jdl
        echo queue >> $jdl
        echo  >> $jdl
        val=`echo "$(($val+1))"`
    done
done

echo $toyd
## pushd $toyd
## condor_submit condor.jdl
## popd
## 
## hold.py
