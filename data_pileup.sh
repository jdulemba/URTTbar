dir=inputs/$jobid/INPUT

echo " ---- Creating PILEUP files combining ---- "

for i in $(ls inputs/$jobid/data*.run.json); do
#for i in $(ls inputs/$jobid/data_*to*.run.json); do
        filename=$(basename "$i")
        name="${filename%.*}"
        name="${name%.*}"
        echo "   $name"
done
echo -e "\n  from inputs/$jobid\n"

mergeJSON.py inputs/$jobid/data_*.run.json  --output=data_all.json
pileupCalc.py -i data_all.json --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69000 --maxPileupBin 100 --numPileupBins 100 $dir/data.meta.pu.root

pileupCalc.py -i data_all.json --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 72450 --maxPileupBin 100 --numPileupBins 100 $dir/data.meta.pu_up.root

pileupCalc.py -i data_all.json --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 65550 --maxPileupBin 100 --numPileupBins 100 $dir/data.meta.pu_down.root

echo -e "\n    data.meta.pu.root\n    data.meta.pu_up.root\n    data.meta.pu_down.root\n have been created and placed in $dir"

rm data_all.json
