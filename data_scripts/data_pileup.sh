dir=../inputs/$jobid/INPUT

echo "---- Creating PILEUP files combining ----"

for i in $(ls ../inputs/$jobid/data_*.final_run.json); do
        filename=$(basename "$i")
        name="${filename%.*}"
        name="${name%.*}"
        echo "   $name"
done
echo -e "\n  from ../inputs/$jobid\n"

pileup_latest=2016legacy_pileup_latest_14Jan2019.txt
echo -e "Pileup file used: $pileup_latest\n"

mergeJSON.py ../inputs/$jobid/data_*.final_run.json  --output=data_all.json

pileupCalc.py -i data_all.json --inputLumiJSON $pileup_latest --calcMode true --minBiasXsec 69000 --maxPileupBin 100 --numPileupBins 100 $dir/data.meta.pu.root
echo -e "\n    data.meta.pu.root\n"

pileupCalc.py -i data_all.json --inputLumiJSON $pileup_latest --calcMode true --minBiasXsec 72450 --maxPileupBin 100 --numPileupBins 100 $dir/data.meta.pu_up.root
echo -e "\n    data.meta.pu_up.root\n"

pileupCalc.py -i data_all.json --inputLumiJSON $pileup_latest --calcMode true --minBiasXsec 65550 --maxPileupBin 100 --numPileupBins 100 $dir/data.meta.pu_down.root
echo -e "\n    data.meta.pu_down.root\n have been created and placed in $dir"

rm data_all.json
