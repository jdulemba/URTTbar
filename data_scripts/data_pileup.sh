outdir=$URA_PROJECT/inputs/$jobid/INPUT

echo "---- Creating PILEUP files combining ----"

for i in $(ls $URA_PROJECT/inputs/$jobid/data_*.run.json); do
        filename=$(basename "$i")
        name="${filename%.*}"
        name="${name%.*}"
        echo "   $name"
done
echo -e "\n  from $URA_PROJECT/inputs/$jobid\n"

pileup_latest=$URA_PROJECT/data_scripts/pileup_latest_2018.txt
echo -e "Pileup file used: $pileup_latest\n"

mergeJSON.py $URA_PROJECT/inputs/$jobid/data_*.run.json  --output=data_all.json

./working_puCalc.py -i data_all.json --inputLumiJSON $pileup_latest --calcMode true --minBiasXsec 69000 --maxPileupBin 100 --numPileupBins 100 $outdir/data.meta.pu.root
echo -e "\n    data.meta.pu.root\n"

./working_puCalc.py -i data_all.json --inputLumiJSON $pileup_latest --calcMode true --minBiasXsec 72450 --maxPileupBin 100 --numPileupBins 100 $outdir/data.meta.pu_up.root
echo -e "\n    data.meta.pu_up.root\n"

./working_puCalc.py -i data_all.json --inputLumiJSON $pileup_latest --calcMode true --minBiasXsec 65550 --maxPileupBin 100 --numPileupBins 100 $outdir/data.meta.pu_down.root
echo -e "\n    data.meta.pu_down.root\n have been created and placed in $outdir"

rm data_all.json
