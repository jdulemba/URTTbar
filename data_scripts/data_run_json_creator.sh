dir=../inputs/$jobid

golden_json=GoldenJSON_Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.json
echo -e "JSON file used: $golden_json\n"

echo -e "---- Creating final run.json files from  ----"

for i in $(ls ../inputs/$jobid/data_*.run.json); do
        filename=$(basename "$i")
        name="${filename%.*}"
        name="${name%.*}"
        outname="$dir/$name.final_run.json"
        echo "   $i"
        #compareJSON.py --and $golden_json $i $outname
        echo -e "   $outname created\n"

done
