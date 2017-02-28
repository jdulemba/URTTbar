#! /bin/env bash
set -o nounset
set -o errexit
cp summary.raw_txt summary.old.raw_txt
cp /afs/cern.ch/work/m/mgul/public/Hto_ttbar/deepCSV8024/CMSSW_8_0_24/src/Tupel/Tupel/synchro_dir/*_event_synchro_step[23].txt ughent/.
cp /afs/cern.ch/work/a/aapopov/public/tmp/HToTTSync/events_v2/*.* lyon/.

pushd ..
scram b -j 8
htt_simple inputs/SYNC/ttJets.txt sync/ttJets.sync.root -c htt_scripts/htt_baseline_j20.cfg --threads 1 --sync > sync/tt_step2.out
htt_simple inputs/SYNC/data_SingleMu.txt sync/data_SingleMu.sync.root -c htt_scripts/htt_baseline_j20.cfg --threads 1 --sync > sync/mu_step2.out
htt_simple inputs/SYNC/data_SingleMuRECO.txt sync/data_SingleMuRECO.sync.root -c htt_scripts/htt_baseline_j20.cfg --threads 1 --sync > sync/muRECO_step2.out
# htt_simple inputs/SYNC/data_SingleElectron.txt sync/data_SingleElectron.sync.root -c htt_scripts/htt_baseline_j20.cfg --threads 1 --sync > sync/el_step2.out
# htt_simple inputs/SYNC/ttJets_noskim.txt sync/ttJets.sync.root -c htt_scripts/htt_baseline_j20.cfg --threads 1 --sync > sync/tt_step2.out
# htt_simple inputs/SYNC/data_SingleMu_noskim.txt sync/data_SingleMu.sync.root -c htt_scripts/htt_baseline_j20.cfg --threads 1 --sync > sync/mu_step2.out
# htt_simple inputs/SYNC/data_SingleElectron_noskim.txt sync/data_SingleElectron.sync.root -c htt_scripts/htt_baseline_j20.cfg --threads 1 --sync > sync/el_step2.out
popd

python <<EOF
from rootpy.io import root_open
with root_open('ttJets.sync.root') as mc:
	with open('mc_mu_sync_s3.list', 'w') as out:
		for i in mc.sync:
			if i.hasMuon == 1:
				out.write('%d:%d:%d\n' % (int(i.Run), int(i.LumiSection), i.Event))
	with open('mc_el_sync_s3.list', 'w') as out:
		for i in mc.sync:
			if i.hasMuon == 0:
				out.write('%d:%d:%d\n' % (int(i.Run), int(i.LumiSection), i.Event))
with root_open('data_SingleMu.sync.root') as mc:
	with open('data_mu_sync_s3.list', 'w') as out:
		for i in mc.sync:
			out.write('%d:%d:%d\n' % (int(i.Run), int(i.LumiSection), i.Event))
with root_open('data_SingleMuRECO.sync.root') as mc:
	Pass = open('data_mu_RECO_pass.list', 'w')
	fail = open('data_mu_RECO_fail.list', 'w')
	with open('data_mu_RECO_s3.list', 'w') as out:
		for i in mc.sync:
			idx = '%d:%d:%d\n' % (int(i.Run), int(i.LumiSection), i.Event)
			out.write(idx)
			if i.RecoSuccess:
				Pass.write(idx)
			else:
				fail.write(idx)
## with root_open('data_SingleElectron.sync.root') as mc:
## 	with open('data_el_sync_s3.list', 'w') as out:
## 		for i in mc.sync:
## 			out.write('%d:%d:%d\n' % (int(i.Run), int(i.LumiSection), i.Event))
EOF

echo '#####################################################' >  summary.raw_txt
echo '                      RECO STEP                      ' >> summary.raw_txt
echo '#####################################################' >> summary.raw_txt

./compare_sync_tree.py data_SingleMuRECO.sync.root lyon/Andrey_v2.root --labels='Rochester,Lyon' --nolist >> summary.raw_txt
echo >> summary.raw_txt

./compare_sync_tree.py data_SingleMuRECO.sync.root ughent/ntuplesynchro_v2.root --labels='Rochester,Ghent' --nolist >> summary.raw_txt
echo >> summary.raw_txt
echo >> summary.raw_txt

echo '#####################################################' >> summary.raw_txt
echo '                    STEP 3 (final)                   ' >> summary.raw_txt
echo '#####################################################' >> summary.raw_txt
echo ' #evt file' >> summary.raw_txt
wc -l *_s3.list | head -n 4 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Muon MC' >> summary.raw_txt
compareEventLists.py ughent/mc_mu_event_synchro_step3.txt mc_mu_sync_s3.list  --labels='UGhent,rochester'  --nolist  >> summary.raw_txt
compareEventLists.py lyon/lyon_sim_mu_EventIDStep3.txt mc_mu_sync_s3.list  --labels='Lyon,rochester'  --nolist  >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Electron MC' >> summary.raw_txt
compareEventLists.py ughent/mc_e_event_synchro_step3.txt mc_el_sync_s3.list  --labels='UGhent,rochester' --nolist >> summary.raw_txt
compareEventLists.py lyon/lyon_sim_e_EventIDStep3.txt mc_el_sync_s3.list  --labels='Lyon,rochester' --nolist >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Muon Data' >> summary.raw_txt
compareEventLists.py ughent/data_mu_event_synchro_step3.txt data_mu_sync_s3.list  --labels='UGhent,rochester' --nolist >> summary.raw_txt
compareEventLists.py lyon/lyon_data_mu_EventIDStep3.txt data_mu_sync_s3.list  --labels='Lyon,rochester' --nolist >> summary.raw_txt

## echo >> summary.raw_txt
## echo >> summary.raw_txt
## echo 'Electron Data' >> summary.raw_txt
## compareEventLists.py lyon_2016v?_data_e_EventIDStep3.txt data_el_sync_s3.list  --labels='lyon,rochester' --nolist >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo '#####################################################' >> summary.raw_txt
echo '                    STEP 2 (lep selection)           ' >> summary.raw_txt
echo '#####################################################' >> summary.raw_txt
echo ' #evt file' >> summary.raw_txt

grep '1 ' tt_step2.out | awk '{print $3}' > mc_mu_sync_s2.list
grep '0 ' tt_step2.out | awk '{print $3}' > mc_el_sync_s2.list
grep '1 ' mu_step2.out | awk '{print $3}' > data_mu_sync_s2.list
grep '1 ' muRECO_step2.out | awk '{print $3}' > data_mu_RECO_s2.list
## grep '0 ' el_step2.out | awk '{print $3}' > data_el_sync_s2.list 
wc -l *_s2.list | head -n 4 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Muon MC' >> summary.raw_txt
compareEventLists.py ughent/mc_mu_event_synchro_step2.txt mc_mu_sync_s2.list  --labels='UGhent,rochester' --nolist >> summary.raw_txt
compareEventLists.py lyon/lyon_sim_mu_EventIDStep2.txt mc_mu_sync_s2.list  --labels='Lyon,rochester' --nolist >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Electron MC' >> summary.raw_txt
compareEventLists.py ughent/mc_e_event_synchro_step2.txt mc_el_sync_s2.list  --labels='UGhent,rochester' --nolist >> summary.raw_txt
compareEventLists.py lyon/lyon_sim_e_EventIDStep2.txt mc_el_sync_s2.list  --labels='Lyon,rochester' --nolist >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Muon Data' >> summary.raw_txt
compareEventLists.py ughent/data_mu_event_synchro_step2.txt data_mu_sync_s2.list  --labels='UGhent,rochester' --nolist >> summary.raw_txt
compareEventLists.py lyon/lyon_data_mu_EventIDStep2.txt data_mu_sync_s2.list  --labels='Lyon,rochester' --nolist >> summary.raw_txt
## 
## echo >> summary.raw_txt
## echo >> summary.raw_txt
## echo 'Electron Data' >> summary.raw_txt
## compareEventLists.py lyon_2016v?_data_e_EventIDStep2.txt data_el_sync_s2.list  --labels='lyon,rochester' --nolist >> summary.raw_txt

cat summary.raw_txt
cp *.list ~/public_html/sync/.
cp data_SingleMuRECO.sync.root ~/public_html/sync/.
cp summary.raw_txt ~/public_html/sync/.
/uscms/home/verzetti/.bin/web.py ~/public_html/sync
