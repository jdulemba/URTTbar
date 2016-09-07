cp summary.raw_txt summary.old.raw_txt

pushd ..
scram b -j 8
# htt_simple inputs/SYNC/ttJets.txt sync/ttJets.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync > sync/tt_step2.out
# htt_simple inputs/SYNC/data_SingleMu.txt sync/data_SingleMu.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync > sync/mu_step2.out
# htt_simple inputs/SYNC/data_SingleElectron.txt sync/data_SingleElectron.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync > sync/el_step2.out
htt_simple inputs/SYNC/ttJets_noskim.txt sync/ttJets.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync > sync/tt_step2.out
htt_simple inputs/SYNC/data_SingleMu_noskim.txt sync/data_SingleMu.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync > sync/mu_step2.out
htt_simple inputs/SYNC/data_SingleElectron_noskim.txt sync/data_SingleElectron.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync > sync/el_step2.out
popd

python <<EOF
from rootpy.io import root_open
with root_open('ttJets.sync.root') as mc:
	with open('mc_mu_sync_s3.list', 'w') as out:
		for i in mc.sync:
			if i.hasMuon == 1:
				out.write('%d:%d:%d\n' % (int(i.run), int(i.lumi), i.evt))
	with open('mc_el_sync_s3.list', 'w') as out:
		for i in mc.sync:
			if i.hasMuon == 0:
				out.write('%d:%d:%d\n' % (int(i.run), int(i.lumi), i.evt))
with root_open('data_SingleMu.sync.root') as mc:
	with open('data_mu_sync_s3.list', 'w') as out:
		for i in mc.sync:
			out.write('%d:%d:%d\n' % (int(i.run), int(i.lumi), i.evt))
with root_open('data_SingleElectron.sync.root') as mc:
	with open('data_el_sync_s3.list', 'w') as out:
		for i in mc.sync:
			out.write('%d:%d:%d\n' % (int(i.run), int(i.lumi), i.evt))
EOF

echo '#####################################################' >  summary.raw_txt
echo '                    STEP 3 (final)                   ' >> summary.raw_txt
echo '#####################################################' >> summary.raw_txt
echo ' #evt file' >> summary.raw_txt
wc -l *_s3.list | head -n 4 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Muon MC' >> summary.raw_txt
compareEventLists.py lyon_2016v?_sim_mu_EventIDStep3.txt mc_mu_sync_s3.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Electron MC' >> summary.raw_txt
compareEventLists.py lyon_2016v?_sim_e_EventIDStep3.txt mc_el_sync_s3.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Muon Data' >> summary.raw_txt
compareEventLists.py lyon_2016v?_data_mu_EventIDStep3.txt data_mu_sync_s3.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Electron Data' >> summary.raw_txt
compareEventLists.py lyon_2016v?_data_e_EventIDStep3.txt data_el_sync_s3.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

echo '#####################################################' >> summary.raw_txt
echo '                    STEP 2 (lep selection)           ' >> summary.raw_txt
echo '#####################################################' >> summary.raw_txt
echo ' #evt file' >> summary.raw_txt

grep '1 ' tt_step2.out | awk '{print $2}' > mc_mu_sync_s2.list
grep '0 ' tt_step2.out | awk '{print $2}' > mc_el_sync_s2.list
grep '1 ' mu_step2.out | awk '{print $2}' > data_mu_sync_s2.list
grep '0 ' el_step2.out | awk '{print $2}' > data_el_sync_s2.list 
wc -l *_s2.list | head -n 4 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Muon MC' >> summary.raw_txt
compareEventLists.py lyon_2016v?_sim_mu_EventIDStep2.txt mc_mu_sync_s2.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Electron MC' >> summary.raw_txt
compareEventLists.py lyon_2016v?_sim_e_EventIDStep2.txt mc_el_sync_s2.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Muon Data' >> summary.raw_txt
compareEventLists.py lyon_2016v?_data_mu_EventIDStep2.txt data_mu_sync_s2.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Electron Data' >> summary.raw_txt
compareEventLists.py lyon_2016v?_data_e_EventIDStep2.txt data_el_sync_s2.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

cat summary.raw_txt
cp *.list ~/public_html/sync/.
cp summary.raw_txt ~/public_html/sync/.
/uscms/home/verzetti/.bin/web.py ~/public_html/sync
