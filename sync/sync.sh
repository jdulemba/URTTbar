pushd ..
scram b -j 8
# htt_simple inputs/SYNC/ttJets.txt sync/ttJets.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync
# htt_simple inputs/SYNC/data_SingleMu.txt sync/data_SingleMu.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync
# htt_simple inputs/SYNC/data_SingleElectron.txt sync/data_SingleElectron.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync
htt_simple inputs/SYNC/ttJets_TRG.txt sync/ttJets.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync
htt_simple inputs/SYNC/data_SingleMu_TRG.txt sync/data_SingleMu.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync
htt_simple inputs/SYNC/data_SingleElectron_TRG.txt sync/data_SingleElectron.sync.root -c htt_scripts/htt_simple.cfg --threads 1 --sync
popd

python <<EOF
from rootpy.io import root_open
with root_open('ttJets.sync.root') as mc:
	with open('mc_mu_sync.list', 'w') as out:
		for i in mc.sync:
			if i.hasMuon == 1:
				out.write('%d:%d:%d\n' % (int(i.run), int(i.lumi), i.evt))
	with open('mc_el_sync.list', 'w') as out:
		for i in mc.sync:
			if i.hasMuon == 0:
				out.write('%d:%d:%d\n' % (int(i.run), int(i.lumi), i.evt))
with root_open('data_SingleMu.sync.root') as mc:
	with open('data_mu_sync.list', 'w') as out:
		for i in mc.sync:
			out.write('%d:%d:%d\n' % (int(i.run), int(i.lumi), i.evt))
with root_open('data_SingleElectron.sync.root') as mc:
	with open('data_el_sync.list', 'w') as out:
		for i in mc.sync:
			out.write('%d:%d:%d\n' % (int(i.run), int(i.lumi), i.evt))
EOF

echo ' #evt file' > summary.raw_txt
wc -l *.list | head -n 4 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Muon MC' >> summary.raw_txt
compareEventLists.py lyon_2016v?_sim_mu_EventIDStep3.txt mc_mu_sync.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Electron MC' >> summary.raw_txt
compareEventLists.py lyon_2016v?_sim_e_EventIDStep3.txt mc_el_sync.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Muon Data' >> summary.raw_txt
compareEventLists.py lyon_2016v?_data_mu_EventIDStep3.txt data_mu_sync.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

echo >> summary.raw_txt
echo >> summary.raw_txt
echo 'Electron Data' >> summary.raw_txt
compareEventLists.py lyon_2016v?_data_e_EventIDStep3.txt data_el_sync.list  --labels='lyon,rochester' | head -n 6 >> summary.raw_txt

cat summary.raw_txt
cp *.list ~/public_html/sync/.
cp *.raw_txt ~/public_html/sync/.
/uscms/home/verzetti/.bin/web.py ~/public_html/sync
