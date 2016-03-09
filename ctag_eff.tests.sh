#! /bin/env bash

rm *.ctag_eff.test.root
rake bin/ctag_eff.exe
NPROCS=0
for input in $(ls inputs/$jobid/*.txt | grep -v ttJets_); do
    sample=`basename $input`
    nfiles=`wc -l $input | awk '{print $1}'`
    echo examining $sample '--' $nfiles
    ctag_eff.exe $input $sample.ctag_eff.test.root -c ctag_eff.cfg --threads 1 --J $nfiles  --nosys 1 &
    NPROC=$(($NPROC+1))
    while [ "$NPROC" -ge 10 ]; do
        sleep 3
        NPROC=`ps au | grep $USER | grep ctag_eff.exe | grep -v grep | wc -l`
    done
done
wait

/uscms/home/verzetti/.bin/rootfind.py *txt.ctag_eff.test.root --type=TH1 --name=*nosys*notag/both_untagged/njets --exec='Integral()' > ctag_eff.test.log

#now switch a little bit to python, is easier
python <<EOF
#build yields dict
f = open('ctag_eff.test.log')
sample=''
yields={}
for line in f:
  if line.strip('\x1b[?1034h').startswith('In'):
    sample=line.split()[1].split('.')[0]
  else:
    sp = line.strip().split()
    sub = sp[0].split('nosys')[0].strip('/')
    yields['%s %s' % (sample, sub)] = float(sp[1])

#compare to previous
import os
if os.path.isfile('ctag_eff.test.json'):
  old = eval(open('ctag_eff.test.json').read())
  for key, val in yields.iteritems():
    if key in old:
      oval = old[key]
      diff = '%.0f%%' % (100*(val-oval)/oval) if oval else '%f --> %f' % (oval, val)
      print '%s %s' % (key, diff)
    else:
      print '%s --> NEW' % key

#replace by current
with open('ctag_eff.test.json', 'w') as out:
  out.write(yields.__repr__())

EOF
rm *.ctag_eff.test.root
rm ctag_eff.test.log
