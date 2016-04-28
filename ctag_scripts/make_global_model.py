import os
from rootpy.io import root_open
from pdb import set_trace
import URAnalysis.Utilities.prettyjson as prettyjson
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('wp')
args = parser.parse_args()

jobid = os.environ['jobid']
wpdir = 'plots/%s/ctageff/mass_discriminant/%s' % (jobid, args.wp)
datacard = [i for i in open('%s/datacard.txt' % wpdir).readlines()]

toys_dir = '%s/toys' % wpdir
if not os.path.isdir(toys_dir):
   os.mkdir(toys_dir)

mceffs = prettyjson.loads(open('plots/%s/ctageff/info.json'% jobid).read())
json = {
   'mceffs' : mceffs[args.wp],   
   }
new_card = []
sig_norm = None
#fix datacard text
for line in datacard:
   if line.startswith('rate '):
      tokens = line.split()[1:]
      for tok in tokens:
         if tok != '0.00000' and tok != '1.00000':
            line = line.replace(tok, '1.00000')      
   elif ' param ' in line:
      if line.startswith('signal_norm'):
         sig_norm = float(line.split()[2])
      ## tokens = line.split()
      ## if tokens[0].endswith('_eff'):
      ##    name = tokens[0][3:-4]
      ##    value = float(tokens[2])
      ##    json['mceffs'][name] = value
      continue   
   elif line.startswith('CTAG'):
      continue #skip ctag systematics
   new_card.append(line)

#clone and normalize shapes, store yields in json
print "writing to %s" % toys_dir
shapes = root_open('%s/datacard.root' % wpdir)
new_shapes = root_open('%s/datacard.root' % toys_dir, 'w')
categories = [i.name for i in shapes.keys()]
for cname in categories:
   category = shapes.Get(cname)
   new_dir = new_shapes.mkdir(cname)
   samples = []
   shifted = []
   for key in category.keys():
      kname = key.name
      if 'Up' in kname or 'Down' in kname:
         shifted.append(kname)
      else:
         samples.append(kname)

   sample_yields = {}
   for sample in samples:
      shape = category.Get(sample)
      if sample not in json:
         json[sample] = {}
         json[sample]['normalization'] = 0         
      integral = shape.Integral()
      sample_yields[sample] = integral
      if sample != 'right_whad':
         json[sample]['normalization'] += integral
      else:
         json[sample]['normalization'] = sig_norm
      if sample != 'data_obs' and integral != 0: shape.Scale(1./integral)
      new_dir.WriteTObject(shape, sample)
   
   for shift in shifted:
      for i in sample_yields:
         if shift.startswith(i):
            sample = i
            break
      shape = category.Get(shift)
      shape.Scale(1./sample_yields[sample])
      new_dir.WriteTObject(shape, shift)
new_shapes.Close()

flavs = [(1, 'light'), (2, 'charm'), (3, 'beauty')]
fmaps = root_open('plots/%s/ctageff/preselection/falvour_maps.root' % jobid)
for k in fmaps.keys():
   fmap = fmaps.Get(k.name) 
   fmap.Scale(1./fmap.Integral())
   flist = []
   for f1id, f1 in flavs:      
      for f2id, f2 in flavs:
         flist.append([f1, f2, fmap[f1id, f2id].value])
   json[k.name]['flavours'] = flist

with open('%s/datacard.json' % toys_dir, 'w') as j:
   j.write(prettyjson.dumps(json))

with open('%s/datacard.txt' % toys_dir, 'w') as j:
   j.write(''.join(new_card))
   
#text2workspace.py datacard.txt --X-allow-no-background -P URAnalysis.AnalysisTools.statistics.CTagEfficiencies:effcomplete --PO constants=datacard.json -o full_model.root
#combine full_model.root -M MaxLikelihoodFit --skipBOnlyFit --saveToys --expectSignal 1 -t 20 --setPhysicsModelParameters charmSF=1,lightSF=1,beautySF=1
#combine fitModel.root -M MaxLikelihoodFit --saveNormalizations --saveWithUncertainties --saveNLL --skipBOnlyFit --minos=all --toysFile toys_src.root --expectSignal 1 -t 20


## for i in range(toy1.numEntries()):
##    aset = rootpy.asrootpy(toy1.get(i))
##    print aset['CMS_th1x'].value, aset['CMS_channel'].index, toy1.weight()
