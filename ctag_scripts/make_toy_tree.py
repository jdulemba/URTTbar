from rootpy.io import root_open
from rootpy.tree import Tree
import os
from rootpy import asrootpy

with root_open('toy_tree.root', 'w') as out:
   tree = Tree('toys')
   tree.create_branches({
         'csf' : 'F',
         'lsf' : 'F',
         'nll' : 'F',
         'csf_pfit' : 'F',
         'lsf_pfit' : 'F',
         'itoy' : 'I'
         })
   with root_open('/uscms_data/d3/verzetti/CMSSW_7_1_5/src/URTTbar/plots/2016Feb23/ctageff/mass_discriminant/ctagLoose/toys_simple_initfuzzy2d/mlfitSTICAZZI.root') as mlfit:
      toys = [i.name for i in mlfit.keys() if i.name.startswith('toy_')]               
      for i in toys:            
         path = os.path.join(i, 'fit_s')
         tree.itoy = int(i.split('_')[1])
         toy = mlfit.Get(path)
         pars = asrootpy(toy.floatParsFinal())
         prefit = asrootpy(toy.floatParsInit())
         #nll = mlfit.Get(os.path.join(i, 'nll'))
         
         tree.csf = pars['charmSF'].value
         if tree.csf > 0.9: print tree.itoy
         tree.lsf = pars['lightSF'].value
         #tree.nll = nll.value
         tree.csf_pfit = prefit['charmSF'].value
         tree.lsf_pfit = prefit['lightSF'].value
         tree.fill()
         toy.IsA().Destructor(toy)
   tree.write()
