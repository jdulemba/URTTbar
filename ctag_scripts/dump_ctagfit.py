import rootpy.io as io
from argparse import ArgumentParser
import rootpy
from pdb import set_trace

parser = ArgumentParser()
parser.add_argument('files', nargs='+')
parser.add_argument('--pulls', action='store_true')
args = parser.parse_args()

def pvar(pars, vname):
   if vname not in pars: return
   var = pars[vname]
   err = var.error
   if hasattr(err, '__len__'):
      val = '%20s   %+.3f   %+.3f/%.3f' % (vname, var.value, err[0], err[1])
   else:
      val = '%20s   %+.3f   +/-%.3f' % (vname, var.value, err)
   return val

class CorrelationMatrix(object):
   def __init__(self, matrix):
      self.matrix = matrix
      self.xmapping = {matrix.xaxis.GetBinLabel(idx) : idx for idx in range(1, matrix.nbins()+1)}
      self.ymapping = {matrix.yaxis.GetBinLabel(idx) : idx for idx in range(1, matrix.nbins()+1)}
      self.nbins = matrix.nbins()

   def __getitem__(self, names):
      v1, v2 = names
      idx1 = self.xmapping[v1]
      idx2 = self.ymapping[v2]
      return self.matrix[idx1,idx2].value

for name in args.files:
   tf = io.root_open(name)
   pars = rootpy.asrootpy(tf.fit_s.floatParsFinal())
   print name
   print pvar(pars, 'charmSF')
   print pvar(pars, 'lightSF')
   
   if args.pulls:      
      correlation = CorrelationMatrix(rootpy.asrootpy(tf.fit_s.correlationHist()))
      print '\n\n--- PULLS ---'
      for par in pars:
         if par.name.endswith('SF'): continue
         print pvar(pars, par.name), '     %+.3f' % correlation['charmSF', par.name]
