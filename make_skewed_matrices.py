import rootpy
import rootpy.plotting as plotting
import rootpy.io as io
import ROOT
rootpy.log["/"].setLevel(rootpy.log.INFO)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)
from argparse import ArgumentParser
import os

parser = ArgumentParser()
parser.add_argument('var', type=str, help='varible to use')
parser.add_argument('bias_id', type=str, help='ID string of bias to apply')
opts = parser.parse_args()

jobid = os.environ['jobid']

matrices_file_name = 'plots/' + jobid + '/ttxsec/migration_matrices.root'
skewed_matrices_file_name = 'plots/' + jobid + '/ttxsec_' + opts.bias_id + '/migration_matrices.root'

matrices_file = io.root_open(matrices_file_name, "UPDATE")
skewed_matrices_file = io.root_open(skewed_matrices_file_name, "READ")

matrices_file.cd()

new_dir = opts.var + '_' + opts.bias_id


#if not 
matrices_file.mkdir(new_dir)

#matrices_file.cd(opts.var)

olddir = getattr(matrices_file,opts.var)
newdir = getattr(matrices_file, new_dir)
skeweddir = getattr(skewed_matrices_file, opts.var)
olddir.cd()

for key in olddir.keys():
    obj = key.ReadObj()
    ROOT.SetOwnership(obj, False)
    if not obj.GetName().startswith("true_distribution"):
        newdir.cd()
        #matrices_file.cd("../" + new_dir)
        #newobj = obj.Clone()
        #ROOT.SetOwnership(newobj, False)
        #newdir.WriteTObject(newobj, newobj.GetName())
        #newdir.WriteTObject(obj, obj.GetName())
        obj.Write()
        #newobj.Write()
        olddir.cd()
        #matrices_file.cd("../" + opts.var)

#skewed_matrices_file.cd(opts.var)


obj = skeweddir.Get('true_distribution')
newdir.cd()
#newobj = obj.clone()
#ROOT.SetOwnership(newobj, False)
#newdir.WriteTObject(newobj, newobj.GetName())
#newdir.WriteTObject(obj, obj.GetName())
obj.Write()
#newobj.Write()
skeweddir.cd()
#skewed_matrices_file.cd("../" + opts.var)

skewed_matrices_file.Close()
matrices_file.Close()
