from rootpy.io import root_open
from rootpy.plotting import Hist, Hist2D
import os, glob, sys
import URAnalysis.Utilities.prettyjson as prettyjson
from pdb import set_trace

jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']

input_dir = '%s/inputs/data/SFMerging' % project
output_dir = '%s/inputs/data' % project

var_names = [
    #('theJSONfile_RunBtoF_Nov17Nov2017', 'trg', 'IsoMu27_PtEtaBins', 'abseta_pt_DATA'),
    ('RunBCDEF_SF_ID', 'id', 'NUM_TightID_DEN_genTracks', 'abseta_pt'),
    ('RunBCDEF_SF_ISO', 'iso', 'NUM_TightRelIso_DEN_TightIDandIPCut', 'abseta_pt'), 
] #hist name, hist title, kind of info

def make_hist(fname, htype, info_type, hname):
    results = prettyjson.loads(open(fname).read())

    xbins = set()
    ybins = set()
    values_list = []
    errors_list = []
    
    for etaKey, values in sorted(results[htype][info_type].iteritems()):
        xbins.add(float(etaKey[etaKey.find(':')+2:len(etaKey)-1].split(',')[0]))
        xbins.add(float(etaKey[etaKey.find(':')+2:len(etaKey)-1].split(',')[1]))
        for ptKey, result in sorted(values.iteritems()) :
            ybins.add(float(ptKey[ptKey.find(':')+2:len(ptKey)-1].split(',')[0]))
            ybins.add(float(ptKey[ptKey.find(':')+2:len(ptKey)-1].split(',')[1]))
            values_list.append(result["value"])
            errors_list.append(result["error"])
   
    #set_trace() 
    xbins = sorted(list(xbins))
    ybins = sorted(list(ybins))
    H2D = Hist2D(ybins, xbins, name=hname, title=htype)#, xtitle='muon |#eta|', ytitle='muon p_{T} (GeV/c)')
    
    for etaKey, values in sorted(results[htype][info_type].iteritems()):
        xlower = float(etaKey[etaKey.find(':')+2:len(etaKey)-1].split(',')[0])
        for ptKey, result in sorted(values.iteritems()) :
            ylower = float(ptKey[ptKey.find(':')+2:len(ptKey)-1].split(',')[0])
            H2D[ybins.index(ylower)+1,xbins.index(xlower)+1].value = result["value"]
            H2D[ybins.index(ylower)+1,xbins.index(xlower)+1].error = result["error"]

    #set_trace() 
    return H2D    


with root_open('ISO_ID_hists.root', 'w') as out:
#with root_open('%s/ISO_ID_hists.root' % input_dir, 'w') as out:
    for file_path, hname, htitle, var_orientation in var_names:
        print file_path, hname, htitle, var_orientation
        hist = make_hist('%s/%s.json' % (input_dir, file_path), htitle, var_orientation, hname)
        hist.Write()

#        set_trace()

