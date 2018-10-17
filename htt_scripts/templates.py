import os, glob
import re
import rootpy.io as io
import rootpy
from pdb import set_trace
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--njets', help='Select between 3 or 4+ jet categories')

args = parser.parse_args()


jobid = os.environ['jobid']
project = os.environ['URA_PROJECT']
outdir = '%s/plots/%s/htt/templates' % (project, jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

njets = '3Jets' if args.njets == '3' else '4PJets'

mu_file = '%s/plots/%s/htt/muons/%s/muons.root' % (project, jobid, njets)
if not os.path.isfile(mu_file):
    raise IOError('File: %s does not exist' % mu_file)

el_file = '%s/plots/%s/htt/electrons/%s/electrons.root' % (project, jobid, njets)
if not os.path.isfile(el_file):
    raise IOError('File: %s does not exist' % el_file)

mu_tfile = io.root_open(mu_file)

mu_sig = [ sys.name for sys in mu_tfile.mujets.keys() if re.match('gg[AH]', sys.name) ]
mu_bkg = [ sys.name for sys in mu_tfile.mujets.keys() if not re.match('gg[AH]', sys.name) ]

set_trace()
el_tfile = io.root_open(el_file)

el_sig = [ sys.name for sys in el_tfile.ejets.keys() if re.match('gg[AH]', sys.name) ]
el_bkg = [ sys.name for sys in el_tfile.ejets.keys() if not re.match('gg[AH]', sys.name) ]

with io.root_open('%s/templates_lj_%s_bkg_%s.root' % (outdir, njets, jobid), 'w') as bkg_out:
    mu_dir = bkg_out.mkdir('mujets')
    mu_dir.cd()

    for hname in mu_bkg:
        hist = mu_tfile.Get('mujets/%s' % hname)
        hist.Write()

    el_dir = bkg_out.mkdir('ejets')
    el_dir.cd()

    for hname in el_bkg:
        hist = el_tfile.Get('ejets/%s' % hname)
        hist.Write()


with io.root_open('%s/templates_lj_%s_sig_%s.root' % (outdir, njets, jobid), 'w') as sig_out:
    mu_dir = sig_out.mkdir('mujets')
    mu_dir.cd()

    for hname in mu_sig:
        hist = mu_tfile.Get('mujets/%s' % hname)
        hist.Write()

    el_dir = sig_out.mkdir('ejets')
    el_dir.cd()

    for hname in el_sig:
        hist = el_tfile.Get('ejets/%s' % hname)
        hist.Write()
