#! /bin/env bash

pushd $URA_PROJECT/plots/$jobid/ctageff/mass_discriminant/ctagMedium/toys_simple/

tag='defaulted'
for src in $(ls toys_src.*.0.5.1.5.root); do
    seed=`echo $src | sed 's|toys_src.||g' | sed 's|.0.5.1.5.root||g'`
    combine ../fitModel.root -M MaxLikelihoodFit --saveToys --saveNLL --skipBOnlyFit -t 200 --minos=all --toysFile $src --name _$tag'_'$seed &> /dev/null 
done

##python <<EOF
##from rootpy.io import root_open
##from rootpy.plotting import Canvas, Hist
##from glob import glob
##import os
##from rootpy import asrootpy
##from rootpy import ROOT
##ROOT.gStyle.SetOptStat(0)
##
##csf = Hist(200, 0, 2)
##csf.markercolor = 'red'
##lsf = Hist(200, 0, 2)
##lsf.markercolor = 'blue'
##
##for name in glob('mlfit.toys_src.*0.5.1.5.root'):
##  with root_open(name) as mlfit:
##    toys = [i.name for i in mlfit.keys() if i.name.startswith('toy_')]
##    for i in toys:
##      path = os.path.join(i, 'fit_s')
##      toy = mlfit.Get(path)
##      pars = asrootpy(toy.floatParsFinal())
##      csf.Fill(float(pars['charmSF'].value))
##      lsf.Fill(float(pars['lightSF'].value))
##      toy.IsA().Destructor(toy)
##
##can = Canvas()
##can.SetLogy()
##pdir = '/uscms/home/verzetti/public_html/2016Feb23/ctageff_2016Apr05_newsf/'
##print max(csf.max(), lsf.max())
##csf.yaxis.range_user = (0.1, max(csf.max(), lsf.max()))
##csf.Draw()
##lsf.Draw('same')
##can.SaveAs('%s/sf_nosys.png' % pdir)     
##
##EOF

popd