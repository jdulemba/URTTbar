# Histogram styles:
# put in styles the stlye you want to apply to each sample you are processing
# in format key : {style_dict}, according to rootpy StyleView https://github.com/rootpy/rootpy/blob/0.7.1/rootpy/plotting/views.py
# prescriptions, any key named title will be strupped off for you and applied in a
# separate TitleView
#
# The keys are allowed to contain POSIX-style regular expressions, if multiple
# matches are found the longest key is used.

import ROOT
styles = {
   '[WZ]Jets*' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : '#FFCC66',
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "V + jets",
      'fillstyle': 'solid',
      },
   '[WZ][WZ]' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : '#f9a505',
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "diboson",
      'fillstyle': 'solid',
      },
   ##'WJets*' : {
   ##   'legendstyle' : 'f',
   ##   'drawstyle' : 'hist',
   ##   'fillcolor' : '#FFCC66',
   ##   'linecolor' : '#FFCC66',
   ##   'name' : "W + jets",
   ##   'fillstyle': 'solid',
   ##   },
   ##'ZJets*' : {
   ##   'legendstyle' : 'f',
   ##   'drawstyle' : 'hist',
   ##   'fillcolor' : '#FFCC66',
   ##   'linecolor' : '#FFCC66',
   ##   'name' : "Z + jets",
   ##   'fillstyle': 'solid',
   ##   },
   'single*' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : ROOT.kMagenta,
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "single top",
      'fillstyle': 'solid',
      },
   'data*' : {
      'legendstyle' : 'p',
      'drawstyle' : 'E0 X0',
      'markerstyle' : 20,
      #'markersize'  : 2,
      'name' : "Observed",
    },
   'tt*' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : ROOT.kOrange + 1,
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "ttbar",
      'fillstyle': 'solid',
      },
   'tt[WZ]*' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : ROOT.kOrange + 1,
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "ttV",
      'fillstyle': 'solid',
      },
   'QCD*' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : ROOT.kGray,
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "QCD",
      'fillstyle': 'solid',
      },
   }

from itertools import product
for bundle in product([400, 500, 600, 750], [5, 10, 25, 50]):
	styles['HtoTT_M%d_%dpc_*' % bundle] = {
		'legendstyle' : 'l',
		'drawstyle' : 'hist',
		'linecolor' : '#2fd00a',
		'linewidth' : 3,
		'name' : "A #rightarrow tt M%d width: %d%%" % bundle,
		'fillstyle': 'hollow',
		}
