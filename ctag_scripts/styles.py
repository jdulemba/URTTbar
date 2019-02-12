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
      'fillcolor' : '#FFD700',
      #'fillcolor' : '#FFD700',
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "V + jets",
      'fillstyle': 'solid',
      },
   #'W[1-5]Jets*' : {
   'W[1-4]Jets*' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : '#FFD700',
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
   ##   'fillcolor' : '#FFD700',
   ##   'linecolor' : '#FFD700',
   ##   'name' : "W + jets",
   ##   'fillstyle': 'solid',
   ##   },
   'ZJets*' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : '#FFD700',
      'linecolor' : 'black',
      'name' : "V + jets",
      'fillstyle': 'solid',
      },
   'single*' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : '#008900',
      #'fillcolor' : '#2aa198',
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "single top",
      'fillstyle': 'solid',
      },
   'data*B*' : {
      'legendstyle' : 'p',
      'drawstyle' : 'E0 X0',
      'markerstyle' : 20,
      #'markersize'  : 2,
      'name' : "Run B",
    },
   'data*CtoE*' : {
      'legendstyle' : 'p',
      'drawstyle' : 'E0 X0',
      'markerstyle' : 20,
      #'markersize'  : 2,
      'name' : "Runs C-E",
    },
   'data*EtoF*' : {
      'legendstyle' : 'p',
      'drawstyle' : 'E0 X0',
      'markerstyle' : 20,
      #'markersize'  : 2,
      'name' : "Runs E-F",
    },
   #'data*BtoF*' : {
   'data*' : {
      'legendstyle' : 'p',
      'drawstyle' : 'E0 X0',
      'markerstyle' : 20,
      #'markersize'  : 2,
      'name' : "2016 Legacy",
      #'name' : "Observed",
    },
   'tt*' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : 'r',
      #'fillcolor' : '#9999CC',
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "t#bar{t}",
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
      'fillcolor' : 'c',
      #'fillcolor' : '#0055ff',
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "QCD",
      'fillstyle': 'solid',
      },
   'wrong_whad *' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : 'darkblue',
      #'fillcolor' : '#ab5555',
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "tt, wrong W_{h}",
      'fillstyle': 'solid',
      },
   'right_whad *' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : 'r',
      #'fillcolor' : '#6666b3',
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "tt, right W_{h}",
      'fillstyle': 'solid',
      },
   'nonsemi_tt *' : {
      'legendstyle' : 'f',
      'drawstyle' : 'hist',
      'fillcolor' : 'm',
      #'fillcolor' : '#668db3',
      'linecolor' : 'black',
			'linewidth' : 1,
      'name' : "tt, non semi-lep.",
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
