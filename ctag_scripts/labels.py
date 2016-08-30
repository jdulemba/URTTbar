
labels = {
   'thadpt'  : 'p_{T}(t_{had}) [GeV]',
   'tleppt'  : 'p_{T}(t_{lep}) [GeV]',
   'thadeta' : '|#eta(t_{had})|',
   'tlepeta' : '|#eta(t_{lep})|',
   'tlepy'   : '|y(t_{lep})|',
   'thady'   : '|y(t_{had})|',
   'tty'     : '|y(tt)|',
   'ttpt'    : 'p_{T}(tt) [GeV]',
   'ttM'     : 'M(tt) [GeV]',
}
def set_pretty_label(variable):
   return labels.get(variable, variable)
