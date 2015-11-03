
labels = {
   'thadpt'  : 'p_{T}(t_{had}) [GeV]',
   'tleppt'  : 'p_{T}(t_{lep}) [GeV]',
   'thadeta' : '|#eta(t_{had})|',
   'tlepeta' : '|#eta(t_{lep})|',
}
def set_pretty_label(variable):
   return labels.get(variable, variable)
