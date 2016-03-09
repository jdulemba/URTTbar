yields = eval(open('ctag_eff.test.json').read())
old = eval(open('tests/ctag_eff.diffxsec.json').read())
for key, val in yields.iteritems():
   if key in old:                                                                                                                                                                                                
      oval = old[key]
      diff = '%.0f%%' % (100*(val-oval)/oval) if oval else '%f --> %f' % (oval, val)
      print '%s %s' % (key, diff)
   else:                                                                                                                                                                                                         
      print '%s --> NEW' % key                   
