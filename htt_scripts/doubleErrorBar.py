# Collection of methods to plot statistical and systematic error bars on the same plot
# Each method takes x data, y data, y statistical error and y systematic error
# Each method also takes optional arguemts to constrol marker style and colors
# The 'inner' error bar is the statistical error
# The 'outer' error bar shows that statistcal+systematic error
# Recommended to copy function you want into your code and modify colors there
# Author: Rhys Taus March 15, 2018
#
# Optional Arguements
# ms = Marker style: Default black circles 'ko'
# color1 = color for stats error: Default blue 'b'
# color2 = color for sys error: Default red 'r'
#
# Optional Arguements to be added in the future
# msize = marker size
#
#
# Will develop to have options for x error bars
# Current recommendation is to have plot anoter plot with just the data to have
# a legend. 
import numpy as np
import matplotlib.pyplot as plt



def method1(x,y,y_stat,y_sys, **kwargs):
    # Simple 2 error bars
    if bool(kwargs.get('ms', None)):
        ms=kwargs.get('ms', None)
    else:
        ms='ko'
    #
    if bool(kwargs.get('color1', None)):
        c1=kwargs.get('color1', None)
    else:
        c1='b'
    #
    if bool(kwargs.get('color2', None)):
        c2=kwargs.get('color2', None)
    else:
        c2='r'
    #
    plt.errorbar(x, y, yerr=y_stat, fmt=ms, ecolor=c1)
    plt.errorbar(x, y, yerr=y_stat+y_sys, fmt=ms, ecolor=c2)


    
def method2(x,y,y_stat,y_sys, **kwargs):
    # Systematic errors are thick
    if bool(kwargs.get('ms', None)):
        ms=kwargs.get('ms', None)
    else:
        ms='ko'
    #
    if bool(kwargs.get('color1', None)):
        c1=kwargs.get('color1', None)
    else:
        c1='b'
    #
    if bool(kwargs.get('color2', None)):
        c2=kwargs.get('color2', None)
    else:
        c2='r'
    #
    plt.errorbar(x, y, yerr=y_stat, fmt=ms, ecolor=c1, capsize=4)
    plt.errorbar(x, y, yerr=y_stat+y_sys, fmt=ms, ecolor=c2, elinewidth=10, alpha=0.6)


    
def method3(x,y,y_stat,y_sys,**kwargs):
    # Both error bars are thick
    if bool(kwargs.get('ms', None)):
        ms=kwargs.get('ms', None)
    else:
        ms='ko'
    #
    if bool(kwargs.get('color1', None)):
        c1=kwargs.get('color1', None)
    else:
        c1='b'
    #
    if bool(kwargs.get('color2', None)):
        c2=kwargs.get('color2', None)
    else:
        c2='r'
    #
    plt.errorbar(x, y, yerr=y_stat, fmt=ms, ecolor=c1, elinewidth=10)
    plt.errorbar(x, y, yerr=y_stat+y_sys, fmt=ms, ecolor=c2, elinewidth=10, alpha=0.6)

##########################
# Test different methods #
##########################

x = [k/10. for k in range(21)]
y = np.exp(x)
y_stat = y/10.
y_sys = y/20.

# Default
plt.figure()
method2(x,y,y_stat,y_sys)
plt.show()


#with arguemets
plt.figure()
method3(x,y,y_stat,y_sys, ms='kd', color1='c', color2='m')
plt.legend()
plt.show()

