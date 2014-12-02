from numpy import *
from matplotlib.pyplot import *
import argparse
import os.path
from scipy.optimize import curve_fit



filenames = [ "2014-11-8_20-17", "2014-11-8_20-18", "2014-11-8_20-20", 
              "2014-11-8_20-23", "2014-11-8_20-26", "2014-11-8_20-30", 
              "2014-11-8_20-35", "2014-11-8_20-39"                     ]
cpufreq   = genfromtxt( "./"+filenames[0]+"_stats.txt", comments='#' )[0]
ardata    = []
for i in range(len(filenames)):
    ardata.append( genfromtxt( "./"+filenames[i]+"_Memcpy_Caching.txt", comments='#' ) )

ibytes  = 0 
iclocks = 1
icerr   = 2
stride  = 5

legendfontsize = 12
cachesizes = [ 32*1024, 256*1024, 20*1024**2 ]

def latencyFunc( NBytes, tLatency, bandwidth):
    return NBytes / (tLatency + NBytes/bandwidth)

def InitOnePlot(width=10):
    fig = figure( figsize=(width,4) )
    ax = subplot(111, xscale='log')
    ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
    xlabel("Message Length / Bytes")
    ylabel("Transfer Speed / (GByte/s)")
    for x in cachesizes:
        axvline( x,color='k',ls='dashed')
    return fig, ax


fig, ax = InitOnePlot()
title("All Threads doing memcpy simultaneously")
for i in range(len(ardata)):
    data = ardata[i][i*stride/len(ardata)::stride]
    # begin with element i, so that not alle points lie other each other
    y    = data[:,ibytes] / data[:,iclocks] * cpufreq / 1e9
    # relative Fehler pflanzt sich relativ fort, wenn y=LOP(x)
    yerr = data[:,icerr]  / data[:,iclocks] * y 
    #errorbar( data[:,ibytes], y, yerr=yerr, label=str(1+i)+" Threads")
    plot( data[:,ibytes], y, '.', label=str(1+i)+" Threads")
legend(prop={'size':legendfontsize}, loc='best' )
tight_layout()
fig.savefig( 'memcpyAll.pdf', format='PDF')




filenames = [ "2014-11-8_20-52", "2014-11-8_20-53", "2014-11-8_20-55", 
              "2014-11-8_20-56", "2014-11-8_20-58", "2014-11-8_21-0" , 
              "2014-11-8_21-2" , "2014-11-8_21-4"                     ]
ardata    = []
for i in range(len(filenames)):
    ardata.append( genfromtxt( "./"+filenames[i]+"_Memcpy_Caching.txt", comments='#' ) )
    
fig, ax = InitOnePlot()
title("Only first Thread doing memcpy, rest waiting")
for i in range(len(ardata)):
    data = ardata[i][i*stride/len(ardata)::stride]
    # begin with element i, so that not alle points lie other each other
    y    = data[:,ibytes] / data[:,iclocks] * cpufreq / 1e9
    # relative Fehler pflanzt sich relativ fort, wenn y=LOP(x)
    yerr = data[:,icerr]  / data[:,iclocks] * y 
    #errorbar( data[:,ibytes], y, yerr=yerr, label=str(1+i)+" Threads")
    plot( data[:,ibytes], y, '.', label=str(1+i)+" Threads")
legend(prop={'size':legendfontsize}, loc='best' )
tight_layout()
fig.savefig( 'memcpyOnlyOne.pdf', format='PDF')




show()
exit()
