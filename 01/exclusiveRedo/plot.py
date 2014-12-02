from numpy import *
from matplotlib.pyplot import *
import argparse
import os.path
from scipy.optimize import curve_fit



file   = "2014-11-8_21-58" # on same node 
fileb  = "2014-11-8_22-12"   # two nodes
    
stride      = 1  # if the plot is too crowded increase this
cpufreq     = genfromtxt( "./"+file +"_stats.txt", comments='#' )[0]
cpufreqb    = genfromtxt( "./"+fileb+"_stats.txt", comments='#' )[0]
cpufreqerr  = genfromtxt( "./"+file +"_stats.txt", comments='#' )[1]
cpufreqerrb = genfromtxt( "./"+fileb+"_stats.txt", comments='#' )[1]

print "CPU Frequency A:", cpufreq /1e9, "+-", cpufreqerr /1e9, "GHz"
print "CPU Frequency B:", cpufreqb/1e9, "+-", cpufreqerrb/1e9, "GHz"

fname11  = "./"+file+"_Caching"
fname14  = "./"+file+"_Caching_Cumulative"
fname21  = "./"+file+"_Memcpy_Caching"
fname11b = "./"+fileb+"_Caching"
fname14b = "./"+fileb+"_Caching_Cumulative"
fname21b = "./"+fileb+"_Memcpy_Caching"

data11  = genfromtxt( fname11+".txt", comments='#' )
data14  = genfromtxt( fname14+".txt", comments='#' )
data21  = genfromtxt( fname21+".txt", comments='#' )
data11b  = genfromtxt( fname11b+".txt", comments='#' )
data14b  = genfromtxt( fname14b+".txt", comments='#' )
data21b  = genfromtxt( fname21b+".txt", comments='#' )

ibytes  = 0 
iclocks = 1
icerr   = 2
iwtime  = 3

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



########################### Two Cores on One Node ##############################

stride  = 10
fig, ax = InitOnePlot()
data = data11[::stride]
# begin with element i, so that not alle points lie other each other
y    = data[:,ibytes] / data[:,iclocks] * cpufreq / 1e9
# relative Fehler pflanzt sich relativ fort, wenn y=LOP(x)
yerr = data[:,icerr]  / data[:,iclocks] * y 
errorbar( data[:,ibytes], y, yerr=yerr, label="GetCPUCycle" )
plot( data14[::stride,ibytes], data14[::stride,ibytes]/data14[::stride,1] / 1e9, 'go', label="MPI_Wtime")
legend(prop={'size':legendfontsize}, loc='best' )
tight_layout()
fig.savefig( 'SameNode.pdf', format='PDF')

stride  = 8
xlim( (130000, 6.5e6) )
ylim( (7, 13) )
tight_layout()
fig.savefig( 'SameNode_Zoom.pdf', format='PDF')


####################### Two Cores on Different Nodes ###########################

def InitTwoPlots():
    fig = figure( figsize=(10,4) )
    axL = subplot(121, xscale='log', xlabel="Message Length / Bytes", 
                                     ylabel="Transfer Speed / (GByte/s)")
    axL.yaxis.get_major_formatter().set_powerlimits((0, 2))
    axR = subplot(122, xscale='log', xlabel="Message Length / Bytes", 
                                     ylabel="Transfer Speed / (GByte/s)")
    axR.yaxis.get_major_formatter().set_powerlimits((0, 2))
    return axL,axR

axL, axR = InitTwoPlots()
stride   = 10
sca(axL)

data = data11b[::stride]
y    = data[:,ibytes] / data[:,iclocks] * cpufreq / 1e9
yerr = data[:,icerr]  / data[:,iclocks] * y 
errorbar( data[:,ibytes], y, yerr=yerr, label="GetCPUCycle" )
plot( data14b[::stride,ibytes], data14b[::stride,ibytes]/data14b[::stride,1] / 1e9, 'g.', label="MPI_Wtime")
sca(axR)
xlim( (1e3, 4e3) )
ylim( (0, 5.5) )

data = data11b[::stride]
y    = data[:,ibytes] / data[:,iclocks] * cpufreq / 1e9
yerr = data[:,icerr]  / data[:,iclocks] * y 
errorbar( data[:,ibytes], y, yerr=yerr, label="GetCPUCycle"  )
plot( data14b[::stride,ibytes], data14b[::stride,ibytes]/data14b[::stride,1] / 1e9, 'g.', label="MPI_Wtime")
legend(prop={'size':legendfontsize}, loc='best' )
tight_layout()
fig.savefig( 'DiffNodes.pdf', format='PDF')


################# Compare local vs. non-local communication ####################

stride=1
mask       = data14b[:,ibytes] > 12100
xdata      = data14b[:,ibytes][mask]
ydata      = data14b[:,ibytes][mask]/data14b[:,1][mask]
popt, pcov = curve_fit(latencyFunc, xdata, ydata, p0=[1e-6,10e9] )
perr  = np.sqrt(np.diag(pcov))
print "Fit in ",file,fname14b,":"
print "    tLatency  =",popt[0]," +- ",perr[0]
print "    bandwidth =",popt[1]," +- ",perr[1]
fig, ax    = InitOnePlot()
plot( data14[::stride,ibytes] , data14[::stride,ibytes] /data14[::stride,1]  / 1e9, 'k.', label="local")
plot( data14b[:,ibytes], data14b[:,ibytes]/data14b[:,1] / 1e9, 'b.', label="remote")

plot( data14b[:,ibytes], latencyFunc( data14b[:,ibytes], popt[0], popt[1] )/1e9, 'r-', label=r"Fit" )
legend(prop={'size':legendfontsize}, loc='best' )
tight_layout()
fig.savefig( 'LocalNonLocal.pdf', format='PDF')


############################### Memcpy on Node #################################

# Fit to Data
mask  = data21[:,0] < 1.2e4
xdata = data21[:,0][mask]
ydata = data21[:,0][mask]/data21[:,1][mask] * cpufreq
popt, pcov = curve_fit(latencyFunc, xdata, ydata, p0=[1e-6,40e8] )
perr  = np.sqrt(np.diag(pcov))
print "Fit in ",file,fname21,":"
print "    tLatency  =",popt[0]," +- ",perr[0]
print "    bandwidth =",popt[1]," +- ",perr[1]

axL, axR = InitTwoPlots()
sca(axL)
plot( data21[::stride,0], data21[::stride,0]/data21[::stride,1] * cpufreq/1e9, 'b.', label="Shortest First" )

sca(axR)
axR.set_xscale('linear')
axR.xaxis.get_major_formatter().set_powerlimits((0, 2))
xlim( (0,3e4) )

plot( data21[::stride,ibytes], data21[::stride,ibytes]/data21[::stride,iclocks] * cpufreq/1e9, 'b.', label="Shortest First" )
plot( xdata, latencyFunc( xdata, popt[0], popt[1] )/1e9, 'r-', label=r"Fit" )
legend(prop={'size':legendfontsize}, loc='lower right' )
tight_layout()
fig.savefig( 'SameNode_Memcpy_NoRandomData.pdf', format='PDF')


fig,ax = InitOnePlot()
plot( data21[::stride,0], data21[::stride,0]/data21[::stride,1] * cpufreq/1e9, 'b.', label="Shortest First" )
legend(prop={'size':legendfontsize}, loc='best' )

tight_layout()
fig.savefig( 'SameNode_Memcpy_ZoomBorderLinLog.pdf', format='PDF')

show()
exit()

fig = figure( figsize=(10,4) )
ax = subplot(121, xscale='log')
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
xlabel("Message Length / Bytes")
ylabel("Transfer Speed / (GByte/s)")
title('Same Node')
plot( data21[::stride,0], data21[::stride,0]/data21[::stride,1] * cpufreq/1e9, 'b.', label="Shortest First"  )
legend(prop={'size':legendfontsize}, loc='upper left' )
ax = subplot(122, xscale='log')
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
xlabel("Message Length / Bytes")
ylabel("Transfer Speed / (GByte/s)")
title('Other Node')
plot( data21b[::stride,0], data21b[::stride,0]/data21b[::stride,1] * cpufreqb/1e9, 'b.', label="Shortest First"  )
tight_layout()
fig.savefig( 'Memcpy_RandomDataComparison.pdf', format='PDF')


########################## Compare Memcpy vs. MPI_Send #########################

fig = figure( figsize=(12,4) )
ax = subplot(131, xscale='log')
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
xlabel("Message Length / Bytes")
ylabel("Transfer Speed / (GByte/s)")
title('Same Node')
plot( data21[::stride,0], data21[::stride,0]/data21[::stride,1] * cpufreq/1e9, 'k.', label="Memcpy"   )
plot( data11[::stride,0], data11[::stride,0]/data11[::stride,1] * cpufreq/1e9, 'r.', label="MPI_Send" )
legend(prop={'size':legendfontsize}, loc='upper left' )
ax = subplot(132, xscale='log')
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
xlabel("Message Length / Bytes")
ylabel("Transfer Speed / (GByte/s)")
title('Other Node')
plot( data21b[::stride,0], data21b[::stride,0]/data21b[::stride,1] * cpufreqb/1e9, 'k.', label="Memcpy"   )
plot( data11b[::stride,0], data11b[::stride,0]/data11b[::stride,1] * cpufreqb/1e9, 'r.', label="MPI_Send" )
legend(prop={'size':legendfontsize}, loc='upper left' )
ax = subplot(133, xscale='log')
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
xlabel("Message Length / Bytes")
ylabel("Transfer Speed / (GByte/s)")
title('Home')
plot( data21c[::stride,0], data21c[::stride,0]/data21c[::stride,1] * cpufreqc/1e9, 'k.', label="Memcpy"   )
plot( data11c[::stride,0], data11c[::stride,0]/data11c[::stride,1] * cpufreqc/1e9, 'r.', label="MPI_Send" )
legend(prop={'size':legendfontsize}, loc='upper left' )
tight_layout()
fig.savefig( 'Memcpy_MPISendComparison.pdf', format='PDF')



show()