from numpy import *
from matplotlib.pyplot import *
import argparse
import os.path
from scipy.optimize import curve_fit

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Directory which contains text files from local benchmark")
#parser.add_argument("fileb", help="Directory which contains text files from non-local benchmark")
parser.add_argument("-t", help="Show statistics for measurements for same length")
args = parser.parse_args()

#python plot.py Taurus_2014-10-30_18-45/2014-10-30_18-45 Taurus_2014-10-30_21-0/2014-10-30_21-0

if args.i == None:
    args.file   = "Taurus_2014-10-30_18-45/2014-10-30_18-45" # on same node 
    args.fileb  = "Taurus_2014-10-30_21-0/2014-10-30_21-0"   # two nodes
    args.filec  = "Home/2014-10-30_18-40"                    # on same node at home
    #args.filec  = "Home/2014-10-25_20-3"                    # at home (older results)

stride      = 1  # if the plot is too crowded increase this
cpufreq     = genfromtxt( "./"+args.file +"_stats.txt", comments='#' )[0]
cpufreqb    = genfromtxt( "./"+args.fileb+"_stats.txt", comments='#' )[0]
cpufreqc    = genfromtxt( "./"+args.filec+"_stats.txt", comments='#' )[0]
cpufreqerr  = genfromtxt( "./"+args.file +"_stats.txt", comments='#' )[1]
cpufreqerrb = genfromtxt( "./"+args.fileb+"_stats.txt", comments='#' )[1]
cpufreqerrc = genfromtxt( "./"+args.filec+"_stats.txt", comments='#' )[1]

print "CPU Frequency A:", cpufreq /1e9, "+-", cpufreqerr /1e9, "GHz"
print "CPU Frequency B:", cpufreqb/1e9, "+-", cpufreqerrb/1e9, "GHz"
print "CPU Frequency C:", cpufreqc/1e9, "+-", cpufreqerrc/1e9, "GHz"

fname11 = "./"+args.file+"_Caching"
fname12 = "./"+args.file+"_Caching_LongestFirst"
fname13 = "./"+args.file+"_NoCache_LongestFirst"
fname14 = "./"+args.file+"_Caching_Cumulative"
fname21 = "./"+args.file+"_Memcpy_Caching"
fname22 = "./"+args.file+"_Memcpy_Caching_LongestFirst"
fname23 = "./"+args.file+"_Memcpy_NoCache_LongestFirst"
fname11b = "./"+args.fileb+"_Caching"
fname12b = "./"+args.fileb+"_Caching_LongestFirst"
fname13b = "./"+args.fileb+"_NoCache_LongestFirst"
fname14b = "./"+args.fileb+"_Caching_Cumulative"
fname21b = "./"+args.fileb+"_Memcpy_Caching"
fname22b = "./"+args.fileb+"_Memcpy_Caching_LongestFirst"
fname23b = "./"+args.fileb+"_Memcpy_NoCache_LongestFirst"
fname11c = "./"+args.filec+"_Caching"
fname12c = "./"+args.filec+"_Caching_LongestFirst"
fname13c = "./"+args.filec+"_NoCache_LongestFirst"
fname14c = "./"+args.filec+"_Caching_Cumulative"
fname21c = "./"+args.filec+"_Memcpy_Caching"
fname22c = "./"+args.filec+"_Memcpy_Caching_LongestFirst"
fname23c = "./"+args.filec+"_Memcpy_NoCache_LongestFirst"

if args.t:
    statdata11  = genfromtxt( fname11+"_Times.txt", comments='#' )
    dataset = 0
    def hist_key_analyzer(event):
        global dataset
        print event.key
        print "Dataset:",dataset
        if (event.key == '+' or event.key == "right"):
            dataset = (dataset+1) % len(statdata11[:,0])
        if (event.key == '-' or event.key == "left"):
            dataset = (dataset-1) % len(statdata11[:,0])
        figure( hist_fig.number )
        clf()
        suptitle( str(statdata11[dataset,0]) + " Bytes" )
        hist( statdata11[dataset,1:], bins=50 )
        draw()

    hist_fig = figure( figsize=(6,4) )
    connect('key_press_event', hist_key_analyzer)
    suptitle( str(statdata11[dataset,0]) + " Bytes" )
    hist( statdata11[dataset,1:], bins=50 )
    tight_layout()


data11  = genfromtxt( fname11+".txt", comments='#' )
data12  = genfromtxt( fname12+".txt", comments='#' )
data13  = genfromtxt( fname13+".txt", comments='#' )
data14  = genfromtxt( fname14+".txt", comments='#' )
data21  = genfromtxt( fname21+".txt", comments='#' )
data22  = genfromtxt( fname22+".txt", comments='#' )
data23  = genfromtxt( fname23+".txt", comments='#' )

data11b  = genfromtxt( fname11b+".txt", comments='#' )
data12b  = genfromtxt( fname12b+".txt", comments='#' )
data13b  = genfromtxt( fname13b+".txt", comments='#' )
data14b  = genfromtxt( fname14b+".txt", comments='#' )
data21b  = genfromtxt( fname21b+".txt", comments='#' )
data22b  = genfromtxt( fname22b+".txt", comments='#' )
data23b  = genfromtxt( fname23b+".txt", comments='#' )

data11c  = genfromtxt( fname11c+".txt", comments='#' )
data12c  = genfromtxt( fname12c+".txt", comments='#' )
data13c  = genfromtxt( fname13c+".txt", comments='#' )
data14c  = genfromtxt( fname14c+".txt", comments='#' )
data21c  = genfromtxt( fname21c+".txt", comments='#' )
data22c  = genfromtxt( fname22c+".txt", comments='#' )
data23c  = genfromtxt( fname23c+".txt", comments='#' )

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


###################### Show some single measurements ###########################

from matplotlib.lines import Line2D
data11times = genfromtxt( fname11+"_Times.txt", comments='#' )
data13times = genfromtxt( fname13+"_Times.txt", comments='#' )
fig = figure( figsize=(10,4) )
ax  = subplot(121, xlabel="Measurement", ylabel="CPU cycles")
plot( data11times[ 0,1:], 'bs', label=str(data11times[0,0])+" Bytes" )
plot( data11times[ 1,1:], 'bo', label=str(data11times[1,0])+" Bytes" )
plot( data11times[ 2,1:], 'b^', label=str(data11times[2,0])+" Bytes" )
plot( data13times[-1,1:], 'rs' )
plot( data13times[-2,1:], 'ro' )
plot( data13times[-3,1:], 'r^' )
legend(prop={'size':legendfontsize})
xlim( (-1,len(data11times[ 0,1:])) )
xmin, xmax, ymin, ymax = axis()
axR = ax.twinx()
ylim( (ymin/cpufreq*1e6, ymax/cpufreq*1e6) )
ylabel( r"Time / $\mu$s" )
ax  = subplot(122, xlabel="Measurement", ylabel="CPU cycles" )
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
plot( data11times[-2,1:], 'bs', label=str(data11times[-2,0])+" Bytes" )
plot( data11times[-3,1:], 'bo', label=str(data11times[-3,0])+" Bytes" )
plot( data11times[-4,1:], 'b^', label=str(data11times[-4,0])+" Bytes" )
plot( data13times[ 0,1:], 'rs' )
plot( data13times[ 1,1:], 'ro' )
plot( data13times[ 2,1:], 'r^' )
legend(prop={'size':legendfontsize})
xlim( (-1,len(data11times[ 0,1:])) )
xmin, xmax, ymin, ymax = axis()
axR = ax.twinx()
ylim( (ymin/cpufreq*1e6, ymax/cpufreq*1e6) )
ylabel( r"Time / $\mu$s" )
fig.suptitle( "Local Benchmark" )
#tight_layout()
subplots_adjust( left=0.08, bottom=0.12,  right=0.93, wspace=0.34 )
fig.savefig( 'SingleTimes.pdf', format='PDF')


########################### Two Cores on One Node ##############################

stride  = 10
fig, ax = InitOnePlot()
ardata  = [data13,data11,data12]
arlabel = ["Random Data","Shortest First","Longest First"]
arfmt   = ['mo','bo','ko']
for i in range(len(ardata)):
    data = ardata[i][i*stride/len(ardata)::stride]
    # begin with element i, so that not alle points lie other each other
    y    = data[:,ibytes] / data[:,iclocks] * cpufreq / 1e9
    # relative Fehler pflanzt sich relativ fort, wenn y=LOP(x)
    yerr = data[:,icerr]  / data[:,iclocks] * y 
    errorbar( data[:,ibytes], y, yerr=yerr, fmt=arfmt[i], label=arlabel[i] )
plot( data14[::stride,ibytes], data14[::stride,ibytes]/data14[::stride,1] / 1e9, 'go', label="MPI_Wtime")
legend(prop={'size':legendfontsize}, loc='best' )
tight_layout()
fig.savefig( 'SameNode.pdf', format='PDF')

stride  = 8
xlim( (130000, 6.5e6) )
ylim( (7, 13) )
tight_layout()
fig.savefig( 'SameNode_Zoom.pdf', format='PDF')


################################# At Home ######################################

stride  = 10
fig, ax   = InitOnePlot()
ardata  = [data13c,data11c,data12c] # file c is home
for i in range(len(ardata)):
    data = ardata[i][i*stride/len(ardata)::stride]
    y    = data[:,ibytes] / data[:,iclocks] * cpufreq / 1e9
    yerr = data[:,icerr]  / data[:,iclocks] * y 
    errorbar( data[:,ibytes], y, yerr=yerr, fmt=arfmt[i], label=arlabel[i] )
stride=1
plot( data14c[::stride,ibytes], data14c[::stride,ibytes]/data14c[::stride,1] / 1e9, 'g.', label="MPI_Wtime")
legend(prop={'size':legendfontsize}, loc='best' )
tight_layout()
fig.savefig( 'Home.pdf', format='PDF')

fig, ax   = InitOnePlot()
plot( data21c[::stride,0], data21c[::stride,0]/data21c[::stride,1] * cpufreqc/1e9, 'b.', label="Shortest First" )
plot( data22c[::stride,0], data22c[::stride,0]/data22c[::stride,1] * cpufreqc/1e9, 'k.', label="Longest First"  ) 
legend(prop={'size':legendfontsize}, loc='upper left' )
tight_layout()
fig.savefig( 'Home_Memcpy.pdf', format='PDF')


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
ardata   = [data13b,data11b,data12b] # file b is remote
stride   = 10
sca(axL)
for i in range(len(ardata)):
    data = ardata[i][i*stride/len(ardata)::stride]
    y    = data[:,ibytes] / data[:,iclocks] * cpufreq / 1e9
    yerr = data[:,icerr]  / data[:,iclocks] * y 
    errorbar( data[:,ibytes], y, yerr=yerr, fmt=arfmt[i], label=arlabel[i] )
plot( data14b[::stride,ibytes], data14b[::stride,ibytes]/data14b[::stride,1] / 1e9, 'g.', label="MPI_Wtime")
sca(axR)
xlim( (1e3, 4e3) )
ylim( (0, 5.5) )
ardata  = [data13b,data11b,data12b] # file b is remote
for i in range(len(ardata)):
    data = ardata[i][i*stride/len(ardata)::stride]
    y    = data[:,ibytes] / data[:,iclocks] * cpufreq / 1e9
    yerr = data[:,icerr]  / data[:,iclocks] * y 
    errorbar( data[:,ibytes], y, yerr=yerr, fmt=arfmt[i], label=arlabel[i] )
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
print "Fit in ",args.file,fname14b,":"
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
mask  = data22[:,0] < 1.2e4
xdata = data22[:,0][mask]
ydata = data22[:,0][mask]/data22[:,1][mask] * cpufreq
popt, pcov = curve_fit(latencyFunc, xdata, ydata, p0=[1e-6,40e8] )
perr  = np.sqrt(np.diag(pcov))
print "Fit in ",args.file,fname22,":"
print "    tLatency  =",popt[0]," +- ",perr[0]
print "    bandwidth =",popt[1]," +- ",perr[1]

axL, axR = InitTwoPlots()
sca(axL)
plot( data21[::stride,0], data21[::stride,0]/data21[::stride,1] * cpufreq/1e9, 'b.', label="Shortest First" )
plot( data22[::stride,0], data22[::stride,0]/data22[::stride,1] * cpufreq/1e9, 'k.', label="Longest First"  ) 

sca(axR)
axR.set_xscale('linear')
axR.xaxis.get_major_formatter().set_powerlimits((0, 2))
xlim( (0,3e4) )

plot( data21[::stride,ibytes], data21[::stride,ibytes]/data21[::stride,iclocks] * cpufreq/1e9, 'b.', label="Shortest First" )
plot( data22[::stride,ibytes], data22[::stride,ibytes]/data22[::stride,iclocks] * cpufreq/1e9, 'k.', label="Longest First"  )  
plot( xdata, latencyFunc( xdata, popt[0], popt[1] )/1e9, 'r-', label=r"Fit" )
legend(prop={'size':legendfontsize}, loc='lower right' )
tight_layout()
fig.savefig( 'SameNode_Memcpy_NoRandomData.pdf', format='PDF')


fig,ax = InitOnePlot()
plot( data21[::stride,0], data21[::stride,0]/data21[::stride,1] * cpufreq/1e9, 'b.', label="Shortest First" )
plot( data22[::stride,0], data22[::stride,0]/data22[::stride,1] * cpufreq/1e9, 'k.', label="Longest First"  ) 
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
plot( data22[::stride,0], data22[::stride,0]/data22[::stride,1] * cpufreq/1e9, 'k.', label="Longest First"   )
plot( data23[::stride,0], data23[::stride,0]/data23[::stride,1] * cpufreq/1e9, 'g.', label="Prevent Caching" )
legend(prop={'size':legendfontsize}, loc='upper left' )
ax = subplot(122, xscale='log')
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
xlabel("Message Length / Bytes")
ylabel("Transfer Speed / (GByte/s)")
title('Other Node')
plot( data21b[::stride,0], data21b[::stride,0]/data21b[::stride,1] * cpufreqb/1e9, 'b.', label="Shortest First"  )
plot( data22b[::stride,0], data22b[::stride,0]/data22b[::stride,1] * cpufreqb/1e9, 'k.', label="Longest First"   )
plot( data23b[::stride,0], data23b[::stride,0]/data23b[::stride,1] * cpufreqb/1e9, 'g.', label="Prevent Caching" )
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





fig = figure(figsize=(12,10))

ax = subplot(231)
xlabel("Message Length / Bytes")
ylabel("Time / s")
ax.xaxis.get_major_formatter().set_powerlimits((0, 2))
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
errorbar( data11[::stride,0], data11[::stride,1] / cpufreq, fmt='b.', yerr=data11[::stride,2] / cpufreq, label="Begin With Shortest Message"    )
errorbar( data12[::stride,0], data12[::stride,1] / cpufreq, fmt='k.', yerr=data12[::stride,2] / cpufreq, label="Longest First"                  )
errorbar( data13[::stride,0], data13[::stride,1] / cpufreq, fmt='g.', yerr=data13[::stride,2] / cpufreq, label="Prevent Caching, Longest First" )
errorbar( data14[::stride,0], data14[::stride,1]          , fmt='r.',                   label="Cumulative Time Meas."          )

subplot(232, xscale='log', yscale='log')
xlabel("Message Length / Bytes")
ylabel("Time / s")
errorbar( data11[::stride,0], data11[::stride,1] / cpufreq, fmt='b.', yerr=data11[::stride,2] / cpufreq, label="Begin With Shortest Message"    )
errorbar( data12[::stride,0], data12[::stride,1] / cpufreq, fmt='k.', yerr=data12[::stride,2] / cpufreq, label="Longest First"                  )
#errorbar( data13[::stride,0], data13[::stride,1] / cpufreq, fmt='g.', yerr=data13[::stride,2] / cpufreq, label="Prevent Caching, Longest First" )
#errorbar( data14[::stride,0], data14[::stride,1]          , fmt='r.',                   label="Cumulative Time Meas."          )

ax = subplot(233, xscale='log')
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
xlabel("Message Length / Bytes")
ylabel("Transfer Speed / (GByte/s)")
plot( data11[::stride,0], data11[::stride,0]/data11[::stride,1] * cpufreq/1e9, 'b.', label="Begin With Shortest Message"    )
plot( data12[::stride,0], data12[::stride,0]/data12[::stride,1] * cpufreq/1e9, 'k.', label="Longest First"                  )
#plot( data13[::stride,0], data13[::stride,0]/data13[::stride,1] * cpufreq/1e9, 'g.', label="Prevent Caching, Longest First" )
#plot( data14[::stride,0], data14[::stride,0]/data14[::stride,1]          /1e9, 'r.', label="Cumulative Time Meas."          )
legend(prop={'size':legendfontsize},)



ax = subplot(234)
xlabel("Message Length / Bytes")
ylabel("Time / s")
ax.xaxis.get_major_formatter().set_powerlimits((0, 2))
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
errorbar( data21[::stride,0], data21[::stride,1] / cpufreq, fmt='b.', yerr=data21[::stride,2] / cpufreq, label="Begin With Shortest Message"    )
errorbar( data22[::stride,0], data22[::stride,1] / cpufreq, fmt='k.', yerr=data22[::stride,2] / cpufreq, label="Longest First"                  )
#errorbar( data23[::stride,0], data23[::stride,1] / cpufreq, fmt='g.', yerr=data23[::stride,2] / cpufreq, label="Prevent Caching, Longest First" )

subplot(235, xscale='log', yscale='log')
xlabel("Message Length / Bytes")
ylabel("Time / s")
errorbar( data21[::stride,0], data21[::stride,1] / cpufreq, fmt='b.', yerr=data21[::stride,2] / cpufreq, label="Begin With Shortest Message"    )
errorbar( data22[::stride,0], data22[::stride,1] / cpufreq, fmt='k.', yerr=data22[::stride,2] / cpufreq, label="Longest First"                  )
#errorbar( data23[::stride,0], data23[::stride,1] / cpufreq, fmt='g.', yerr=data23[::stride,2] / cpufreq, label="Prevent Caching, Longest First" )

ax = subplot(236, xscale='log')
ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
xlabel("Message Length / Bytes")
ylabel("Transfer Speed / (GByte/s)")
plot( data21[::stride,0], data21[::stride,0]/data21[::stride,1] * cpufreq/1e9, 'b.', label="Begin With Shortest Message"    )
plot( data22[::stride,0], data22[::stride,0]/data22[::stride,1] * cpufreq/1e9, 'k.', label="Longest First"                  )
#plot( data23[::stride,0], data23[::stride,0]/data23[::stride,1] * cpufreq/1e9, 'g.', label="Prevent Caching, Longest First" )
legend(prop={'size':legendfontsize},)

tight_layout()

show()