from numpy import *
from matplotlib.pyplot import *

stride  = 1
ardata = []

ardata.append( genfromtxt( "./2014-11-9_20-44_data.txt", comments='#' ) )
ardata.append( genfromtxt( "./2014-11-9_20-52_data.txt", comments='#' ) )
ardata.append( genfromtxt( "./2014-11-10_18-51_data.txt", comments='#' ) )
ardata.append( genfromtxt( "./2014-11-10_19-13_data.txt", comments='#' ) )
ardata.append( genfromtxt( "./2014-11-10_20-4_data.txt", comments='#' ) )

ardata.append( genfromtxt( "./2014-10-30_23-45_data.txt", comments='#' ) )
ardata.append( genfromtxt( "./2014-11-9_19-41_data.txt", comments='#' ) )

ardata.append( genfromtxt( "./2014-11-11_20-46_data.txt", comments='#' ) )
ardata.append( genfromtxt( "./2014-11-11_21-48_data.txt", comments='#' ) )
arlabels = [ "On Taurusi 1004",
             "Home, i3 3320",
             "On Taurus 1189, 8 OpenMP Threads with AVX",
             "Home, i3 3320, 4 OpenMP Threads with AVX",
             "On Taurus 1189, 8 OpenMP Threads with AVX",
             "1 - old measurements home",
             "2 - old measurements home",
             "more points home",
             "more points Taurusi1254",
           ]
datamemcpy = genfromtxt( "../01/Taurus_2014-10-30_18-45/2014-10-30_18-45_Memcpy_Caching.txt", comments='#' )
cpufreqmemcpy = genfromtxt( "../01/Taurus_2014-10-30_18-45/2014-10-30_18-45_stats.txt", comments='#' )[0]
datamemcpynew = genfromtxt( "./2014-11-10_19-39_data.txt", comments='#' )

fig = figure(figsize=(8,6))

# integration intervals	processors used	time / CPU clocks	sigma CPU clocks	result	relative Error
imsize = 0
idsize = 1
it     = 2
iterr  = 3

print "Array of Datatypesize in Bytes (should be all the same!): "
for k in range(len(ardata)):
    print ardata[0][0,idsize]

ax = subplot(111, xscale='log', xlabel="Matrix Size / Bytes",
                  yscale='linear', ylabel="Gflops" )
cachesizes = [ 32*1024, 256*1024, 20*1024**2 ]
for x in cachesizes:
    axvline( x,color='k',ls='dashed')
for k in range(len(ardata)):
    y = 2*ardata[k][:,imsize]**3 / ardata[k][:,it]
    errorbar( ardata[k][:,imsize]**2*ardata[k][:,idsize], y/1e9, yerr=ardata[k][:,iterr]/ardata[k][:,it]*y/1e9, label=arlabels[k] )
stride = 5
plot( datamemcpy[::stride,0], datamemcpy[::stride,0]/datamemcpy[::stride,1] * cpufreqmemcpy/1e9 / ardata[0][0,idsize], 'b.', label="Memcpy Taurus (1 flop=8 Byte/s)" )

#size = datamemcpy[:,imsize]**2*datamemcpy[:,idsize]
#errorbar( size, size/datamemcpy[:,it]/1e9, yerr=datamemcpy[:,iterr]/datamemcpy[:,it] * size/datamemcpy[:,it]/1e9, fmt='b.', label="Memcpy Taurus (1 flop=8 Byte/s) 1 Thread" )
legend(loc='best')

tight_layout()
show()
exit()

