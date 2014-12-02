from numpy import *
from matplotlib.pyplot import *
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t", help="Show statistics for measurements for same length", action="store_true")
args = parser.parse_args()

stride  = 1
ardata = []

#ardata.append( genfromtxt( "./2014-11-12_20-43_data.txt", comments='#' ) )
ardata.append( genfromtxt( "./2014-11-12_21-35_data.txt", comments='#' ) )
times = genfromtxt( "./2014-11-12_21-35_Times.txt", comments='#' )
arlabels = [ "Home",
           ]


if args.t:
    dataset = 0
    def hist_key_analyzer(event):
        global dataset
        print event.key
        print "Dataset:",dataset
        if (event.key == '+' or event.key == "right"):
            dataset = (dataset+1) % len(times[:,0])
        if (event.key == '-' or event.key == "left"):
            dataset = (dataset-1) % len(times[:,0])
        figure( hist_fig.number )
        clf()
        suptitle( str(times[dataset,0]) + " Bytes" )
        hist( times[dataset,2:], bins=50 )
        draw()

    hist_fig = figure( figsize=(6,4) )
    connect('key_press_event', hist_key_analyzer)
    suptitle( str(times[dataset,0]) + " Bytes" )
    hist( times[dataset,2:], range=(986000,992000), bins=50 )
    tight_layout()
    show()
    exit()


# integration intervals	processors used	time / CPU clocks	sigma CPU clocks	result	relative Error
imsize  = 0
idsize  = 1
it      = 2
iterr   = 3
istride = 4

fig = figure(figsize=(8,6))
ax = subplot(111, xscale='log', xlabel="Matrix Size / Bytes",
                  yscale='linear', ylabel="Gflops" )

cachesizes = [ 32*1024, 256*1024, 20*1024**2 ]
for x in cachesizes:
    axvline( x,color='k',ls='dashed')

for stride in range(3):
    for k in range(len(ardata)):
        #x = ardata[k][:,imsize]**2*ardata[k][:,idsize]
        data    = ardata[k]
        mask    = data[:,istride] == 64*stride
        xlength = data[mask,imsize]
        #mask    = data[:,imsize] == 32*1+1
        #xstride = data[mask,istride]
        y = 2*data[mask,imsize] / data[mask,it]
        errorbar( xlength, y/1e9, yerr=data[mask,iterr]/data[mask,it]*y/1e9, label="Stride:"+str(64*stride) )

legend(loc='best')

tight_layout()
show()
exit()

