from numpy import *
from matplotlib.pyplot import *
import argparse
from scipy.optimize import curve_fit
import os.path

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--file", help="Directory which contains 'simOutput'")
parser.add_argument("-ib", "--fileb", help="Trapeze")
parser.add_argument("-ic", "--filec", help="Center Float")
parser.add_argument("-isp1", help="")
parser.add_argument("-isp2", help="")
parser.add_argument("-isp3", help="")
parser.add_argument("--stride", help="Stride")
parser.add_argument("--task"    , action='store_true')
parser.add_argument("--relerr"  , action='store_true')
parser.add_argument("--overhead", action='store_true')
parser.add_argument("--speedup" , action='store_true')
parser.add_argument("--overview", action='store_true')
args = parser.parse_args()

if not args.task and not args.relerr and not args.overhead and not args.speedup and not args.orverview:
    args.task     = True
    args.relerr   = True
    args.overhead = True
    args.speedup  = True
    args.overview = True

if args.file == None:
    args.file  = 'TaurusCenter/2014-11-2_20-19'
    args.fileb = 'TaurusTrapeze/2014-11-2_19-52'
    args.filec = 'TaurusCenterFloat/2014-11-2_20-46'
else:
    if args.fileb == None:
        args.fileb = ''
    if args.filec == None:
        args.filec = ''

if args.stride:
    stride = args.stride
else:
    stride = 1
    
fname = "./"+args.file+"_data.txt"
if os.path.isfile(fname):
    data = genfromtxt( fname, comments='#' )
else:
    data = None
fname = "./"+args.fileb+"_data.txt"
if os.path.isfile(fname):
    datab = genfromtxt( fname, comments='#' )
else:
    datab = None
fname = "./"+args.filec+"_data.txt"
if os.path.isfile(fname):
    datac = genfromtxt( fname, comments='#' )
else:
    datac = None


# integration intervals	processors used	time / CPU clocks	sigma CPU clocks	result	relative Error
iint    = 0
inp     = 1
it      = 2
iterr   = 3
ires    = 4
irelerr = 5

max_nint = int( max( data[:,iint] ) )
max_np   = int( max( data[:,inp ] ) )

Sp1,Sp2,Sp3    = [],[],[]   # speedup 
Ep1,Ep2,Ep3    = [],[],[]   # efficiency
if args.isp1 == None:
    iSp1,iSp2,iSp3 = 999,1250,1300
else:
    iSp1,iSp2,iSp3 = args.isp1, args.isp2, args.isp3
iEp1,iEp2,iEp3 = iSp1,iSp2,iSp3





####################### Speedup and Parallel Efficiency ########################

if args.speedup:
    fig = figure( figsize=(10,4) )
    ax = subplot(121, xscale='log')
    xlabel("Integration intervals")
    ylabel("Speedup")
    ylim([0,1.5*max_np])
    mask1 = data[:,inp] == 1
    for NProcessors in range( 1, max_np+1 ): # last point is not in range itself!!!
        mask  = data[:,inp] == NProcessors
        xdata = data[mask, iint]
        ydata = data[mask1,it][:mask.sum()] / data[mask,it]
        plot( xdata[::stride], ydata[::stride], '.', label=str(NProcessors)+" Threads" )
        
        if NProcessors > 1:
            if max(iSp1,iSp2,iSp3) < len(xdata):
                Sp1.append( ydata[iSp1] )
                Sp2.append( ydata[iSp2] )
                Sp3.append( ydata[iSp3] )
            
    legend(loc='upper left', prop={'size':8})
    ax = subplot(122, xscale='log')
    xlabel("Integration intervals")
    ylabel("Parallel Efficency")
    ylim([0,1.5])
    mask1 = data[:,inp] == 1
    def fitfunc( N, r ):
        return 1.*N/(N+r)
    nfit = []
    rfit = []
    rerr = []
    for NProcessors in range( 1, max_np+1 ):
        mask  = data[:,inp] == NProcessors
        xdata = data[mask, iint]
        ydata = data[mask1,it][:mask.sum()] / data[mask,it] / NProcessors
        plot( xdata[::stride], ydata[::stride], '.', label=str(NProcessors)+" Threads" )
        
        if NProcessors > 1:
            # Save efficiency and speedup for 3 intervall numbers 
            if max(iEp1,iEp2,iEp3) < len(xdata):
                Ep1.append( ydata[iEp1] )
                Ep2.append( ydata[iEp2] )
                Ep3.append( ydata[iEp3] )
        
            # Fit N/(N+r) to parallel efficiency 
            popt, pcov = curve_fit(fitfunc, xdata[147:], ydata[147:], p0=[0.5] )
            perr  = np.sqrt(np.diag(pcov))
            nfit.append(NProcessors)
            rfit.append(popt[0])
            rerr.append(perr[0])
            print "Fit for ",NProcessors," Threads"
            print "    r  =",popt[0]," +- ",perr[0]
            #plot( xdata[147:], fitfunc( xdata[147:], popt[0]), 'r-' ) # for debug purposes to see if fit converged
        
    tight_layout()
    fig.savefig( "CentralSpeedupEfficiency.pdf" )
    fig.savefig( "CentralSpeedupEfficiency.png" )


# Eigentliche Aufgabenstellung -.-"
if args.task:
    if max(iSp1,iSp2,iSp3) < len(xdata):
        fig = figure( figsize=(10,4) )
        ax = subplot(121)
        xlabel("Number of Threads")
        ylabel("Speedup")
        xlim( (2,max_np) )
        plot( nfit, Sp1, 'o-', label=str(xdata[iSp1])+" Intervals" )
        plot( nfit, Sp2, 'o-', label=str(xdata[iSp2])+" Intervals" )
        plot( nfit, Sp3, 'o-', label=str(xdata[iSp3])+" Intervals" )
        legend( loc='upper left')
        ax = subplot(122)
        xlabel("Number of Threads")
        ylabel("Parallel Efficency")
        xlim( (2,max_np) )
        plot( nfit, Ep1, 'o-', label=str(xdata[iEp1])+" Intervals" )
        plot( nfit, Ep2, 'o-', label=str(xdata[iEp2])+" Intervals" )
        plot( nfit, Ep3, 'o-', label=str(xdata[iEp3])+" Intervals" )
        tight_layout()
        fig.savefig( "CentralSpeedupEfficiencyOverThreads.pdf" )
        fig.savefig( "CentralSpeedupEfficiencyOverThreads.png" )


# Plot overhead, der bei ersten plot for parallel efficiency durch fitten bestimmt wurde
if args.overhead:
    fig = figure( figsize=(6,4) )
    ax = subplot(111)
    xlabel("Number of Threads")
    ylabel(r"$\frac{C_\mathrm{Overhead}}{C_\mathrm{Integration}}$")
    xlim( (1,max_np+1) )
    plot( nfit, rfit, 'ro' )
    tight_layout()
    fig.savefig( "Overheadverhalten.pdf" )
    fig.savefig( "Overheadverhalten.png" )


############################### Relative Error #################################

if args.relerr:
    """
    fig = figure( figsize=(16,12) )
    nplotsx = int(sqrt(max_np+1))
    nplotsy = int(sqrt(max_np+1))
    if nplotsx * nplotsy < max_np:
        nplotsx += 1
    if nplotsx * nplotsy < max_np:
        nplotsy += 1
    for NProcessors in range( 1, max_np+1 ):
        ax  = subplot( nplotsx, nplotsy, NProcessors, xscale='log', yscale='log')
        title( str(NProcessors)+" Threads" )
        xlabel("Integration intervals")
        ylabel("Relative Error")
        mask  = data[:,inp] == NProcessors
        plot( data[mask, iint][::stride], data[mask,irelerr][:mask.sum():stride], '.' )
    tight_layout()
    fig.savefig( "RelativeErrorCenter.pdf")
    fig.savefig( "RelativeErrorCenter.png")

    if datab != None:
        fig = figure( figsize=(16,12) )
        for NProcessors in range( 1, 16+1 ):
            ax  = subplot( 4,4,NProcessors, xscale='log', yscale='log')
            title( str(NProcessors)+" Threads" )
            xlabel("Integration intervals")
            ylabel("Relative Error")
            mask  = datab[:,inp] == NProcessors
            plot( datab[mask, iint][::stride], datab[mask,irelerr][:mask.sum():stride], '.' )
        tight_layout()
        fig.savefig( "RelativeErrorTrapeze.pdf")
        fig.savefig( "RelativeErrorTrapeze.png")

    if datac != None:
        fig = figure( figsize=(16,12) )
        for NProcessors in range( 1, 16+1 ):
            ax  = subplot( 4,4,NProcessors, xscale='log', yscale='log')
            title( str(NProcessors)+" Threads" )
            xlabel("Integration intervals")
            ylabel("Relative Error")
            mask  = datac[:,inp] == NProcessors
            plot( datac[mask, iint][::stride], datac[mask,irelerr][:mask.sum():stride], '.' )
        tight_layout()
        fig.savefig( "RelativeErrorCenterFloat.pdf")
        fig.savefig( "RelativeErrorCenterFloat.png")
    """

    stride = 10
    fig = figure( figsize=(6,4) )
    ax = subplot(111, xscale='log', yscale='log')
    xlabel("Integration intervals")
    ylabel("Relative Error")
    mask1 = data[:,inp] == 1
    for NProcessors in range( 1, max_np+1 ):
        mask = data[:,inp] == NProcessors
        plot( data[mask, iint][::stride], data[mask,irelerr][:mask.sum()][::stride], '.' )
    #plot( data[mask1, iint][::stride], data[mask1,irelerr][:mask1.sum():stride], '.' )
        
    #Scaling Behavior
    N = 10**linspace(0.5,2.5,100) 
    plot( N,N**(-2)/1e3 ,'k--',label=r"$\propto h^{-2}$")
    N = 10**linspace(4,6.5,100) 
    plot( N,N/1e12      , 'k--',label=r"$\propto h$")

    # Beschriftung der Geraden
    text(10 , 1e-8, r"$\propto N^{-2}$", fontsize=18, color="k")
    text(2e5, 1e-8, r"$\propto N     $", fontsize=18, color="k")

    tight_layout()
    fig.savefig( "RelativeError.pdf" )
    fig.savefig( "RelativeError.png" )


##################### Some old overview. Not being saved #######################

if args.overview:
    fig = figure(figsize=(14,10))

    ax = subplot(231)
    xlabel("Integration Intervals")
    ylabel("Time / s")
    ax.xaxis.get_major_formatter().set_powerlimits((0, 2))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
    for NProcessors in range( 1, max_np+1 ):
        mask = data[:,inp] == NProcessors
        errorbar( data[mask, iint][::stride], data[mask,it][::stride], yerr=data[mask,iterr][::stride], fmt='.', label=str(NProcessors)+" Thread(s)" )
    legend()


    subplot(232)
    xlabel("Threads")
    ylabel("Time / s")
    ax.xaxis.get_major_formatter().set_powerlimits((0, 2))
    ax.yaxis.get_major_formatter().set_powerlimits((0, 2))
    xlim( [0.5,0.5+max_np] )
    mask = data[:,iint] == max_nint
    print "nint:", max_nint, " processors:",data[mask,inp]
    errorbar( data[mask,inp], data[mask,it], yerr=data[mask,iterr], fmt='b.', label="999 Integration Intervals" )


    ax = subplot(233, xscale='log')
    xlabel("Integration intervals")
    ylabel("Speedup")
    ylim([0,1.5*max_np])
    mask1 = data[:,inp] == 1
    for NProcessors in range( 2, max_np+1 ):
        mask = data[:,inp] == NProcessors
        plot( data[mask, iint][::stride], data[mask1,it][:mask.sum():stride] / data[mask,it][::stride], '.', label=str(NProcessors)+" Thread(s)" )
    legend()



    ax = subplot(234, xscale='log')
    xlabel("Integration intervals")
    ylabel("Parallel Efficency")
    ylim([0,1.5])
    mask1 = data[:,inp] == 1
    for NProcessors in range( 2, max_np+1 ):
        mask = data[:,inp] == NProcessors
        plot( data[mask, iint][::stride], data[mask1,it][:mask.sum():stride] / data[mask,it][::stride] / NProcessors, '.', label=str(NProcessors)+" Thread(s)" )
    legend()


    ax = subplot(235, xscale='log', yscale='log')
    xlabel("Integration intervals")
    ylabel("Relative Error")
    mask1 = data[:,inp] == 1
    for NProcessors in range( 2, max_np+1 ):
        mask = data[:,inp] == NProcessors
        plot( data[mask, iint][::stride], data[mask1,irelerr][:mask.sum():stride], '.', label=str(NProcessors)+" Thread(s)" )
    legend()


tight_layout()
show()
