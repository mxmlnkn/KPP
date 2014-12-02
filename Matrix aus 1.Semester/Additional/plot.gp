m=-1
n=1000
g(x)=m*x+n
set fit errorvariables
fit g(x) "Pa234.dat" using 1:2 via m,n

print "FIT_NDF = Number of degrees of freedom : ", FIT_NDF
print "FIT_WSSR = Weighted sum-of-squares residual : ", FIT_WSSR
print "FIT_STDFIT = sqrt(WSSR/NDF) : ", FIT_STDFIT
print "m_err : ",m_err
print "n_err : ",n_err
	  
mw=-1
nw=1000
gw(x)=mw*x+nw
fit gw(x) "Pa234.dat" using 1:2:(sqrt($2)) via mw,nw

print "FIT_NDF = Number of degrees of freedom : ", FIT_NDF
print "FIT_WSSR = Weighted sum-of-squares residual : ", FIT_WSSR
print "FIT_STDFIT = sqrt(WSSR/NDF) : ", FIT_STDFIT
print "m_err : ",mw_err
print "n_err : ",nw_err


ml=-1
nl=1000
gl(x)=ml*x+nl
fit gl(x) "Pa234.dat" using 1:(log($2)):(sqrt($2)/$2) via ml,nl

print "FIT_NDF = Number of degrees of freedom : ", FIT_NDF
print "FIT_WSSR = Weighted sum-of-squares residual : ", FIT_WSSR
print "FIT_STDFIT = sqrt(WSSR/NDF) : ", FIT_STDFIT
print "m_err : ",ml_err
print "n_err : ",nl_err



N0=1000
T12=1000
#f(x)=log(N0)-x*log(2)/T12
f(x)=N0*exp( -x*log(2)/T12 )
fit f(x) "Pa234.dat" using 1:2 via N0,T12

print "FIT_NDF = Number of degrees of freedom : ", FIT_NDF
print "FIT_WSSR = Weighted sum-of-squares residual : ", FIT_WSSR
print "FIT_STDFIT = sqrt(WSSR/NDF) : ", FIT_STDFIT
print "m_err : ",N0_err
print "n_err : ",T12_err





set multiplot layout 1,2

plot "Pa234.dat" using 1:2:(sqrt($2)) with errorbars,\
	 -0.902998*x+458.974 title "Own Program weighted: -0.902998*x+458.974" linewidth 4,\
	 -1.9299*x+833.264 title "Own Program unweighted: -1.9299*x+833.264" linewidth 4,\
	 gw(x) title "Gnuplot weighted",\
	 g(x) title "Gnuplot unweighted",\
	 1359.75*2**(-x/85.9025),\
	 f(x)

plot "Pa234.dat" using 1:(log($2)):(sqrt($2)/$2) with errorbars,\
	 "Pa234.dat" using 1:(log($2)),\
	 -0.008069*x+7.21506 title "Own Program unweighted: -0.008069*x+7.21506" linewidth 4,\
     gl(x)

unset multiplot

pause -1
