set fit errorvariables
ml=-1
nl=1000
gl(x)=ml*x+nl
fit gl(x) "Pa234.dat" using 1:(log($2)):(sqrt($2)/$2) via ml,nl

print "FIT_NDF = Number of degrees of freedom : ", FIT_NDF
print "FIT_WSSR = Weighted sum-of-squares residual : ", FIT_WSSR
print "FIT_STDFIT = sqrt(WSSR/NDF) : ", FIT_STDFIT
print "m_err : ",ml_err
print "n_err : ",nl_err

set terminal png size 640,480
set output "plot.png"
plot "Pa234.dat" using 1:(log($2)):(sqrt($2)/$2) title "Linearized Data" with errorbars,\
	 -0.008069*x+7.21506 title "Own Program: -0.902998*x+458.974" linewidth 4,\
	 gl(x)
	 
reset
set terminal wxt
plot "Pa234.dat" using 1:(log($2)):(sqrt($2)/$2) title "Linearized Data" with errorbars,\
	 -0.008069*x+7.21506 title "Own Program: -0.902998*x+458.974" linewidth 4,\
	 gl(x)
	 
pause -1
