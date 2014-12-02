from numpy import *
from matplotlib.pyplot import *

a=1
b=4
def f1a(x):
    return  3*sqrt(x)
def f1b(x):
    return -3*sqrt(x)
def f2(x):
    return x**2-4*x+6
def f(x):
    return f1a(x)-f2(x)
def fnum(x,N=3):
    dx = (b-a)/N
    res = 0*x
    for kres in range(len(res)):
        for i in range(N):
            xa = a + dx*i
            xb = a + dx*(i+1)
            if xa <= x[kres] and x[kres] < xb:
                res[kres] = f( (xa+xb)/2. )
            #print "xa:",xa," xb:",xb," x[kres]:",x[kres]," x0:",(xa+xb)/2.," f(x0):",f((xa+xb)/2)
    return res

    
fig = figure( figsize=(10,4) )
x = linspace(0,5,1001)
ax = subplot(121)
plot( x, f1a(x), 'k-', label=r"$\pm 3 \sqrt{x}$" )
plot( x, f1b(x), 'k-' )
ax.fill_between(x, f1a(x), f2(x), where=f1a(x)>=f2(x), facecolor='blue', alpha=0.5, interpolate=True)
plot( x, f2(x) , 'k--', label=r"$x^2+4x+6$" )
ylim( (1,8) )
legend( loc='lower right' )
ax = subplot(122)
plot( x, f(x), 'k-', label=r"$3 \sqrt{x} - (x^2+4x+6)$" )

plot( x, fnum(x), 'b-' )
ax.fill_between(x, fnum(x), 0, where=fnum(x)>=0, facecolor='blue', alpha=0.5, interpolate=True)
plot( [1.5,2.5,3.5], f(array([1.5,2.5,3.5])), 'bo' )

ax.text( 1.5, f(1.5), r"$f(x_0)$", verticalalignment='bottom', horizontalalignment='center' )
ax.text( 2.5, f(2.5), r"$f(x_1)$", verticalalignment='bottom', horizontalalignment='center' )
ax.text( 3.5, f(3.5), r"$f(x_2)$", verticalalignment='bottom', horizontalalignment='center' )

ax.text( 1.5, 0, r"$a+\frac{h}{2}+0\cdot h$", verticalalignment='top', horizontalalignment='center' )
ax.text( 2.5, 0, r"$a+\frac{h}{2}+1\cdot h$", verticalalignment='top', horizontalalignment='center' )
ax.text( 3.5, 0, r"$a+\frac{h}{2}+2\cdot h$", verticalalignment='top', horizontalalignment='center' )

ax.text( 0.4, 2.8, r"$a=1$"+"\n"+r"$b=4$"+"\n"+r"$N=3$"+"\n"+r"$h=\frac{b-a}{N}$", verticalalignment='top', horizontalalignment='left' )
legend( loc='lower right' )

ylim( (-1,3) )
tight_layout()
fig.savefig("CentralIntegrationExplanation.pdf")

show()
exit()
