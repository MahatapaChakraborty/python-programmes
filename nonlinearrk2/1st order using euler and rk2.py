#dy/dx=(x^2)-y, y(-3.9204)=0,y(3.01)=5.015, original solution=(1-x)^2 +1 - 0.5*exp(-x)
#y(x0+h)-y(x0)=h*((x0)^2 - y0)
#the iterative formula is
#                       y(n+1)=y(n)+h*(xn^2 - yn)
import numpy as np
import matplotlib.pyplot as plt

x0=-3.9204
y0=0
xn=5.1
n=90#total no. of points-1
h=(xn-x0)/n
xx=[]
for i in range(n+1):
    x=x0+i*h
    xx.append(x)

#print(xx)

def originalsoln(zz):
    originalsoln=[]
    for i in range(len(zz)):
        originalsol1=((1-zz[i])**2)+1-(0.5*np.exp(-zz[i]))
        originalsoln.append(originalsol1)
    return(originalsoln)


def eulersoln(xx,y0):
    eulersoln=[]
    yinitial=y0
    yold=yinitial

    for i in range(len(xx)):
        ynew=yold+h*(((xx[i])**2)-yold)
        eulersoln.append(ynew)
        yold=ynew

    return(eulersoln)

#print(eulersoln(xx,y0))

def rk2soln(eulersoln,xx,y0):
    rksoln=[]
    yinitial=y0
    yrkold=yinitial

    for i in range(len(eulersoln)-1):
        yrknew=yrkold+(h/2)*(((((xx[i])**2)-yrkold)+(((xx[i+1])**2)-eulersoln[i])))
        rksoln.append(yrknew)
        yrkold=yrknew

    return(rksoln)

#print(rk2soln(eulersoln(xx,y0),xx,y0))

xxnew=[]
yyeulernew=[]
for i in range(len(xx)-1):
    xnew=xx[i]
    yeulernew=eulersoln(xx,y0)[i]
    yyeulernew.append(yeulernew)
    xxnew.append(xnew)

#print(xxnew)
#print(yyeulernew)

plt.xlabel('x')
plt.ylabel('y')
plt.plot(xxnew,originalsoln(xxnew),color='r',label='original soln')
plt.plot(xxnew,yyeulernew,color='k',label='eulers soln.')
plt.plot(xxnew,rk2soln(eulersoln(xx,y0),xx,y0),color='b',label='rk2 soln')
plt.title('dy/dx=(x^2)-y, y(-3.92)=0,h=0.100226')
plt.legend()
plt.grid()
plt.show()





    
    
