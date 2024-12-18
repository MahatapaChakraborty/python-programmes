import numpy as np
import matplotlib.pyplot as plt

#x=5

#y=np.exp(5)

#print(y)

def xaxispoints(x0,xn,n):#x0=initial point xn=final point n=total no. of points-1
    xx=[]
    h=(xn-x0)/n
    x=x0

    for i in range(n+1):
        x=x0+(i*h)
        xx.append(x)

    return(xx)

#print(xaxispoints(0,1,10))#<-gives 11 points in total

def yaxispoints(xx):
    yy=[]
    for i in range(len(xx)):
        x1=np.exp(-(xx[i])**2)
        yy.append(x1)
    
    return(yy)

x0=-5
xn=5
n=10000

#print(yaxispoints(xaxispoints(x0,xn,n)))

#plt.plot(xaxispoints(x0,xn,n),yaxispoints(xaxispoints(x0,xn,n)))
#plt.xlabel("x")
#plt.ylabel("exp(-x^2)")
#plt.grid()
#plt.show()

def trapizoidal(xx,yy):
    h=(xx[1]-xx[0])
    integrationval1=0
    for i in range(1,len(xx)):
        integrationval1+=(yy[i]*h)

    integralval2=(h/2)*(yy[0]+yy[len(yy)-1])
    integrationval=integrationval1+integralval2
    return(integrationval)

trapintval=trapizoidal(xaxispoints(x0,xn,n),yaxispoints(xaxispoints(x0,xn,n)))
print("integral of exp(-x**2) from x=-5 to x=5 is",trapintval)
        
        
        


