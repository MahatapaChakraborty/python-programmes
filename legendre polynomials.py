#for numerical purpose, the recurrence relation among legendre's polynomials is
            #p(n+1)=2*x*p(n)-p(n-1)-(1/(n+1))*(x*p(n)-p(n-1)) p(x,n)=n-th legendre polynomial

#write a python program to generate p(m,x) for any arbitrary m>=2 on N number of discrete point x(i)
#for i=0,1,2,...N-1 with x0=-1 and x(N-1)=1.initial condition to use the recurrence relation is given
#by specifying p(0,x)=1 and p(1,x)=x

import numpy as np
import matplotlib.pyplot as plt

xx=[]
x0=-1
xn=1
n=500#n=total no. of points-1
h=(xn-x0)/n
for i in range(n+1):
    x1=x0+i*h
    xx.append(x1)

#print(xx)

def legendrepolynomial(n,x):
    if n==0:
        return(1)
    else:
        if n==1:
            return(x)
        else:
            return(2*x*legendrepolynomial(n-1,x)-legendrepolynomial(n-2,x)-(1/n)*(x*legendrepolynomial(n-1,x)-legendrepolynomial(n-2,x)))


def legendrepolynomials(n,xx):
    legendrepolynomials=[]
    for i in range(len(xx)):
        legendrepolynomials1=legendrepolynomial(n,xx[i])
        legendrepolynomials.append(legendrepolynomials1)

    return(legendrepolynomials)
plt.xlabel('x')
plt.ylabel('function value')
plt.plot(xx,legendrepolynomials(0,xx),color='r',label='p0')
plt.plot(xx,legendrepolynomials(1,xx),color='b',label='p1')
plt.plot(xx,legendrepolynomials(2,xx),color='g',label='p2')
plt.plot(xx,legendrepolynomials(3,xx),color='y',label='p3')
plt.plot(xx,legendrepolynomials(4,xx),color='c',label='p4')
plt.plot(xx,legendrepolynomials(5,xx),color='m',label='p5')
plt.plot(xx,legendrepolynomials(6,xx),color='k',label='p6')
plt.title('legendre polynomials')
plt.legend()
plt.grid()
plt.show()
    
































#wri
            
