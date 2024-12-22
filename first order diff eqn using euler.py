#solving differential equaton using euler's method
#dy/dx =cos(x), y(0)=0

import numpy as np
import matplotlib.pyplot as plt

xx=[]
x0=0
xn=2*(np.pi)
n=95
y0=0

h=(xn-x0)/n
for i in range(n+1):
    x1=x0+i*h
    xx.append(x1)

#print(xx)

cosx=[]
for i in range(len(xx)):
    cosx1=np.cos(xx[i])
    cosx.append(cosx1)

#print(expx)

yydiffres=[]
yydiffres1=y0

for i in range(len(xx)):
    yydiffres1+=h*cosx[i]
    yydiffres.append(yydiffres1)

#print(yydiffres)

plt.xlabel('x')
plt.ylabel('function')
plt.plot(xx,np.sin(xx),color='r',label='original soln.')
plt.plot(xx,yydiffres,color='b',label='eulers soln')
plt.title('soln. of dy/dx=cosx')
plt.legend()
plt.grid()
plt.show()
