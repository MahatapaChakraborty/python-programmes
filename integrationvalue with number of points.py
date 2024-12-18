import numpy as np
import matplotlib.pyplot as plt

def integrationvalue(x0,xn,n):
    h=(xn-x0)/n
    xx=[]
    
    for i in range(n+1):
        x1=x0+i*h
        xx.append(x1)

    #return(xx)
#print(integrationvalue(0,1,10))

    yy=[]
    for i in range(len(xx)):
        y=np.exp(-(xx[i])**2)
        yy.append(y)


    integrationval1=0
    for i in range(1,len(xx)):
        integrationval1+=yy[i]*h

    integrationval2=(h/2)*(yy[0]+yy[len(yy)-1])

    integrationval=integrationval1+integrationval2

    return(integrationval)

#print(integrationvalue(-5,5,100000)

xx1=[]
for i in range(1,60):
    x1=i
    xx1.append(x1)

#print(xx1)

yy1=[]
for i in range(1,60):
    y1=integrationvalue(-5,5,i)
    yy1.append(y1)

plt.plot(xx1,yy1)
plt.xlabel("no. of points/rectangles taken")
plt.ylabel("value of the integration")
plt.grid()
plt.show()
