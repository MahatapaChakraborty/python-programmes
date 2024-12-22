import numpy as np
import matplotlib.pyplot as plt
x0=0
xn=(np.pi)*2
n=200#n=total no. of points-1
h=((xn-x0)/n)
xx=[]
for i in range(n+1):
    x1=x0+i*h
    xx.append(x1)

#print(xx)

yy=[]
for i in range(len(xx)):
    y1=np.sin(xx[i])
    yy.append(y1)

yyoriginal=[]
for i in range(len(yy)-2):
    yyoriginal1=yy[i]
    yyoriginal.append(yyoriginal1)

    

#plt.plot(xx,yy,color='g',label='sin')
#plt.legend()
#plt.show()

#forward difference

yyforward=[]

for i in range(1,len(xx)-1):
    yf1=(yy[i+1]-yy[i])/h
    yyforward.append(yf1)

#print(yyforward)

xxf=[]
for i in range(len(xx)-2):
    xxf1=xx[i]
    xxf.append(xxf1)

#print(xxf)
plt.xlabel('angle')
plt.ylabel('function')
plt.plot(xxf,yyforward,color='g',label='dsinx/dx')
plt.plot(xxf,yyoriginal,color='r',label='sin')
plt.legend()
plt.show()

    
