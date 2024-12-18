import numpy as np
import matplotlib.pyplot as plt
x0=0
xn=1
n=30#n=total no. of points-1
h=((xn-x0)/n)
xx=[]
for i in range(n+1):
    x1=x0+i*h
    xx.append(x1)

#print(xx)

yy=[]
for i in range(len(xx)):
    y1=np.exp(-xx[i])
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
#plt.xlabel('angle')
#plt.ylabel('function')
#plt.plot(xxf,yyforward,color='g',label='dsinx/dx forward')
#plt.plot(xxf,yyoriginal,color='r',label='sin')
#plt.legend()
#plt.show()

#backward difference

yybackward=[]
for i in range(1,len(yy)-1):
    yybackward1=(yy[i]-yy[i-1])/h
    yybackward.append(yybackward1)

#print(yybackward)
#plt.xlabel('angle')
#plt.ylabel('function')
#plt.plot(xxf,yyforward,color='g',label='dsinx/dxforward')
#plt.plot(xxf,yybackward,color='b',label='dsinx/dxbackward')
#plt.plot(xxf,yyoriginal,color='r',label='sinx')
#plt.grid()
#plt.legend()
#plt.show()

#central difference
yycentral=[]
for i in range(1,len(yy)-1):
    yycentral1=(yy[i+1]-yy[i-1])/(2*h)
    yycentral.append(yycentral1)

#print(yycentral)
#print(yybackward)
plt.xlabel('x')
plt.ylabel('function')
plt.plot(xxf,yycentral,color='y',label='dexp-x/dxcentral')
plt.plot(xxf,yyforward,color='g',label='dexp-x/dxforward')
plt.plot(xxf,yybackward,color='b',label='dexp-x/dxbackward')
plt.plot(xxf,yyoriginal,color='r',label='exp-x')
plt.grid()
plt.legend()
plt.show()


    

    
