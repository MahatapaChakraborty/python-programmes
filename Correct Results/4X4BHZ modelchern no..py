import numpy as np
import matplotlib.pyplot as plt
#working formula
#berry phase around one plaquette=Imlndet(M12M23M34M41)
#Mij_(mn)=<um(ki)|un(kj)> m,n->Band index (occupied bands)
def bigBHZCno(N,m):
    Kx=np.linspace(0,2*np.pi,N,endpoint=False)
    Ky=np.linspace(0,2*np.pi,N,endpoint=False)

    def eigenvectors(kx,ky):
        H12=np.sin(kx)- (1j)* (np.sin(ky))
        H21=np.sin(kx)+ (1j)* (np.sin(ky))
        H11=m + np.cos(kx) + np.cos(ky)
        H22=-H11
        H13=0
        H14=0
        H23=0
        H24=0
        H31=0
        H32=0
        H41=0
        H42=0
        H34=np.sin(kx)- (1j)* (np.sin(ky))
        H43=np.sin(kx)+ (1j)* (np.sin(ky))
        H33=m + np.cos(kx) + np.cos(ky)
        H44=-H33

        H=np.array([[H11,H12,H13,H14],
                    [H21,H22,H23,H24],
                    [H31,H32,H33,H34],
                    [H41,H42,H43,H44]])

        E,v=np.linalg.eig(H)
        idx=np.argsort(E)#sorts E and gives their positions from smallest to largest
        V=v[:,idx]
        V0=V[:,0]
        V2=V[:,1]#two lower bands taken
        return(V0/np.linalg.norm(V0) , V2/np.linalg.norm(V2))
    #return(eigenvectors(3.79,4)[0])
#print(bigBHZCno(3,2))

    #def Mmatrix(vector1,vector2):
        #Mmatrix=np.array([[np.vdot(vector1,vector1),np.vdot(vector1,vector2)],[np.vdot(vector2,vector1),np.vdot(vector2,vector2)]])
        #return(Mmatrix)

    sumsq=0.0

    for i in range(N):
        for j in range(N):
            A0 = eigenvectors(Kx[i], Ky[j])[0]#A,B,C,D diff. points around the plaquette
            A1 = eigenvectors(Kx[i], Ky[j])[1]#0,1 different band state vectors at the same (kx,ky) point
            B0 = eigenvectors(Kx[i], Ky[(j+1)%N])[0]
            B1 = eigenvectors(Kx[i], Ky[(j+1)%N])[1]
            C0 = eigenvectors(Kx[(i+1)%N], Ky[(j+1)%N])[0]
            C1 = eigenvectors(Kx[(i+1)%N], Ky[(j+1)%N])[1]
            D0 = eigenvectors(Kx[(i+1)%N], Ky[j])[0]
            D1 = eigenvectors(Kx[(i+1)%N], Ky[j])[1]

            M12=np.array([[np.vdot(A0,B0),np.vdot(A0,B1)],[np.vdot(A1,B0),np.vdot(A1,B1)]])
            M23=np.array([[np.vdot(B0,C0),np.vdot(B0,C1)],[np.vdot(B1,C0),np.vdot(B1,C1)]])
            M34=np.array([[np.vdot(C0,D0),np.vdot(C0,D1)],[np.vdot(C1,D0),np.vdot(C1,D1)]])
            M41=np.array([[np.vdot(D0,A0),np.vdot(D0,A1)],[np.vdot(D1,A0),np.vdot(D1,A1)]])

            M12M23=np.matmul(M12,M23)
            M34M41=np.matmul(M34,M41)
            M12M23M34M41=np.matmul(M12M23,M34M41)
            detoneplaquette=np.linalg.det(M12M23M34M41)
            lndet=np.log(detoneplaquette/np.abs(detoneplaquette))
            phase=lndet.imag
            sumsq+=-phase
    return(np.round(sumsq/(2*np.pi)))

#print(bigBHZCno(50,-1))
#result 1.9999999999999971

M=np.linspace(-5,5,150)
C=[]
for i in range(len(M)):
    c=bigBHZCno(50,M[i])
    C.append(c)


plt.xlabel('m (150 points between +-5)')
plt.ylabel('Chern Number')
plt.title('4X4 BHZ Model Variation of Chern No. (50X50 grid taken)(two lower bands)')
plt.axhline(y=2,color='b',linestyle='--')
plt.axhline(y=-2,color='b',linestyle='--')
plt.plot(M,C)

plt.show()

            
    
            
            
            
            
            
            

    


           


        
            
