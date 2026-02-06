import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi,voronoi_plot_2d

#lattice vectors
a1=np.array([1,0])
a2=np.array([1/2,np.sqrt(3)/2])
#reciprocal lattice vectors
b1=np.array([2*np.pi,-(2*np.pi)/np.sqrt(3)])
b2=np.array([0,4*np.pi/np.sqrt(3)])
#nearest neighbors
n1=np.array([0.0,1/(np.sqrt(3))])
n2=np.array([0.5,-1/(2*np.sqrt(3))])
n3=np.array([-0.5,-1/(2*np.sqrt(3))])
#next nearest neighbors
nn1=n2-n3
nn2=n3-n1
nn3=n1-n2
#print(n1)
N=3
M1=np.linspace(1,N,N)
voronoiseeds=[]
for i in range(len(M1)):
    for j in range(len(M1)):
        seed=i*b1 + j*b2
        voronoiseeds.append(seed)

#print(voronoiseeds)
vor=Voronoi(voronoiseeds)
#fig=voronoi_plot_2d(vor)
#plt.show()
#print(vor.vertices)
Nk=200
kxrange=np.linspace(np.min(vor.vertices[:,0]),np.max(vor.vertices[:,0]),Nk)
kyrange=np.linspace(np.min(vor.vertices[:,1]),np.max(vor.vertices[:,1]),Nk)
kx,ky=np.meshgrid(kxrange,kyrange)

t1=1
t2=1
M=np.sqrt(3)*3
phi=np.pi/2

def energies(kx,ky,t1,t2,M,phi):
    unitmat=np.array([[1,0],[0,1]])
    sigmaz=np.array([[1,0],[0,-1]])
    sigmax=np.array([[0,1],[1,0]])
    sigmay=np.array([[0,-1j],[1j,0]])
    kdotn1=kx*n1[0] + ky*n1[1]
    kdotn2=kx*n2[0] + ky*n2[1]
    kdotn3=kx*n3[0] + ky*n3[1]
    kdotnn1=kx*nn1[0] + ky*nn1[1]
    kdotnn2=kx*nn2[0] + ky*nn2[1]
    kdotnn3=kx*nn3[0] + ky*nn3[1]
    coskdotnextnearestneigh=np.cos(kdotnn1)+np.cos(kdotnn2)+np.cos(kdotnn3)
    sinkdotnextnearestneigh=np.sin(kdotnn1)+np.sin(kdotnn2)+np.sin(kdotnn3)

    coskdotnearestneigh=np.cos(kdotn1)+np.cos(kdotn2)+np.cos(kdotn3)
    sinkdotnearestneigh=np.sin(kdotn1)+np.sin(kdotn2)+np.sin(kdotn3)

    H0=2*t2*np.cos(phi)*coskdotnextnearestneigh*unitmat
    H1=t1*(coskdotnearestneigh*sigmax + sinkdotnearestneigh*sigmay)
    H2=(M-2*t2*np.sin(phi)*sinkdotnextnearestneigh)*sigmaz

    H=H0+H1+H2
    E=np.linalg.eigvals(H)
    Ereal=[]
    for i in range(len(E)):
        Ereal.append(E[i].real)
    Ereal.sort()
    return(Ereal)

#print(energies(0,0,0,0,0,0))

lowerenergies=np.zeros((Nk,Nk),dtype=complex)
for i in range(Nk):
    for j in range(Nk):
        lowerenergies[j,i]=energies(kxrange[i],kyrange[j],t1,t2,M,phi)[0]
        # in python convention A[i,j] i is row index or in math convention the y co-ordinate
#in python convention A[i,j] j is column index or in math convention the x co-ordinate
#pcolormesh uses kx as x axis and ky as y axis 

upperenergies=np.zeros((Nk,Nk),dtype=complex)
for i in range(Nk):
    for j in range(Nk):
        upperenergies[j,i]=energies(kxrange[i],kyrange[j],t1,t2,M,phi)[1]

'''fig=plt.figure(figsize=(10,7))
ax=fig.add_subplot(111,projection='3d')

ax.plot_surface(kx,ky,lowerenergies,color='red')
ax.plot_surface(kx,ky,upperenergies,color='blue')

ax.set_title('2X2 Haldane model t1=t2=1 M=3root3 phi=pi/2')
#ax.set_zlim(-0.1,0.1)
ax.set_xlabel('Kx')
ax.set_ylabel('Ky')
ax.set_zlabel('Ek')


plt.show()'''

'''plt.figure(figsize=(10,7))
plt.pcolormesh(kx,ky,(upperenergies-lowerenergies).real,shading='auto',cmap='coolwarm')#diff between 1st and second band
plt.colorbar(label="Difference between two bands")
plt.xlabel("k_x values")
plt.ylabel("k_y values")
plt.title("2X2 Haldane model energy difference t1=t2=1 M=3*(root3) phi=pi/2")
plt.show()'''

'''Gamma = np.array([0.0, 0.0])
K = (2*b1 + b2) / 3
Kp=(b1+2*b2)/3
Mpt = (b1+b2) / 2#high symmetry points


def k_path(points, n_per_segment=20):
    klist = []
    tick_positions = [0]
    total_len = 0

    for i in range(len(points)-1): #points is a list of 5 high symmetry points
        k_start = points[i]
        k_end = points[i+1]

        segment = np.linspace(k_start, k_end, n_per_segment, endpoint=False)#gives list of points
        # starting from k_start. it's an array of 2D points like [[_,_],[_,_],...]
        klist.append(segment)

        total_len += len(segment)#counts total no. of 2D points
        tick_positions.append(total_len)#tick position has 5 entries each
        #gives total no. of 2D points before a certain symmetry point

    klist.append(np.array([points[-1]]))# appends the gamma point at end of the
    #list. we started from gamma point and end also at gamma point. while returning to
    #gamma point from Mpoint we used endpoint=false so we have to by hand put gamma at the end
    kpath = np.vstack(klist)#klist=[segment1=[[_,_],[_,_],...],segment2=[[_,_],[_,_],...],..]
    #vstack makes kpath=[[_,_],[_,_],..[_,_]]

    return kpath, tick_positions



t1 = 1.0
t2 = 0.1347
M = 0.7
phi = np.pi/2    


# Generating energies along high symmetry path


points = [Gamma, K, Mpt, Kp, Gamma]
kpath, tick_positions = k_path(points, n_per_segment=200)
#print(kpath)
E1 = []#E1 and E2 are lists
E2 = []

for i in range(len(kpath)):
    kx=kpath[i][0]
    ky=kpath[i][1]
    
    e = energies(kx, ky, t1, t2, M, phi)
    E1.append(e[0])
    E2.append(e[1])

E1 = np.array(E1)
E2 = np.array(E2)


plt.figure(figsize=(7,5))

plt.plot(E1, 'k', lw=1.5)
plt.plot(E2, 'k', lw=1.5)

plt.xticks(tick_positions, [r'$\Gamma$', r'$K$', r'$M$',r"$K'$" ,r'$\Gamma$'])
plt.ylabel("Energy")
plt.title("Haldane model: high-symmetry band structure t1=1,t2=0.1347 M=0.7 phi=pi/2")
plt.grid(alpha=0.3)

plt.show()'''




    
    
    
    
        

