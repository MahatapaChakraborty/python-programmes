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
N=5
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
t2=0.1347#1
M=0.7#np.sqrt(3)*3
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
        #python has different convention than math matrix here A[i,j]
        #means ith row and j-th column with makes j the x indeces and i is
        #the y indeces

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
plt.title("2X2 Haldane model energy difference t1=1 t2=0.1347 M=0.7 phi=pi/2")
plt.show()'''

def k_path(points, n_per_segment):
    klist = []
    totallength = [0.0]
    tick_positions = [0.0]

    total_len = 0.0

    for i in range(len(points) - 1):
        segment = np.linspace(points[i], points[i+1],
                               n_per_segment, endpoint=False)

        for k in segment:
            if len(klist) > 0:
                dk = np.linalg.norm(k - klist[-1])
                total_len += dk
                totallength.append(total_len)
            klist.append(k)

        tick_positions.append(total_len)

    # final point
    dk = np.linalg.norm(points[-1] - klist[-1])
    total_len += dk
    klist.append(points[-1])
    totallength.append(total_len)

    return np.array(klist), np.array(totallength), tick_positions

t1  = 1.0
t2  = 0#0.1347
M   = 0.7
phi = np.pi/2
pointspersegment=50

Gamma = np.array([0.0, 0.0])
K     = (2*b1 + b2) / 3
Kp    = (b1 + 2*b2) / 3
Mpt   = (b1 + b2) / 2
highsympoints = [Gamma, K, Mpt, Kp, Gamma]

kpath, totallength, ticks = k_path(highsympoints, pointspersegment)

Elower, Eupper = [], []
for k in kpath:
    e = energies(k[0], k[1], t1, t2, M, phi)
    Elower.append(e[0])
    Eupper.append(e[1])

Elower = np.array(Elower)
Eupper = np.array(Eupper)

plt.figure(figsize=(7,5))

plt.scatter(totallength, Elower, s=2, color='black')
plt.scatter(totallength, Eupper, s=2, color='black')
#plt.plot(totallength, Elower, lw=1, color='black')
#plt.plot(totallength, Eupper, lw=1, color='black')

plt.xticks(
    ticks,
    [r'$\Gamma$', r'$K$', r'$M$', r"$K'$", r'$\Gamma$']
)

plt.xlabel("k-path length")
plt.ylabel("Energy")
plt.title("Haldane model band structure M=0.7 t2=0 phi=pi/2")
plt.grid(alpha=0.3)

plt.show()
