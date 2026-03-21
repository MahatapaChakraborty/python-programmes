import numpy as np
import matplotlib.pyplot as plt


t  = 1.0
t2 = 0.2

Nk = 30             # k-mesh (Nk x Nk)
M_vals   = np.linspace(-3, 3, 100)
phi_vals = np.linspace(-np.pi, np.pi, 100)


b1 = np.array([2*np.pi/np.sqrt(3),  2*np.pi/3])
b2 = np.array([-2*np.pi/np.sqrt(3), 2*np.pi/3])


delta = np.array([
    [0, 1],
    [ np.sqrt(3)/2, -1/2],
    [-np.sqrt(3)/2, -1/2]
])

bvec = np.array([
    delta[1] - delta[2],
    delta[2] - delta[0],
    delta[0] - delta[1]
])


def Haldane_H(kx, ky, M, phi):
    k = np.array([kx, ky])

    dx = -t * np.sum(np.cos(delta @ k))
    dy = -t * np.sum(np.sin(delta @ k))
    dz = M - 2*t2*np.sin(phi) * np.sum(np.sin(bvec @ k))

    H = np.array([[ dz, dx - 1j*dy],
                  [ dx + 1j*dy, -dz]])
    return H


def chern_number(M, phi):
    # k-grid
    k1 = np.linspace(0, 1, Nk, endpoint=False)
    k2 = np.linspace(0, 1, Nk, endpoint=False)

    # eigenvectors
    u = np.zeros((Nk, Nk, 2), dtype=complex)

    for i, a in enumerate(k1):
        for j, b in enumerate(k2):
            k = a*b1 + b*b2
            _, vecs = np.linalg.eigh(Haldane_H(k[0], k[1], M, phi))
            u[i, j] = vecs[:, 0]   # lower band

    F = 0.0
    for i in range(Nk):
        for j in range(Nk):
            ip = (i+1) % Nk
            jp = (j+1) % Nk

            U1 = np.vdot(u[i, j],   u[ip, j])
            U2 = np.vdot(u[ip, j],  u[ip, jp])
            U3 = np.vdot(u[ip, jp], u[i, jp])
            U4 = np.vdot(u[i, jp],  u[i, j])

            F += np.log(U1*U2*U3*U4).imag

    return F / (2*np.pi)


chern_map = np.zeros((len(M_vals), len(phi_vals)))

for i, M in enumerate(M_vals):
    for j, phi in enumerate(phi_vals):
        chern_map[i, j] = np.round(chern_number(M, phi))


plt.figure(figsize=(7,5))
plt.imshow(chern_map,
           extent=[phi_vals[0], phi_vals[-1], M_vals[0], M_vals[-1]],
           origin='lower',
           aspect='auto',
           cmap='RdBu')

plt.colorbar(label="Chern number")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$M$")
plt.grid(alpha=0.5)
plt.title("Haldane Model Chern Number (Plaquette Method) 30X30 K-grid 100 points for M and phi t1=1 t2=0.2")
plt.show()
