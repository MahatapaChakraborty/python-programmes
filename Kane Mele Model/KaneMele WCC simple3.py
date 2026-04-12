import numpy as np
import matplotlib.pyplot as plt


# PARAMETERS

t = 1.0
lambda_SO = 0.24
lambda_R = 0.05
lambda_v = 0.7
a = 1.0

Nk1 = 280
Nk2 = 250


# PAULI MATRICES

sigma_x = np.array([[0,1],[1,0]], dtype=complex)
sigma_y = np.array([[0,-1j],[1j,0]], dtype=complex)
sigma_z = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2)


# HAMILTONIAN

def H_k(kx, ky):
    x = kx * a / 2
    y = np.sqrt(3) * ky * a / 2

    d1 = t * (1 + 2*np.cos(x)*np.cos(y))
    d2 = lambda_v
    d3 = lambda_R * (1 - np.cos(x)*np.cos(y))
    d4 = -np.sqrt(3)*lambda_R*np.sin(x)*np.sin(y)

    d12 = -2*t*np.cos(x)*np.sin(y)
    d15 = 2*lambda_SO*(np.sin(2*x) - 2*np.sin(x)*np.cos(y))
    d23 = -lambda_R*np.cos(x)*np.sin(y)
    d24 = np.sqrt(3)*lambda_R*np.sin(x)*np.cos(y)

    Gamma1 = np.kron(I2, sigma_x)
    Gamma2 = np.kron(I2, sigma_z)
    Gamma3 = np.kron(sigma_x, sigma_y)
    Gamma4 = np.kron(sigma_y, sigma_y)
    Gamma5 = np.kron(sigma_z, sigma_y)

    def comm(A,B):
        return (A@B - B@A)/(2j)

    H = (d1*Gamma1 + d2*Gamma2 + d3*Gamma3 + d4*Gamma4 +
         d12*comm(Gamma1,Gamma2) +
         d15*comm(Gamma1,Gamma5) +
         d23*comm(Gamma2,Gamma3) +
         d24*comm(Gamma2,Gamma4))

    return H


# BZ VECTORS

b1 = (2*np.pi/a) * np.array([1, 1/np.sqrt(3)])
b2 = (2*np.pi/a) * np.array([-1, 1/np.sqrt(3)])


# WCC FUNCTION

def wcc(n):

    def k(i):
        return (i/Nk1)*b1 + (n/Nk2)*b2

    def occupied_vecs(k):
        eigvals, eigvecs = np.linalg.eigh(H_k(k[0], k[1]))
        idx = np.argsort(eigvals)
        u1 = eigvecs[:, idx[0]]
        u2 = eigvecs[:, idx[1]]
        return u1, u2

    M = np.eye(2, dtype=complex)

    for i in range(Nk1):

        u1_k, u2_k = occupied_vecs(k(i))
        u1_kp, u2_kp = occupied_vecs(k((i+1) % Nk1))

        # Explicit overlap matrix (your style)
        M11 = np.vdot(u1_k, u1_kp)
        M12 = np.vdot(u1_k, u2_kp)
        M21 = np.vdot(u2_k, u1_kp)
        M22 = np.vdot(u2_k, u2_kp)

        Mcurrent = np.array([[M11, M12],
                             [M21, M22]])

        
        #Ulink, _, Vh = np.linalg.svd(Mcurrent)
        #Mcurrent = Ulink @ Vh

        M = np.matmul(M, Mcurrent)

    W = M

    eigvals, _ = np.linalg.eig(W)
    phases = np.sort(np.angle(eigvals))

    return phases


# COMPUTE WCC FLOW

b1portions = np.linspace(0,1,Nk2)

wccupper = []
wcclower = []

for n in range(Nk2):
    phases = wcc(n)
    wcclower.append(phases[0])
    wccupper.append(phases[1])

wccupper = np.array(wccupper)
wcclower = np.array(wcclower)


# PLOT

plt.figure(figsize=(7,5))

plt.plot(b1portions, wccupper, color='red')
plt.plot(b1portions, wcclower, color='red')

plt.xlabel("b2 direction (ky index)")
plt.ylabel("WCC (phase)")
plt.ylim(-np.pi,np.pi)
plt.title(f"WCC flow λ_SO={lambda_SO}, λ_v={lambda_v}, λ_R={lambda_R}")
plt.grid(alpha=0.3)

plt.show()
