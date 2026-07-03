# using koopman operator to detetct phase transition
import random
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.signal import savgol_filter


def energy_2d(s, J):
    l = len(s)
    E_ij = 0
    for i in range(l):
        for j in range(l):
            E_ij -= J * s[i, j] * (s[(i + 1) % l, j] + s[i, (j + 1) % l] + s[(i - 1) % l, j] + s[i, (j - 1) % l])
    return E_ij / 2

def random_lattice(num_sites):
    s = np.ones([num_sites, num_sites])
    # randomising spins so that net magnetic moment is roughly 0
    for i in range(num_sites):
        for j in range(num_sites):
            if random.random() > 0.5:
                s[i, j] = 1
            else:
                s[i, j] = -1
    return s

def correlations(s):
    n=len(s)
    sum_1=0
    sum_2=0
    for i in range(n):
        for j in range(n):
            #C1
            sum_1+= s[i, j] * (s[(i + 1) % n, j] + s[i, (j + 1) % n])
            #C2
            sum_2+= s[i, j] * (s[(i + 2) % n, j] + s[i, (j + 2) % n])
    # average of all configs
    C1=sum_1/(2*n*n)
    C2=sum_2/(2*n*n)
    return C1, C2

def glauber(s, T):
    s_2d = s.copy()
    E_initial = energy_2d(s_2d,J)
    E_glauber= []
    C1=[]
    C2=[]
    M = []  # magnetisation
    n = len(s_2d)
    steps = 1000000
    # generating time series
    for i in range(steps):
        r = random.randrange(n)  # random row
        c = random.randrange(n)  # random column
        s_new = s_2d[r, c] * (-1)  # flipping the spin and storing new spin
        neighbours = s_2d[r, (c - 1) % n] + s_2d[r, (c + 1) % n] + s_2d[(r - 1) % n, c] + s_2d[
            (r + 1) % n, c]  # modulo taken to enforce periodic boundary condition ((num_sites+1)th site=0)
        dE = 2 * J * s_2d[r, c] * neighbours  # change in energy on spin flip
        E_new = E_initial + int(dE)  # new value of energy
        # P=1 if gj<=gi and g(Ej)/g(Ei) if gj>gi =e^ln(g(Ej))/e^ln(g(Ei))=e^(ln(g(Ej))-ln(g(Ei)))
        if np.random.rand() < 1 / (1 + np.exp(dE / T)):
            E_initial = E_new  # accepting flip
            s_2d[r, c] = s_new
        # correlations for edmd observables (just E and M are non-markovian as multiple configs can have same M and E)
        c1,c2 = correlations(s_2d)
        C1.append(c1)
        C2.append(c2)
        #M.append(np.sum(s_2d))
        #E_glauber.append(E_initial)
        #lattice.append(s_2d.copy())
    return np.array(C1), np.array(C2)

def observables(C1,C2):
    #m = M.reshape(-1, 1) / (num_sites**2)   # normalise to [-1,1]
    #e = E.reshape(-1, 1) / (2 * num_sites**2)  # normalise to [-1,1]
    C1=C1.reshape(-1, 1)
    C2=C2.reshape(-1, 1)
    psi = np.hstack([C1,C2])   # shape (n, 5)
    return psi

def edmd(X1,X2,Y1,Y2):
    psi_X = observables(X1,X2)
    psi_Y = observables(Y1,Y2)
    M = X1.shape[0]  # number of snapshots
    print(psi_X.shape)
    G = (psi_X.T @ psi_X) / M  # edmd gram matrix (covariance matrix)
    A = (psi_X.T @ psi_Y) / M  # edmd matrix
    K = np.linalg.lstsq(G + 1e-6 * np.eye(2), A, rcond=None)[0]  # matrix was ill-conditioned and division by zero errors arose. had to use matrix norm
    eigvals = np.linalg.eigvals(K)
    eigvals = np.sort(np.abs(eigvals))[::-1]  # sorted eigenvalues (descending)
    return eigvals


num_sites = 16  # number of lattice sites
n=num_sites
num_spins= num_sites*num_sites
n_2=num_spins
J = 1  # ferromagnetic
kb = 1  # boltzman constant in units

# discrete time evolution of state of the model using Koopman operator (edmd)

# defining parameter (temperature range)
T = np.linspace(2,3,50)
T_relax=[]
gaps_1=[]
gaps_2=[]
relax = []
for k in range(len(T)):
    t = T[k]
    s_2d_1=random_lattice(num_sites)
    s_2d_2=random_lattice(num_sites*2)
    C11,C12= glauber(s_2d_1, t)
    C21,C22= glauber(s_2d_2, t)
    # pairs of successive magnetizations, energies and correlations
    #X1 = M[:-1]
    #Y1 = M[1:]
    #X2= E_glauber[:-1]
    #Y2= E_glauber[1:]
    X11=C11[:-1]
    Y11=C11[1:]
    X12=C12[:-1]
    Y12=C12[1:]
    X21=C21[:-1]
    Y21=C21[1:]
    X22=C22[:-1]
    Y22=C22[1:]
    koop_eigvals_1 = edmd(X11,X12,Y11,Y12)
    koop_eigvals_2 = edmd(X21,X22,Y21,Y22)
    spec_gap_1 = 1 - np.abs(koop_eigvals_1[1])/np.abs(koop_eigvals_1[0])
    spec_gap_2 = 1 - np.abs(koop_eigvals_2[1])/np.abs(koop_eigvals_2[0])
    gaps_1.append(spec_gap_1)
    gaps_2.append(spec_gap_2)
    # finding Koopman operator of discrete-time evolution of DoS
    # time-series pairs such that X[i]=(E[i],M[i]) and Y[i]=(E[i+1],M[i+1]) (state vector)
    # X =np.column_stack([t[:-1],M[:-1]])
    # Y= np.column_stack([t[1:],M[1:]])
    '''
    if spec_gap_1!=0:
        relax.append(1/spec_gap) #relaxation time
        T_relax.append(t)
    print('done ',k)
    '''

#plotting spectral gap
gaps_smooth_1=savgol_filter(gaps_1, 11, 3)
gaps_smooth_2=savgol_filter(gaps_2, 11, 3)
plt.figure()
plt.plot(T, gaps_1, label='16x16 (raw)')
plt.plot(T, gaps_2, 'g', label='32x32 (raw)')
plt.plot(T , gaps_smooth_1, 'r--', label='16x16 (smoothed)')
plt.plot(T , gaps_smooth_2, 'b:', label='32x32 (smoothed)')
plt.axvline(x=2 / np.log(1 + np.sqrt(2)),linestyle='--',label='Tc (Onsager)')
plt.xlabel('Temperature')
plt.ylabel('Spectral Gap ')
plt.grid(True)
plt.title('Koopman Spectral Gap vs Temperature')
plt.legend()
plt.show()

'''
#plotting relaxation time
plt.figure()
plt.plot(T_relax, relax)
plt.axvline(x=2 / np.log(1 + np.sqrt(2)), linestyle='--', label=r'Tc (Onsager)')
plt.xlabel('Temperature')
plt.ylabel('Relaxation Time ')
plt.grid(True)
plt.title('Koopman Relaxation time vs Temperature')
plt.legend()
plt.show()
'''
