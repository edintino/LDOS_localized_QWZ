import numpy as np
from numpy import cos,sin,tan,arccos,pi
import sys
import scipy.sparse as ss
import scipy.sparse.linalg as sl
import time
import pickle
import math

# Pauli matrices:
s0,s1,s2,s3 = np.array([[1.,0.j],[0.,1.]]),np.array([[0.j,1.],[1.,0.]])\
    ,np.array([[0.,-1.j],[1.j,0.]]),np.array([[1.,0.j],[0.,-1.]])

# Important note, print has to be saved into a file from the command line.
# Pickle could be used also, but that precision is not required.

def QWZ_Hamiltonian(LL=200,uu=3.,disorder_weight=0., periodicity=1):
    """ Constructs the Hamiltonian matrix of the QWZ model. 
    CSR sparse matrix is used. Square of size LL x LL, 
    with or without periodic boundary conditions."""
    
    # Simple matrices needed for the construction of the Hamiltonian: 
    mat1 = ss.eye(LL)
    mat2 = ss.diags(np.ones(LL),1,(LL,LL)).tocsr()
    rand_mat = ss.diags(np.random.rand(LL*LL)-0.5)

    # Defines periodicity
    if periodicity == 1:   mat2[-1,0] = 1.

    # Matrices for the hoppings and onsite potentials:
    xhop = ss.kron(mat2,mat1)
    yhop = ss.kron(mat1,mat2)
    diag = ss.eye(LL*LL)

    # Nearest neighbor hopping and onsite matrices of the QWZ model:
    Tx = (s3+1.j*s1)/2.
    Ty = (s3+1.j*s2)/2.
    UU = uu*s3
    WW = disorder_weight*s0
    
    # Construction of the Hamiltonian, using kronecker product:
    HH = ss.kron(xhop, Tx) + ss.kron(yhop, Ty)
    HH += HH.H
    HH += uu*ss.kron(diag, s3) + ss.kron(rand_mat, WW)
    
    return HH

def getmu(rr,HH,num_moments):
    """Constructs the moments required for the Chebyshev polynomial.
    Recursion relation is used to obtain all the required wavefunctions,
    from which twice as many moments can be calculated via recursion
    relations, this decreases the computation time."""

    psi_list = []
    psi = np.matrix(rr).T
    psi_list.append(psi)
    psi = HH*psi
    psi_list.append(psi)

    # First two moments for recursion:
    moment0 = (psi_list[0].H*psi_list[0])[0,0]
    moment1 = (psi_list[0].H*psi_list[1])[0,0]

    # Calculates all the required wavefunctions for the moments using recursion:
    while len(psi_list) < int(num_moments/2):
        psi = 2.*HH*psi_list[len(psi_list)-1] - psi_list[len(psi_list)-2]
        psi_list.append(psi)

    # Data manipulation for recursion:
    psi_list = np.array(psi_list)[:,:,0]
    psi_rolled = psi_list[1:]

    # Moments recursion relation:
    even_moments = list(2*np.sum(psi_list.conj()*psi_list, axis=1) - moment0)
    odd_moments = list(2*np.sum(psi_rolled*psi_list[:-1], axis=1) - moment1)

    moments = [None]*(len(even_moments) + len(odd_moments))

    # Puts together odd and even moments into one list:
    moments[::2] = even_moments
    moments[1::2] = odd_moments
    
    return moments

# Parameters of the QWZ model:
uu = 3.
rseed = 100
LL = 200
disorder_weight = 3.
num_moments = 300

# Passing the parameters from the command line:
nparams = int((len(sys.argv)-1)/2)
for nparam in range(nparams):
    if sys.argv[1+nparam*2] == "-s":
        print("# random seed set to ",sys.argv[2+nparam*2])
        rseed = int(sys.argv[2+nparam*2])
    if sys.argv[1+nparam*2] == "-L":
        print("# system size set to ",sys.argv[2+nparam*2])
        LL = int(sys.argv[2+nparam*2])
    if sys.argv[1+nparam*2] == "-u":
        print("# QWZ onsite set to ",sys.argv[2+nparam*2])
        uu = float(sys.argv[2+nparam*2])
    if sys.argv[1+nparam*2] == "-d":
        print("# disorder weight set to ",sys.argv[2+nparam*2])
        disorder_weight = float(sys.argv[2+nparam*2])
    if sys.argv[1+nparam*2] == "-m":
        print("# number of moments set to ",sys.argv[2+nparam*2])
        num_moments = int(sys.argv[2+nparam*2])


"""LDOS calculation: Calculate Hamiltonian, than shrink it into
the [-1,1] domain for the Chebyshev polynomial, and sample the
energy spectrum for every disorder at the same energies."""

np.random.seed(rseed)
start = time.time()
HH = QWZ_Hamiltonian(LL=LL,uu=uu,disorder_weight=disorder_weight)

"""Imports a correct shrinking factor, I sampled at every integer
disorder values for both topologically trivial a non-trivial case
u=3 and u=1 respectively. Which I extrapolated linearly for decimal
values also, with W in [0.0,13.0]"""
shrink_factor = pickle.load( open( "shrink_factors.p", "rb" ) )
if disorder_weight < 1.0:
    HH /= shrink_factor[(uu,1.0)]*1.05
    Energies = np.linspace(0.02,0.9,100)
else:
    HH /= shrink_factor[(uu,disorder_weight)]*1.05
    Energies = np.linspace(0.02,0.9,75) * shrink_factor[(uu,1.0)] / shrink_factor[(uu,disorder_weight)]
end = time.time()
print("# time taken for hamiltonian= ", end-start, " sec")

# Choose the [0,0] lattice site and calculates the LDOS
start = time.time()
rr = np.zeros(LL*LL*2,dtype=complex)
rr[0] = 1.
rr[1] = 1.
rr /= math.sqrt(sum(rr*rr))
moments = np.array(getmu(rr,HH,num_moments))
moments_len = len(moments)
moments_arange = np.arange(moments_len)
end = time.time()
print("# time taken for moments= ", end-start, " sec")
start = time.time()
LDOS = []
for EE in Energies:
    # Use the jackson kernel to improve convergence
    g_jackson = ((moments_len - moments_arange)* \
                 cos(moments_arange*pi/moments_len) + \
                 sin(moments_arange*pi/moments_len) / tan(pi/moments_len)) \
                 / moments_len 
    # Evaluate the Chebyshev expansion:
    chebt = cos(np.arange(len(moments))*arccos(EE)) * g_jackson
    chebt[0] /= 2.
    LDOS.append(np.dot(chebt,moments))
LDOS = np.array(LDOS).real / pi*np.sqrt(1.-Energies**2)
end = time.time()
print("# time taken for LDOS= ", end-start, " sec")
for densities in LDOS:
    print(densities)
