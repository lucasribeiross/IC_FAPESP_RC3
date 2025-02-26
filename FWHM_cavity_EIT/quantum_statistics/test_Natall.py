from qutip import *
import time
import numpy as np
from scipy.interpolate import interp1d,CubicSpline
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
qutip.settings.num_cpus = 128

# Parameters

kappa = 1.0 #cavity mode dissipation rate (this is why the parameters should be passed as arguments to the function func1)
Oc = 1.0*kappa #Rabi frequency 
delta = 0.5*kappa   # Delta_P
Nat = 3    
N = 12 #Fock space dimension
#g = (1.0/np.sqrt(Nat))*kappa # List atom-field coupling strenght
EM = np.sqrt(1.0)*kappa #pump field strength
norm = 4*(EM/(kappa))**2 #normalization
#decay rates
G31 = 0.5*kappa  #atom
G32 = 0.5*kappa  #atom
#detunings
D1 = 0.0
D2 = 0.0
# coupling strenght
glist = np.linspace(0, 20, 201)*kappa
# time list
t_list = np.linspace(0, 2100, 1050)*kappa


# Field operator
op_list = []
for m in range(Nat+1):
    op_list.append(qeye(3))
op_list[0] = destroy(N)
a = tensor(op_list) 

# Atomic operators
S11_list = []
S22_list = []
S33_list = []
S13_list = []
S23_list = []
    
for n in range(Nat):
    op_list = []
    for m in range(Nat+1):
        op_list.append(qeye(3))

    op_list[0] = qeye(N)
        
    op_list[n+1] = basis(3,0)*basis(3,0).dag()
    S11_list.append(tensor(op_list))

    op_list[n+1] = basis(3,1)*basis(3,1).dag()
    S22_list.append(tensor(op_list))

    op_list[n+1] = basis(3,2)*basis(3,2).dag()
    S33_list.append(tensor(op_list))
        
    op_list[n+1] = basis(3,0)*basis(3,2).dag()
    S13_list.append(tensor(op_list))
        
    op_list[n+1] = basis(3,1)*basis(3,2).dag()
    S23_list.append(tensor(op_list))        
       
#total atomic operators
S11 = 0
S22 = 0
S33 = 0
S13 = 0
S23 = 0
    
for n in range(Nat):
    S11 = S11 + S11_list[n]
    S22 = S22 + S22_list[n]
    S33 = S33 + S33_list[n]
    S13 = S13 + S13_list[n]
    S23 = S23 + S23_list[n]


#Projectors (TALVEZ AQUI POSSA TER UM ERRO)

p1_list = []
p2_list = []
p3_list = []
for m in range(Nat+1):
    p1_list.append(qeye(3))
    p2_list.append(qeye(3))
    p3_list.append(qeye(3))
p1_list[0] = projection(N, 1, 1, offset = None)
p2_list[0] = projection(N, 2, 2, offset = None)
p3_list[0] = projection(N, 3, 3, offset = None)

p1 = 1 * tensor(p1_list)
p2 = 1 * tensor(p2_list)
p3 = 1 * tensor(p3_list)

# Hamiltonian

def calculate_H(g):
    H = D1*S33 + (D1-D2)*S22 + delta*S11 - delta*a.dag()*a + (g*a*S13.dag() + Oc*S23.dag() + EM*a + g*a.dag()*S13 + Oc*S23 + EM*a.dag())
    return H

# Collapse operators

C = np.sqrt(kappa)*a   #cavity  
C31_list = []
C32_list = []
    
for n in range(Nat):
    C31_list.append(np.sqrt(G31)*S13_list[n])
    C32_list.append(np.sqrt(G32)*S23_list[n])
    
Clist = [C] + C31_list + C32_list

# Initial state

aux_list = []
for m in range(Nat+1):
    aux_list.append(basis(3,1))
aux_list[0] = basis(N,0)
psi0 = tensor(aux_list)

# Create empty lists

nc1_list = []      #Mean number of photon inside the cavity
nc2_list = []      #to calculate g2
n11_list = []       #Mean number of occupation in excited state |1>
n22_list = []       #Mean number of occupation in excited state |2>
n33_list = []       #Mean number of occupation in excited state |3>
g2_0_list = []     #Second-order correlation function 
p_1_list = []
p_2_list = []
p_3_list = []

# Parameters for MC

arg = [a.dag()*a, (a.dag()**2)*(a**2), S11, S22, S33, p1, p2, p3]
opt = Options(ntraj=2560,num_cpus=128)

for _g in glist:
    print(_g)
    
    data = mcsolve(calculate_H(_g), psi0, t_list, Clist, arg, progress_bar=False, options=opt)
    
    nc1 = data.expect[0]/norm
    nc1_m = np.mean(nc1[-100:])
    nc2 = data.expect[1]/norm
    nc2_m = np.mean(nc2[-100:])
    
    g2_0 = nc2_m / nc1_m ** 2
    
    # results lists
    nc1_list.append(nc1_m)
    n11_list.append(data.expect[2])
    n22_list.append(data.expect[3])
    n33_list.append(data.expect[4])
    g2_0_list.append(g2_0)

    p_1_list.append(data.expect[5])
    p_2_list.append(data.expect[6])
    p_3_list.append(data.expect[7])
    

# Saving

output_data = np.vstack((glist, nc1_list, n11_list, n22_list, n33_list, g2_0_list, p_1_list, p_2_list, p_3_list))   
file_data_store('MeanNumberxg_Nat3_ep1_Oc1_Dp05_MC.csv', output_data.T, numtype="real", numformat="decimal", sep= ",") 