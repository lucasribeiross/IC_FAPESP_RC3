from qutip import *
import time
import numpy as np
from scipy.interpolate import interp1d,CubicSpline
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
qutip.settings.num_cpus = 128

def mean_number_photons(Dp,input_list,times,psi0):

    """
    This function handles the time evolution and returns the mean number 
    of photons for each value of detuning.
    """
    
    Nat = input_list[0]
    N = input_list[1] 
    kappa = input_list[2]
    g = input_list[3]
    EM = input_list[4]
    Oc = input_list[5]
    norm = input_list[6]
    
    G31 = input_list[7]
    G32 = input_list[8]
    D1 = input_list[9]
    D2 = input_list[10]
    
    #field operator
    op_list = []
    for m in range(Nat+1):
        op_list.append(qeye(3))
    op_list[0] = destroy(N)
    a = tensor(op_list) #muda conforme o numero de atoomos
    
    #atomic operators
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

    #colapse Operators
    C = np.sqrt(2*kappa)*a   #cavity  
    C31_list = []
    C32_list = []
    
    for n in range(Nat):
        C31_list.append(np.sqrt(G31)*S13_list[n])
        C32_list.append(np.sqrt(G32)*S23_list[n])
    
    Clist = [C] + C31_list + C32_list

    Ht = D1*S33 + (D1-D2)*S22 + Dp*S11 - Dp*a.dag()*a + (g*a*S13.dag() + Oc*S23.dag() + EM*a + g*a.dag()*S13 + Oc*S23 + EM*a.dag())
        
    # monte carlo solver
    opt = Options(ntraj=2560,num_cpus=128)
    data = mcsolve(Ht, psi0, times, Clist, [a.dag() * a, (a.dag()**2)*(a**2)], progress_bar=False, options=opt)
    nm_at = data.expect[0]/norm
    nm_a = np.mean(nm_at[-100:])
    
    
    #Calculating the expect values (mean number of photons):
    n_1_pre = data.expect[0]
    n_1 = np.mean(n_1_pre[-100:])
    n_2_pre = data.expect[1] 
    n_2 = np.mean(n_2_pre[-100:])
    
    #Correlation function (g^(2)(0) = a.dag()*a.dag()*a*a/(a.dag()*a)^2):

    g2_0 = (n_2)/(n_1**2) 

    
    #print("%f,%f" %(delta,nm_a)) 

    #steady state method (demands a high amount of memory as we increase the system's dimension)
    #rho_a = steadystate(Ht, Clist, method='direct', use_precond=False, use_rcm=True, maxiter=1000, tol=1e-9)
    #nm_a = expect(a.dag()*a, rho_a)/norm  # expectation value of nm -> normalized

    
    #Calculating the expect values (mean number of photons):
        
    #n_1 = nm_a*norm
    #n_2 = expect((a.dag()**2)*(a**2), rho_a) #/norm #norm**2 
    
    #Correlation function (g^(2)(0) = a.dag()*a.dag()*a*a/(a.dag()*a)^2):

    #g2_0 = (n_2)/(n_1**2) 

    
    return [n_1, g2_0]


if __name__ == '__main__':   
    
    ini = time.time() #starts counting time
    print(time.ctime()) #print initial time
    
    Nat = 3    
    N = 6 #Fock space dimension
    kappa = 1.0 #cavity mode dissipation rate (this is why the parameters should be passed as arguments to the function func1)
    nsteps1 = 300 #number of steps in the parameter list (300)
    Oc = 1.0*kappa   #rabi frequency of control field
    #g = (1.0/np.sqrt(Nat))*kappa #atom-field coupling strenght
    g = 10.0*kappa #retirando a renormalização
    EM = np.sqrt(0.1)*kappa #pump field strength
    norm = 4*(EM/(kappa))**2 #normalization
    #decay rates
    G31 = 0.5*kappa  #atom
    G32 = 0.5*kappa  #atom
    #detunings
    D1 = 0.0
    D2 = 0.0
    Dp_list = np.linspace(-2.5,2.5,300)/kappa          #detunig: (omega-omega_P) -> resonance (max transmission)
    index_j = range(0,len(Dp_list)) #length of the parameter list
        
    times = np.linspace(0.0, 2100, 1050) #evolution time list
    aux_list = []
    for m in range(Nat+1):
        aux_list.append(basis(3,1))
    aux_list[0] = basis(N,0)
    psi0 = tensor(aux_list)
    
    #Loop for the values of DeltaP and calculate g2(0):
    
    n_1_list = []
    g2_0_list = []

    for j in index_j:
        
        #transmission curve
        Dp = Dp_list[j]
        
        input_list = [Nat,N,kappa,g,EM,Oc,norm,G31,G32,D1,D2]
        
        results = mean_number_photons(Dp, input_list, times, psi0)
        
        n_1 = results[0] 
        n_1_list.append(n_1)
             
        g2_0 = results[1] 
        g2_0_list.append(g2_0)  
     
        #progress tracking
        fin = time.time() 
        print("%f,%f,%f,%d,%.2f" % (Dp,n_1,g2_0,j,(fin-ini)/3600))

    output_data = np.vstack((Dp_list, n_1_list, g2_0_list))   # join Oc and expt data
    file_data_store('g2_0_x_Dp_and_n1.csv', output_data.T, numtype="real", numformat="decimal", sep= ",") 
    
