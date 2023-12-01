from qutip import *
import time
import numpy as np
from scipy.interpolate import interp1d,CubicSpline
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline
qutip.settings.num_cpus = 128

def mean_number_photons(delta,input_list,times,psi0):

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
    a = tensor(op_list) #muda conforme o numero de atomos
    
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
    C = np.sqrt(kappa)*a   #cavity  
    C31_list = []
    C32_list = []
    
    for n in range(Nat):
        C31_list.append(np.sqrt(G31)*S13_list[n])
        C32_list.append(np.sqrt(G32)*S23_list[n])
    
    Clist = [C] + C31_list + C32_list

    Ht = D1*S33 + (D1-D2)*S22 + delta*S11 - delta*a.dag()*a + (g*a*S13.dag() + Oc*S23.dag() + EM*a + g*a.dag()*S13 + Oc*S23 + EM*a.dag())
        

    # master equation solver
    #data = mesolve(Ht, psi0, times, Clist, [a.dag() * a])
    #nm_at = data.expect[0]/norm
    #nm_a = np.mean(nm_at[-100:])
        
    # monte carlo solver
    opt = Options(ntraj=2560,num_cpus=128)
    data = mcsolve(Ht, psi0, times, Clist, [a.dag() * a], progress_bar=False, options=opt)
    nm_at = data.expect[0]/norm
    nm_a = np.mean(nm_at[-100:])
    
    #print("%f,%f" %(delta,nm_a)) 

    #steady state method (demands a high amount of memory as we increase the system's dimension)
    #rho_a = steadystate(Ht, Clist, method='iterative-gmres', use_precond=True, use_rcm=True, maxiter=1000, tol=1e-12)
    #nm_a = expect(a.dag()*a, rho_a)/norm  # expectation value of nm -> normalized

    return nm_a


def width_interpol(NMa,delta_list,Oc):
    
    x_in = delta_list
    y_in = NMa
    HMy = 0.5 #(np.max(y)-np.min(y))/2
    
    cs = CubicSpline(x_in,y_in)
    wid = (cs.solve(y=HMy,extrapolate=False)*2)[0]
    
    #saving data
    # Transmission_data = np.vstack((delta_list,NMa))
    # file_data_store('Oc=%.3f.csv'%Oc, Transmission_data.T, numtype="real", numformat="decimal", sep= ",") 
    
    return wid


if __name__ == '__main__':   
    
    ini = time.time() #starts counting time
    print(time.ctime()) #print initial time
    nsteps1 = 300 #number of steps in the parameter list (300)
    Oc_list = np.concatenate([np.linspace(0.05, 1.0, 250),np.linspace(1.0,2.0,50)]) #Rabi frequency list 
    index_j = range(0,len(Oc_list)) #length of the parameter list
    Nat = 1    
    N = 6 #Fock space dimension
    kappa = 1.0 #cavity mode dissipation rate (this is why the parameters should be passed as arguments to the function func1)
    g = (10.0/np.sqrt(Nat))*kappa #atom-field coupling strenght
    EM = np.sqrt(0.001)*kappa #pump field strength
    norm = 4*(EM/(kappa))**2 #normalization
    #decay rates
    G31 = 0.5*kappa  #atom
    G32 = 0.5*kappa  #atom
    #detunings
    D1 = 0.0
    D2 = 0.0
        
    times = np.linspace(0.0, 2100, 1050) #evolution time list
    aux_list = []
    for m in range(Nat+1):
        aux_list.append(basis(3,1))
    aux_list[0] = basis(N,0)
    psi0 = tensor(aux_list)
    
    #create empty lists to store the calculated FWHMs
    Wa = []
    width_a = kappa
    wid_guess = width_a/2

    for j in index_j:
        
        guess_idx = 1
        #transmission curve
        Oc = Oc_list[j]
        input_list = [Nat,N,kappa,g,EM,Oc,norm,G31,G32,D1,D2]
        nsteps = 5  #number of steps for the detuning list

        while width_a == 2*wid_guess:
            
            delta_min = wid_guess-0.01*guess_idx if wid_guess-0.01*guess_idx >= 0 else 0
            delta_max = wid_guess+0.01*guess_idx

            delta_list = np.linspace(delta_min,delta_max,nsteps) #detuning list
            # Obtains "nsteps" points of nm x Delta and calculates the FWHM using a Cubic Spline intepolation
            NMa = []
            for delta in delta_list:
                NMa.append(mean_number_photons(delta, input_list, times, psi0))
            try:
                width_a = width_interpol(NMa, delta_list, Oc)
            except:
                guess_idx += 1

        Wa.append(width_a)
        wid_guess = width_a/2      
        #progress tracking
        fin = time.time() 
        print("%f,%f,%d,%.2f" % (Oc,width_a,j,(fin-ini)/3600))

    output_data = np.vstack((Oc_list, Wa))   # join Oc and expt data
    file_data_store('test.csv', output_data.T, numtype="real", numformat="decimal", sep= ",") 
    


