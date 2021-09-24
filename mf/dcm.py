import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore", message="invalid value encountered in true_divide")

def MatrixR(rho,c0,k0,A,rij):
    import scipy.special as sp
    #Compute the Acoustical Resistance Radiation Matrix via Hashimoto (discrete calculation method)
    ai  = np.sqrt(A/np.pi)
    ak  = np.sqrt(A.T/np.pi)
    AA  = A*A.T
    
    R = 2*rho*c0*(k0**2)*AA /(np.pi)* ( sp.jv(1, (k0*ai) )/ (k0*ai))* ( sp.jv(1, (k0*ak) )/ (k0*ak))*  np.sin(k0*rij)/(k0*rij)
    Rii = rho*c0*A* ( 1 - sp.jv(1, (2*k0*ai) ) /(k0*ai) )
    R[np.isnan(R)] = 0
    np.fill_diagonal(R, Rii)
    R = np.real(R)
    return R

def distance3(cs,FieldPoint):
    r = np.zeros([cs.shape[0], FieldPoint.shape[0]], dtype=np.float32)
    for ii in range(FieldPoint.shape[0]):
        r[:,ii] = np.sqrt( np.sum( (cs-FieldPoint[ii,:])**2 ,1) )
    return r

def RadModes(R, n_f, n_modes):
    import scipy.linalg as la            

    N = int(np.sqrt( R.shape[0] ))
    
    ll, Q = la.eig(R)                                          # RADIATION MODES 
    mode_i = [0]*n_modes
    efficiency_i = [0]*n_modes
    
    for irm in range(0, n_modes):
        # mQ[irm] = Q.T[irm,:]
        mode_i[irm] = np.real( Q.T[irm,:] ).reshape(N,N) 
        efficiency_i[irm] = ll[irm]
    return mode_i, efficiency_i

