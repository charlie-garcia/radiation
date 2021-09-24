#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 11:18:43 2021

@author: carlos
"""
import sys, os
import numpy as np

def radiation(self):
    sigma_force,  W, mean_v = (np.zeros([self.Nf,1], dtype='complex') for i in range(3))
    px0  = self.px
    py0  = self.py

    a0   = np.sin(self.km*px0)
    b0   = np.sin(self.kn*py0)
    
    se   = np.sum(self.A) 

    om   = self.omega.reshape(1,1, self.Nf)  #(self.Nm,self.Nn,1)
    wmn4 = self.wmn.reshape(self.Nm,self.Nn,1)
    A0   =  (1j*4*om*self.F/self.mass * (a0.T*b0).reshape(self.Nm, self.Nn, 1) / (wmn4**2*(1+1j*self.eta) - om**2) ).reshape(self.Nm, self.Nn, self.Nf)

    for iif in range(0,self.Nf):
        vn          = np.sum((self.alpha.T.dot(A0[:,:,iif])) * self.beta.T, 1, keepdims=True)
        mean_v[iif] =  a0.dot(A0[:,:,iif]).dot(b0.T)
        
        v2m = np.mean(np.abs(vn)**2)/2    
        W[iif]            = 1/2* vn.conj().T @ self.R[iif] @ vn     # W = 1/2 * vn'*Q*L*Q'*vn; W = 1/2 * L2'*(abs(y).^2);                                      

        sigma_force[iif]  = W[iif]  / (self.rho0*self.c0*se*v2m)                        # Radiation Efficiency    

    return sigma_force, W, mean_v

def radiation_efficiency(self):
    # Extract variables
    N  = self.N
    M  = self.M
    Nf = self.Nf

    omega = self.omega
    rho0  = self.rho0
    c0    = self.c0
    R     = self.R
    f     = self.f
    se    = np.sum(self.A)    

    sigma_mode  = np.zeros([Nf,1], dtype='complex')

    for iif in range(0,Nf):
        
        vn_modes         = np.array(self.mode.reshape(-1))
        v2m              = np.mean(np.abs(vn_modes)**2)/2    
        W                = 1/2* (vn_modes.conj().T.dot(R[iif])).dot(vn_modes )          # W = 1/2 * vn'*Q*L*Q'*vn; W = 1/2 * L2'*(abs(y).^2);                                      
        sigma_mode[iif]  = W / (rho0*c0*se*v2m)                                         # Radiation Efficiency    

    return sigma_mode

def mode_excitation(self):
    px0  = self.px
    py0  = self.py
    a0   = np.sin(self.km*px0)
    b0   = np.sin(self.kn*py0)
    se   = np.sum(self.A) 

    om   = self.omega.reshape(1,1, self.Nf)  #(self.Nm,self.Nn,1)
    wmn4 = self.wmn.reshape(self.Nm,self.Nn,1)
    A0   =  (1j*4*om*self.F/self.mass * (a0.T*b0).reshape(self.Nm, self.Nn, 1) / (wmn4**2*(1+1j*self.eta) - om**2) ).reshape(self.Nm, self.Nn, self.Nf)

    iif  = self.scroll_frequency.value()
    vn   = np.sum((self.alpha.T.dot(A0[:,:,iif])) * self.beta.T, 1, keepdims=True)

    return vn    