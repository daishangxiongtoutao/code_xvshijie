# -*- coding: utf-8 -*-
"""
Created on Mon Nov 10 14:32:43 2014

@author: gome_pa
"""
from __future__ import division
import scipy.special 
from numpy import *
#import pandas
#==============================================================================
# Calulation of the Tenti s6 profile. Follorwing  Pan et al. "Coherent Rayleigh-
# Brillouin scattering in molecular gases" (2004)
#==============================================================================
def w_0_func(Z):
    w=pi*exp(-Z**2)*(scipy.special.erfi(Z)-1j)    
    return w
    
    
def S6_spectrum_reduced_variables(xi,y,c_int,R_int,f_u):
    y_7=y
#==============================================================================
# J elements (including n_0/(k*v_0) prefactor)     
#==============================================================================
    J_020=-y
    J_030=J_020
    J_100=-y/R_int*c_int/(3./2. + c_int)
    J_001=3./(2.*c_int)*J_100
    J_001_100=sqrt(3./(2.*c_int))*J_100
    J_100_001=sqrt(3./(2.*c_int))*J_100
    J_110=-2./3.*y-5./6.*(y/R_int)*c_int/(3./2. + c_int)
    J_011_110=-sqrt(5./(8.*c_int))*(y/R_int)*c_int/(3./2. + c_int)
    
    J_011=-2./3.*y*c_int/(3./2. + c_int)* (2./5.*(3./2. + c_int)+(3.+c_int)/(2.*R_int)+9*f_u/(16*R_int**2))/(-1.+4./15.*f_u*(3./2.+c_int)+c_int*f_u/(3.*R_int))
#==============================================================================
# w elements
#==============================================================================
    Z=xi+1j*y_7
    w_0=w_0_func(Z)
    w_1=-sqrt(pi)+Z*w_0
    w_2=Z*w_1
    w_3=-sqrt(pi)/2.+Z*w_2
    w_4=Z*w_3
    w_5=-3./4.*sqrt(pi)+Z*w_4
    w_6=Z*w_5
    
#==============================================================================
# I elements 
#==============================================================================
    I_00_00=1./sqrt(pi)*w_0
    I_02_02=2./(3.*sqrt(pi))*(w_0-2.*w_2+2.*w_4)
    I_00_01=sqrt(2./pi)*w_1
    I_02_10=1./(3.*sqrt(2.*pi))*(w_0+4.*w_2-4.*w_4)
    I_00_02=1./sqrt(3.*pi)*(-w_0+2.*w_2)
    I_00_10=1./sqrt(6.*pi)*(-w_0+2.*w_2)
    I_01_01=2./sqrt(pi)*w_2
    I_02_11=1./sqrt(15.*pi)*(-w_1+8.*w_3-4.*w_5)    
    I_10_10=1./(6*sqrt(pi))*(5.*w_0-4.*w_2+4.*w_4)
    I_00_11=1./sqrt(5.*pi)*(-3.*w_1+2*w_3)
    I_01_02=sqrt(2./(3.*pi))*(-w_1+2.*w_3)
    I_10_11=1./sqrt(30.*pi)*(7.*w_1-8.*w_3+4.*w_5)
    I_01_10=1./sqrt(3.*pi)*(-w_1+2.*w_3)
    I_01_11=sqrt(2./(5.*pi))*(-3.*w_2+2.*w_4)
    I_11_11=1./(5.*sqrt(pi))*(13.*w_2-12.*w_4+4.*w_6)
    
    I_01_00=I_00_01
    I_10_02=I_02_10
    I_02_00=I_00_02
    I_10_00=I_00_10
    I_11_02=I_02_11
    I_11_00=I_00_11
    I_02_01=I_01_02
    I_11_10=I_10_11
    I_10_01=I_01_10
    I_11_01=I_01_11
    
    

#==============================================================================
#Calclates A matrix        
#==============================================================================
    A=array(zeros([len(xi),6,6]),dtype=complex)
    A[:,0,0]=-J_030*I_00_00+1j
    A[:,0,1]=-J_030*I_01_00
    A[:,0,2]=(J_030-J_110)*I_11_00
    A[:,0,3]=(J_030-J_100)*I_10_00
    A[:,0,4]=J_100_001*I_10_00
    A[:,0,5]=J_011_110*I_11_00
    
    A[:,1,0]=-J_030*I_00_01
    A[:,1,1]=-J_030*I_01_01+1j
    A[:,1,2]=(J_030-J_110)*I_11_01 
    A[:,1,3]=(J_030-J_100)*I_10_01
    A[:,1,4]=J_100_001*I_10_01
    A[:,1,5]=J_011_110*I_11_01
    
    A[:,2,0]=-J_030*I_00_11
    A[:,2,1]=-J_030*I_01_11
    A[:,2,2]=(J_030-J_110)*I_11_11-1j
    A[:,2,3]=(J_030-J_100)*I_10_11
    A[:,2,4]=J_100_001*I_10_11
    A[:,2,5]=J_011_110*I_11_11
    
    A[:,3,0]=-J_030*I_00_10
    A[:,3,1]=-J_030*I_01_10
    A[:,3,2]=(J_030-J_110)*I_11_10
    A[:,3,3]=(J_030-J_100)*I_10_10-1j
    A[:,3,4]=J_100_001*I_10_10
    A[:,3,5]=J_011_110*I_11_10
    
    A[:,4,0]=0
    A[:,4,1]=0
    A[:,4,2]=-J_011_110*I_01_00
    A[:,4,3]=-J_100_001*I_00_00
    A[:,4,4]=(J_001-J_030)*I_00_00+1j
    A[:,4,5]=(J_011-J_030)*I_01_00
    
    A[:,5,0]=0
    A[:,5,1]=0
    A[:,5,2]=-J_011_110*I_01_01
    A[:,5,3]=-J_100_001*I_00_01
    A[:,5,4]=(J_001-J_030)*I_00_01
    A[:,5,5]=(J_011-J_030)*I_01_01+1j
    
#==============================================================================
# Generating B column (skipping overall constant term a=1:
#==============================================================================
    B=-array([I_00_00,I_00_01,I_00_11,I_00_10,zeros(len(I_00_00)),zeros(len(I_00_00))],dtype=complex)
    B=B.T
#==============================================================================
# Calculates X
#==============================================================================
    
    X=linalg.solve(A,B)
    return real(X[:,0])
 
 
#==============================================================================
# Tenti S6_model in tems of:
#   T: absolute temperature
#   P: pressure array.
#   shear_viscosity
#   bulk_viscosity
#   m: mass
#   thermal_conductivity
#   c_int: dimensionless internal specific heat capacity.
#   k_light: wavevector of the incident laser light 2*pi/wavelength
#   freq: detuning (frequency difference) with respect to freq_light, where the Tenti S6 spectrum is evaluated.
#   scatering_angle: Angle between incoming lase beam and scattered light (in our config. scatering_angle=pi) 
#   kb: Boltzmann constant 

#==============================================================================
 
def S6_spectrum(T,P,shear_viscosity,bulk_viscosity,m,thermal_conductivity,c_int,k_light,freq,scattering_angle,kb):
    v_0=sqrt(2*kb*T/m) # Thermal velocity
    k=2*k_light*sin(scattering_angle/2.) # Wave-vector of the density wave for elastic scattering.
    y=P/(shear_viscosity*k*v_0) # y parameter
    R_int=3./2.*bulk_viscosity/shear_viscosity*(3./2.+c_int)/c_int
    f_u=m*thermal_conductivity/(shear_viscosity*kb*(3./2.+c_int))
    w=2*pi*freq
    xi=w/(k*v_0)

    return S6_spectrum_reduced_variables(xi,y,c_int,R_int,f_u)
 
 
 
 
    
def bulk_viscosity_func(T):
    eta_b=1.61e-7*T-3.1e-5
    return eta_b
    
def shear_viscosity_func(T):
    eta_0=1.846e-5
    T_0=300.
    T_eta=110.4
    eta=eta_0*(T/T_0)**(3./2.)*(T_0+T_eta)/(T+T_eta)
    return eta

def thermal_conductivity_func(T):   
    kappa_0=26.24e-3
    T_0=300.
    T_A=245.4
    T_B=27.6
    kappa=kappa_0*(T/T_0)**(3./2.)*(T_0+T_A*exp(-T_B/T_0))/(T+T_A*exp(-T_B/T))
    return kappa
    
    
#==============================================================================
# Resize array. If normalization is required it has to be done afterwards.
# Variables:
#   original_array
#   size_new_array: size of the new array
#==============================================================================
def resize_array(original_array,size_new_array):
    size_original_array=len(original_array)
    aux_array=meshgrid(original_array,zeros(size_new_array))[0]
    aux_array=aux_array.T.flatten()
    new_array=zeros(size_new_array)
    for i in range (0,size_new_array):
        new_array[i]=sum(aux_array[i*size_original_array:(i+1)*size_original_array])
    new_array/=double(size_new_array)
    return new_array
