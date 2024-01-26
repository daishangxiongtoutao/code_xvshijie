# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:44:37 2015

@author: gome_pa
"""
from __future__ import division
import Tenti_s6_retrieval as s6
import numpy as np
import matplotlib.pyplot as plt

#==============================================================================
#  Thermodynamic properties:
#==============================================================================
m=(1.66053886e-27)*28.97
c_int=1
kb=1.3806503e-23
P=1010e2
wavelength= 532.e-9
k_light=2*np.pi/wavelength
scattering_angle=2.*np.pi/360.*(180)
T=273.  
shear_viscosity=s6.shear_viscosity_func(T)
bulk_viscosity=s6.bulk_viscosity_func(T)
thermal_conductivity=s6.thermal_conductivity_func(T)

#freq=np.linspace(-4e9,4e9,100)
#
#spectrum=s6.S6_spectrum(T,P,shear_viscosity,bulk_viscosity,m,thermal_conductivity,c_int,k_light,freq,scattering_angle,kb)
#
#close("all")
#figure()
#plot(freq,spectrum)



#==============================================================================
#Test with discrete frequency references 
#==============================================================================
N_pix=32
N_freq=1000

FSR_freq=8e9
FSR_px=N_pix

freq_step=FSR_freq/FSR_px

px_array=np.arange(N_pix)
freq_array=np.arange(-N_freq/2,N_freq/2)/N_freq*8e9


def Airy_func(freq,px,FSR,center=0,F=20):
    I=1/(1.+F*np.cos(np.pi*(freq-px*freq_step-center)/FSR)**2)
#    Finesse=pi/(2*arcsin(1/sqrt(F)))
    return I    
    

freq_mesh,px_mesh=np.meshgrid(freq_array,px_array)

instr_matrix= Airy_func(freq_mesh,px_mesh,FSR_freq)   

spectrum=s6.S6_spectrum(T,P,shear_viscosity,bulk_viscosity,m,thermal_conductivity,c_int,k_light,freq_array,scattering_angle,kb)
spectrum=np.zeros(N_freq)
spectrum[N_freq/2]=1





spectrum/=max(spectrum)
measured_spectrum=np.dot(instr_matrix,spectrum)
measured_spectrum/=max(measured_spectrum)

close("all")
figure()
plot(px_array,instr_matrix[:,N_freq/2],linewidth=2,label="Instr. function")
plot(np.arange(N_freq)*N_pix/N_freq,spectrum,linewidth=2,label="Mie ")
plot(px_array,measured_spectrum,"--",linewidth=2,label="PMT signal")
xlabel("px")
ylabel("Intensity [a.u.]")
xlim(0,32)
legend()

figure()
imshow(instr_matrix,aspect="auto",extent=[min(freq_array)/1e9,max(freq_array)/1e9,1,32])
xlabel(r"$\Delta f_j$ [GHz]")
ylabel("Pixel")