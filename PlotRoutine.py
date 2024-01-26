# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 19:09:34 2015

@author: gome_pa
"""
import numpy as np
import matplotlib.pyplot as plt
plt.close("all")

#
#fig = plt.figure("Fit")
#ax = fig.gca()
#ax.plot(freqSorted,fitFunction(freqSorted,*popt),"r-",label=("T=%.2f +/- %.2f K"%(T,T_err) ))
#ax.plot(freqSorted,dataRayleighSorted,"go")
#ax.legend()
##for dataRayleighSorted in RayleighMatrix[::10]:
##    ax.plot(freqSorted,dataRayleighSorted,"b-")
#
#fig = plt.figure("Temperature Profile")
#ax = fig.gca()
#ax.plot(TArray,indexRayleighArray)
#
#plt.show()

fig = plt.figure("Illumination")
ax = fig.gca(projection='3d')
Y = np.arange(32)
X = np.arange(670)
X=1e-3*gates_2_altitude(2e-6,X)
X, Y = np.meshgrid(X, Y)
ax.plot_surface(X, Y, dataIllumination , rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=True)
ax.set_xlabel("Altitude to ground (km)")
ax.set_ylabel("Pins")
plt.show()

fig = plt.figure("RefFringe")
ax = fig.gca(projection='3d')
Y = np.arange(32)
X = np.arange(670)
X=1e-3*gates_2_altitude(2e-6,X)
X, Y = np.meshgrid(X, Y)
ax.plot_surface(X, Y, dataRefFringe , rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=True)
ax.set_xlabel("Altitude to ground (km)")
ax.set_ylabel("Pins")
plt.show()


dataRefFringeCorrected=dataRefFringe[ROI_channels[0]:ROI_channels[1],ROI_bins[0]:ROI_bins[1]]/dataIllumination[ROI_channels[0]:ROI_channels[1],ROI_bins[0]:ROI_bins[1]]
fig = plt.figure("RefFringeCorrected")
ax = fig.gca(projection='3d')
Y = np.arange(ROI_channels[0],ROI_channels[1])
X = np.arange(ROI_bins[0],ROI_bins[1])
X=1e-3*gates_2_altitude(2e-6,X)
X, Y = np.meshgrid(X, Y)
ax.plot_surface(X, Y, dataRefFringeCorrected , rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=True)
ax.set_xlabel("Altitude to ground (km)")
ax.set_ylabel("Pins")
plt.show()


fig = plt.figure("Cirrus Signal")
ax = fig.gca()
ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataCirrusCorrected,"g-")
ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataCirrusCorrected,"go",label="Cirrus Signal")
ax.legend()
plt.show()


fig = plt.figure("Rayleigh Signal")
ax = fig.gca()
ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataRayleighFSR,"b-")
ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataRayleighFSR,"bo",label="Rayleigh Signal")
ax.plot([Rayleigh_analysis[0][0],Rayleigh_analysis[0][0]],[0,Rayleigh_analysis[0][1]],"r-",label="Center")
ax.plot([Rayleigh_analysis[1][0],Rayleigh_analysis[1][0]],[0,Rayleigh_analysis[1][1]],"r-")
ax.plot([Rayleigh_analysis[0][0],Rayleigh_analysis[1][0]],[Rayleigh_analysis[1][1]/2,Rayleigh_analysis[1][1]/2],"rx",label="FSR=%.2f"%FSR)
ax.plot([Rayleigh_analysis[0][0],Rayleigh_analysis[1][0]],[Rayleigh_analysis[1][1]/2,Rayleigh_analysis[1][1]/2],"r--")
ax.legend(loc="best")
plt.show()

fig = plt.figure("Cirrus Signal 1FSR")
ax = fig.gca()
ax.plot(xSorted,dataCirrusSorted,"g-")
ax.plot(xSorted,dataCirrusSorted,"go",label="Cirrus Signal")
ax.legend()
plt.show()

fig = plt.figure("Rayleigh Signal 1FSR")
ax = fig.gca()
ax.plot(xSorted,dataRayleighFSRSorted,"b-")
ax.plot(xSorted,dataRayleighFSRSorted,"bo",label="Rayleigh Signal")
ax.legend()
plt.show()

fig = plt.figure("Instrument function")
ax = fig.gca()
ax.plot(xExtended,dataCirrusExtended,"go",label="Instr. func.")
ax.plot(xFine,dataCirrusFine,"g-",label="Oversampled instr. func.")

ax.legend()
plt.show()


fig = plt.figure("Instrument function freq")
ax = fig.gca()

freq_Extended=(xExtended-np.mean(xFine))*FSR_freq/FSR
freq_Fine=(xFine-np.mean(xFine))*FSR_freq/FSR

ax.plot(freq_Extended*1e-9,dataCirrusExtended,"go",label="Instr. func.")
ax.plot(freq_Fine*1e-9,dataCirrusFine,"g-",label="Oversampled instr. func.")
ax.set_xlabel("GHz")
ax.legend()
plt.show()


fig = plt.figure("Tenti S6 spectrum")
ax = fig.gca()
freq_Fine=(xFine-np.mean(xFine))*FSR_freq/FSR
T=273.
P=1010e2
shear_viscosity=s6.shear_viscosity_func(T)
bulk_viscosity=s6.bulk_viscosity_func(T)
thermal_conductivity=s6.thermal_conductivity_func(T)
spectrum=s6.S6_spectrum(T,P,shear_viscosity,bulk_viscosity,m,thermal_conductivity,c_int,k_light,freq_array,scattering_angle,kb)
spectrum/=max(spectrum)    
spectrum[freq_array==0]+=0.25


ax.plot(1e-9*freq_Fine,spectrum,"b-",label="Tenti-s6 spectrum \n(including Mie peak)")
ax.set_xlabel("GHz")
ax.legend()
plt.show()



fig = plt.figure("Convolution")
ax = fig.gca()
freq_Fine=(xFine-np.mean(xFine))*FSR_freq/FSR
T=273.
P=1010e2
shear_viscosity=s6.shear_viscosity_func(T)
bulk_viscosity=s6.bulk_viscosity_func(T)
thermal_conductivity=s6.thermal_conductivity_func(T)
spectrum=s6.S6_spectrum(T,P,shear_viscosity,bulk_viscosity,m,thermal_conductivity,c_int,k_light,freq_array,scattering_angle,kb)
spectrum/=max(spectrum)    
spectrum[freq_array==0]+=0.25

convolution=fitFunction(np.linspace(-4e9,4e9),offset,1,0.25,T)


ax.plot(1e-9*np.linspace(-4e9,4e9),convolution,"r-",label="Convolution")
ax.set_xlabel("GHz")
ax.legend()
plt.show()