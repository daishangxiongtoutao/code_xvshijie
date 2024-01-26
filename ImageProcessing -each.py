# -*- coding: utf-8 -*-
"""
Created on Fri Aug 07 12:52:40 2015

@author: gome_pa
"""
from __future__ import division
import LoadData
import numpy as np
import os
import binascii
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root, curve_fit
from scipy.interpolate import interp1d
from copy import deepcopy
import Tenti_s6_retrieval as s6
import scipy
import scipy.io as sio
import pandas
#==============================================================================
# Functions
#==============================================================================

def read_radiosonde(filePath):
    f=open(filePath)
    data=f.read()
    f.close()
    
    ind0=data.index("---\n")
    ind0+=4
    ind0=data.index("---\n",ind0)
    ind0+=4
    
    ind1=data.index("\n\n",ind0)
    table=data[ind0:ind1]
    table=table.split("\n")
    
   
    PArray=[] 
    HArray=[]
    TArray=[]
    for line in table:
        try:
            P=float(line[0:7])
            H=float(line[7:14])
            T=float(line[14:21])
            
            PArray+=[P]
            HArray+=[H]
            TArray+=[T]
        except:
            print("Line has been omitted")
    
    ind0=data.find("Station elevation:")
    ind0+=18
    
    ind1=data.find("\n",ind0)
    
    H0=float(data[ind0:ind1]) # Ground altitude (respect to sea level).   

    PArray=1e2*np.array(PArray) #Pressure in Pa
    HArray=np.array(HArray)-H0 #Altitude in meters from GROUND level (for sea level add again the H0)
    TArray=np.array(TArray)+273.15 #Temperature in Kelvin
    
    return [PArray,HArray,TArray, H0]


#def altitude_2_pressure(p0,h):
#    p=p0*np.exp(-h/7990,)
#    return p


def gates_2_altitude(t0,gate):
    t=t0+gate*200.e-9
    h=3e8*t/2.
    return h

def segments_length(xArray,yArray): #Returns the sum of the lengths of a given set of segments
    xArray=np.array(xArray)
    yArray=np.array(yArray)
    lengths=np.sqrt((xArray[1:]-xArray[:-1])**2+(yArray[1:]-yArray[:-1])**2)
    return sum(lengths)
    

def parabola(x,x0,a,b):
    p=a*(x-x0)**2+b
    return p

def peak_finder(data_array,tolerance):
    """Finds the index of all local maxima inside the data array
    
    data_array: array of y data
    tolerance: ratio to be fulfilled between the found peak and the mximum of the array"""
    
    peaks=[]
    max_peak=max(data_array)
    for i in range(1,len(data_array)-1):
        if data_array[i-1]<=data_array[i] and data_array[i]>=data_array[i+1] and data_array[i]/max_peak>=tolerance:
            peaks.append(i)
    return peaks
        
def peak_finder_optimized(data_x,data_y,px_peak,tolerance):
    """Optimzes the peak position using a parabolic fit around the peak
    
    data_x: x coordinate

    data_y: y coordinate and where the local maxima want to be found   
    
    px_peak: number of pixels before and after the peaks and which are included in order to perform the fit function.
    
    Returns:
        x,y,peak   Where "peak" is the iterator found by peak_finder()
    
    """
    peaks=peak_finder(data_y,tolerance)

    peaks_optimized=[]
    
    for i in range(0,len(peaks)):

        try:
            
            peak=peaks[i]
            x=data_x[peak-px_peak:peak+px_peak+1]
            y=data_y[peak-px_peak:peak+px_peak+1]
    
            #print(y)
            x0_guess=data_x[peak]
           #print(x0_guess)
            delta_x=data_x[peak]-data_x[peak-1] #Assumes constant spacing in order to determine the guess value for the second derivative
            a_guess=np.diff(data_y[peak-1:peak+2],n=2)/(delta_x**2)
            #print(a_guess)
            b_guess=data_y[peak]
            #print(b_guess)
            #popt,pcov=curve_fit(parabola,x,y,p0=[x0_guess,a_guess,b_guess])
            #print(popt)
            popt,pcov=curve_fit(parabola,x,y,bounds=([x0_guess-1,-100,0],[x0_guess+1,0,2000]))
            x_max=popt[0]
            y_max=popt[2]
            peaks_optimized.append([x_max,float(y_max),peak])
                     
        except:
            
            print ("peak_finder_optimized ERROR: Fit routine not possible. Local maximum discarted.")
#    print(peaks_optimized)       
    return peaks_optimized


def FWHM(data_x,data_y,px_peak,tolerance):
    print("data")
    print(data_x)
    print(data_y)
    peaks_optimized=peak_finder_optimized(data_x,data_y,px_peak,tolerance)    
    print("data_1")
    print(peaks_optimized)
    return_array=[]
    for i in range (0,len(peaks_optimized)):
        try:
            x_max,y_max,peak= peaks_optimized[i]       
            upper_index=peak
            lower_index=peak
            print("data")
            print(data_y)
            while(data_y[upper_index]>=0.45*y_max):
                upper_index+=1
            while(data_y[lower_index]<=0.45*y_max):
                lower_index-=1
            print("lower_index")
            print(lower_index)
            #Finding the lower point
            print(upper_index)
            x=data_x[lower_index:peak+1]
            y=data_y[lower_index:peak+1]
            f=interp1d(x, y, kind='linear')
            sol=root(lambda x:f(x)-0.5*y_max,data_x[lower_index+1])
            print (sol.message)
            x1=sol.x[0]
            #Finding the upper point
            x=data_x[peak:upper_index+1]
            y=data_y[peak:upper_index+1]
            f=interp1d(x, y, kind='linear')
            sol=root(lambda x:f(x)-0.5*y_max,data_x[upper_index-1])
            print (sol.message)
            x2=sol.x[0]
            FWHM=x2-x1
            return_array.append([x_max,y_max,FWHM,peak,x1,x2])
        except:
            print ("FWHM ERROR: Fit routine not possible. Local maximum discarted.")
            return_array.append([x_max,y_max,None,peak,None,None])
    return return_array


#==============================================================================
# Global parameters
#==============================================================================
ROI_channels=[6,24]
ROI_bins=[0,150]


#==============================================================================
# Load Radiosonde
#==============================================================================

PRadioArray,HRadioArray,TRadioArray,H0=read_radiosonde("Radiosonde.txt")

altitude_2_pressure=interp1d(HRadioArray,PRadioArray, kind=1)
TRadioInterpol=interp1d(HRadioArray,TRadioArray, kind=1)


#==============================================================================
# Illumination
#==============================================================================

FileNameArray=[]
for Root,Dir,Files in os.walk("Illumination"):  
    FileNameArray=Files

dataIllumination=np.zeros([32,670])

for i in range(0,1):
    FullPathFringe=os.path.join("Illumination",FileNameArray[i])
    dataIllumination+=LoadData.LoadLicelData(FullPathFringe)



#==============================================================================
#Cirrus layer 
#==============================================================================
FileNameArray=[]
for Root,Dir,Files in os.walk("Fringe"):  
    FileNameArray=Files

dataRefFringe=np.zeros([32,670])

for i in range(0,1):
    FullPathFringe=os.path.join("Fringe",FileNameArray[i])
    dataRefFringe+=LoadData.LoadLicelData(FullPathFringe)




lowerBound=350
upperBound=450

sumData=sum(dataRefFringe[:,lowerBound:upperBound],0)
indexCirrus=lowerBound+list(sumData).index(max(sumData))

dataCirrusCorrected=dataRefFringe/dataIllumination
dataCirrusCorrected=dataCirrusCorrected[ROI_channels[0]:ROI_channels[1],indexCirrus]

#==============================================================================
# FSR
#==============================================================================
indexRayleigh=120 #Just for determining FSR
width=10
dataRayleighFSR=dataRefFringe/dataIllumination


dataRayleighFSR=np.sum(dataRayleighFSR[ROI_channels[0]:ROI_channels[1],indexRayleigh-width:indexRayleigh+width+1],axis=1)

FSR_analysis= FWHM(np.arange(ROI_channels[0],ROI_channels[1]),dataRayleighFSR,1,0.7)

FSR_guess=FSR_analysis[1][0]-FSR_analysis[0][0]






#def miminizeFunction(FSR):
#    xFull=np.arange(ROI_channels[0],ROI_channels[1])
#    
#    xFull=np.array(xFull,dtype=float)
#    x0=xFull[0]
#    xSorted=deepcopy(xFull)
#    xSorted[xSorted>(x0+FSR)]-=FSR #Folding the second fringe back to a single FSR
#    dataRayleighSorted=np.array([x for (y,x) in sorted(zip(xSorted,dataRayleighFSR))])
##    dataCirrusSorted=np.array([x for (y,x) in sorted(zip(xSorted,dataCirrusCorrected))])
#    xSorted=np.array(sorted(xSorted))
#    
#    return segments_length(xSorted,dataRayleighSorted)
#
#FSR= minimize(miminizeFunction,x0=FSR_guess).x
FSR=FSR_guess


    
xFull=np.arange(ROI_channels[0],ROI_channels[1])
xFull=np.array(xFull,dtype=float)
x0=xFull[0]
xSorted=deepcopy(xFull)
xSorted[xSorted>(x0+FSR)]-=FSR #Translating the second fringe back to a single FSR
   



dataRayleighFSRSorted=np.array([x for (y,x) in sorted(zip(xSorted,dataRayleighFSR))])
dataCirrusSorted=np.array([x for (y,x) in sorted(zip(xSorted,dataCirrusCorrected))])
xSorted=np.array(sorted(xSorted))



#==============================================================================
# Fringes
#==============================================================================
FileNameArray=[]
for Root,Dir,Files in os.walk("Long time series"):  
    FileNameArray=Files

dataFringe=np.zeros([32,670])

#for i in range(0,30):
#FullPathFringe=os.path.join("Long time series",FileNameArray[i])
#dataFringe+=LoadData.LoadLicelData(FullPathFringe)
FullPathRefFringe=os.path.join("Long time series","a1580701.351069")
dataFringe=LoadData.LoadLicelData(FullPathRefFringe)




dataFringeCorrected=dataFringe/dataIllumination
#dataFringeCorrected=dataFringeCorrected[ROI_channels[0]:ROI_channels[1],ROI_bins[0]:ROI_bins[1]]

#==============================================================================
# LOOP:
#==============================================================================

MieArray=[]
MieErrArray=[]
OffsetArray=[]
OffsetErrArray=[]
indexRayleighArray=[]
RayleighMatrix=[]
FitMatrix=[]
scale_factorArray=[]

for indexRayleigh in np.arange(50,150):
    #print (indexRayleigh)


    #==============================================================================
    # Determine complete Instrument Function over exactly 1 FSR
    #==============================================================================
    dataRayleigh=np.sum(dataFringeCorrected[ROI_channels[0]:ROI_channels[1],indexRayleigh:indexRayleigh+1],axis=1)    
    
    xFull=np.arange(ROI_channels[0],ROI_channels[1])
    xFull=np.array(xFull,dtype=float)
    x0=xFull[0]
    xSorted=deepcopy(xFull)
    xSorted[xSorted>(x0+FSR)]-=FSR #Translating the second fringe back to a single FSR
   
#    dataRayleigh=dataRayleigh[int(FSR):]
#    xSorted=xFull[int(FSR):]
#     
    
    
    
    dataRayleighSorted=np.array([x for (y,x) in sorted(zip(xSorted,dataRayleigh))])
    dataCirrusSorted=np.array([x for (y,x) in sorted(zip(xSorted,dataCirrusCorrected))])
    xSorted=np.array(sorted(xSorted))
    
  
    
    
    
#    dataRayleighSorted=dataRayleighSorted[:int(FSR)]
#    dataCirrusSorted=dataCirrusSorted[:int(FSR)]
#    xSorted=xFull[:int(FSR)]   
    
    #==============================================================================
    # Polynomial fit
    #==============================================================================
    
    xExtended=list(xSorted)+list(xSorted+FSR)+list(xSorted+2*FSR)+[xSorted[0]+3*FSR]
    dataCirrusExtended=3*list(dataCirrusSorted)+[dataCirrusSorted[0]]
    f=interp1d(xExtended, dataCirrusExtended, kind=3)
    
    fineResolution=100
    xFine=np.linspace(xExtended[0],xExtended[0]+3*FSR,3*fineResolution+1)
    dataCirrusFine=f(xFine)
    
    
    #==============================================================================
    # Tensi s6 spectrum fit
    #==============================================================================
    FSR_freq=8.387e9
    
    delta_freq=FSR_freq/FSR # freq setps as a function of 
    xFineSymmetric=xFine-xFine[int(len(xFine)/2)]
    freq_array=xFineSymmetric*delta_freq
    
    m=(1.66053886e-27)*28.97
    c_int=1
    kb=1.3806503e-23
    P0=1010e2
    t0=2e-6
    h=gates_2_altitude(t0,indexRayleigh)  
    P=altitude_2_pressure(h)
    T=TRadioInterpol(h)
    wavelength= 532.e-9
    k_light=2*np.pi/wavelength
    scattering_angle=2.*np.pi/360.*(180)
     
    
    
    
    
    
    def fitFunction(freq,offset,scale_factor,mie_intensity):
        T=TRadioInterpol(h)
        shear_viscosity=s6.shear_viscosity_func(T)
        bulk_viscosity=s6.bulk_viscosity_func(T)
        thermal_conductivity=s6.thermal_conductivity_func(T)
        spectrum=s6.S6_spectrum(T,P,shear_viscosity,bulk_viscosity,m,thermal_conductivity,c_int,k_light,freq_array,scattering_angle,kb)
        spectrum/=max(spectrum)    
        spectrum[freq_array==0]+=abs(mie_intensity)
        
        instrumentFunction=dataCirrusFine    
        spectrumConvolved=np.convolve(spectrum,instrumentFunction,mode="same")
        f=interp1d(freq_array,spectrumConvolved, kind=1)    
        return scale_factor*f(freq-offset)
    
    xSymmetric=xSorted-xSorted[int(len(xSorted)/2)]
    freqSorted=(FSR_freq/FSR)*xSymmetric
    offset_guess=1e9
    scale_factor_guess=max(dataRayleighSorted)/fitFunction(0,0,1,0)
    mie_intensity_guess=1
    
    
    
    popt,pcov=scipy.optimize.curve_fit(fitFunction,freqSorted,dataRayleighSorted,p0=[offset_guess,scale_factor_guess,mie_intensity_guess],maxfev=5000)
    [offset,scale_factor,mie_intensity]=popt
    [offset_err,scale_fator_err,mie_intensity_err]=np.sqrt(np.diag(pcov))
    
    
    
    scale_factorArray+=[scale_factor]
    MieArray+=[mie_intensity]
    MieErrArray+=[mie_intensity_err]
    OffsetArray+=[offset]
    OffsetErrArray+=[offset_err]
    indexRayleighArray+=[indexRayleigh]
    RayleighMatrix+=[dataRayleighSorted]
    FitMatrix+=[fitFunction(freqSorted,*popt)]


#==============================================================================
# Sliding average:
#==============================================================================
scale_factorArray=np.array(scale_factorArray)
MieArray=np.array(MieArray)
MieErrArray=np.array(MieErrArray)
OffsetArray=np.array(OffsetArray)
OffsetErrArray=np.array(OffsetErrArray)
indexRayleighArray=np.array(indexRayleighArray)
window_size=1
scale_factorAverage=[]
TErrAverage=[]
indexAverage=[]
MieAverage=[]
MieErrAverage=[]
OffsetAverage=[]
OffsetErrAverage=[]
for i in range (0,len(MieArray)-window_size+1):
    scale_factorAverage+=[np.mean(scale_factorArray[i:i+window_size])]
    indexAverage+=[np.mean(indexRayleighArray[i:i+window_size])]
    MieAverage+=[np.mean(MieArray[i:i+window_size])]
    MieErrAverage+=[np.sqrt(sum(MieErrArray[i:i+window_size]**2))/window_size]
    OffsetAverage+=[np.mean(OffsetArray[i:i+window_size])]
    OffsetErrAverage+=[np.sqrt(sum(OffsetErrArray[i:i+window_size]**2))/window_size]

scale_factorAverage=np.array(scale_factorAverage)
indexAverage=np.array(indexAverage)
MieAverage=np.array(MieAverage)
MieErrAverage=np.array(MieErrAverage)
OffsetAverage=np.array(OffsetAverage)
OffsetErrAverage=np.array(OffsetErrAverage)
#==============================================================================
# Plot routine
#==============================================================================

#plt.close("all")

#fig = plt.figure("dataFringe")
#ax = fig.gca(projection='3d')
#Y = np.arange(32)
#X = np.arange(670)
#X, Y = np.meshgrid(X, Y)
#ax.plot_surface(X, Y, dataFringe , rstride=1, cstride=1, cmap=cm.coolwarm,
#        linewidth=0, antialiased=True)
#
#
#
#fig = plt.figure("dataRefFringe")
#ax = fig.gca(projection='3d')
#Y = np.arange(32)
#X = np.arange(670)
#X, Y = np.meshgrid(X, Y)
#ax.plot_surface(X, Y, dataRefFringe , rstride=1, cstride=1, cmap=cm.coolwarm,
#        linewidth=0, antialiased=True)


#fig = plt.figure("dataFringeCorrected 3D")
#ax = fig.gca(projection='3d')
#Y = np.arange(ROI_channels[0],ROI_channels[1])
#X = np.arange(ROI_bins[0],ROI_bins[1])
#X, Y = np.meshgrid(X, Y)
#ax.plot_surface(X, Y, dataFringeCorrected , rstride=1, cstride=1, cmap=cm.coolwarm,
#        linewidth=0, antialiased=True)

#fig = plt.figure("Corrected Spectra")
#ax = fig.gca()
#ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataCirrusCorrected,"g-")
#ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataCirrusCorrected,"go",label="Cirrus")
#
#ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataRayleigh,"b-")
#ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataRayleigh,"bo",label="Rayleigh")
#
#ax.plot([Rayleigh_analysis[0][0],Rayleigh_analysis[0][0]],[0,Rayleigh_analysis[0][1]],"r-",label="Center")
#ax.plot([Rayleigh_analysis[1][0],Rayleigh_analysis[1][0]],[0,Rayleigh_analysis[1][1]],"r-")
#ax.legend()
#
#
#fig = plt.figure("Sorted Spectra over 1 FSR")
#ax = fig.gca()
#ax.plot(xSorted,dataCirrusSorted,"g-")
#ax.plot(xSorted,dataCirrusSorted,"go",label="Cirrus")

#ax.plot(xSorted,dataRayleighSorted,"b-")
#ax.plot(xSorted,dataRayleighSorted,"bo",label="Rayleigh")
#ax.legend()
#
#
#fig = plt.figure("Fine Cirrus Spectrum")
#ax = fig.gca()
#ax.plot(xExtended,dataCirrusExtended,"go")
#ax.plot(xFine,dataCirrusFine,"g-")
#
#
#
#fig = plt.figure("Fit")
#ax = fig.gca()
#ax.plot(freqSorted,fitFunction(freqSorted,*popt),"r-",label=("T=%.2f +/- %.2f K"%(T,T_err) ))
#ax.plot(freqSorted,dataRayleighSorted,"go")
#ax.legend()




#h=1e-3*gates_2_altitude(2e-6,indexRayleighArray)



#fig = plt.figure("Temperature Profile")
#ax = fig.gca()
#ax.plot(TArray,h,"b-")
#ax.errorbar(TArray, h, xerr=TErrArray, fmt='ro')
#ax.set_xlabel("Temp. (K)")
#ax.set_ylabel("Altitude (km)")
#ax.set_xlim([180,340])



#h=1e-3*gates_2_altitude(2e-6,indexAverage)

#fig = plt.figure("Temperature Profile Average")
#ax = fig.gca()
#ax.errorbar(TAverage, h, xerr=TErrAverage, fmt='r-',label="T-Lidar")
#xlim=ax.get_xlim()
#ylim=ax.get_ylim()

#ax.plot(TRadioArray,1e-3*HRadioArray,"b-",label="Radiosonde")

#yGrid=np.linspace(ylim[0],ylim[1])
#xGrid=-7*(yGrid-ylim[0])
#for dx in np.arange(xlim[0],xlim[1]+7*(ylim[1]-ylim[0]),10):
 #   ax.plot(xGrid+dx,yGrid,"--",color="Black")
    
#ax.set_xlim(xlim)
#ax.set_ylim(ylim)

#ax.set_xlabel("Temp. (K)")
#ax.set_ylabel("Altitude to ground (km)")
#ax.legend()

plt.close("all")

fig = plt.figure("dataFringe")
ax = fig.gca(projection='3d')
Y = np.arange(32)
X = np.arange(670)
X, Y = np.meshgrid(X, Y)
ax.plot_surface(X, Y, dataFringe , rstride=1, cstride=1, cmap=plt.cm.coolwarm,
        linewidth=0, antialiased=True)



fig = plt.figure("dataRefFringe")
ax = fig.gca(projection='3d')
Y = np.arange(32)
X = np.arange(670)
X, Y = np.meshgrid(X, Y)
ax.plot_surface(X, Y, dataRefFringe , rstride=1, cstride=1, cmap=plt.cm.coolwarm,
        linewidth=0, antialiased=True)

fig = plt.figure("dataIllumination")
ax = fig.gca(projection='3d')
Y = np.arange(32)
X = np.arange(670)
X, Y = np.meshgrid(X, Y)
ax.plot_surface(X, Y, dataIllumination , rstride=1, cstride=1, cmap=plt.cm.coolwarm,
        linewidth=0, antialiased=True)
#fig = plt.figure("dataFringeCorrected 3D")
#ax = fig.gca(projection='3d')
#Y = np.arange(ROI_channels[0],ROI_channels[1])
#X = np.arange(ROI_bins[0],ROI_bins[1])
#X, Y = np.meshgrid(X, Y)
#ax.plot_surface(X, Y, dataFringeCorrected , rstride=1, cstride=1, cmap=cm.coolwarm,
#        linewidth=0, antialiased=True)

#fig = plt.figure("Corrected Spectra")
#ax = fig.gca()
#ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataCirrusCorrected,"g-")
#ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataCirrusCorrected,"go",label="Cirrus")
#
#ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataRayleigh,"b-")
#ax.plot(np.arange(ROI_channels[0],ROI_channels[1]),dataRayleigh,"bo",label="Rayleigh")
#
#ax.plot([Rayleigh_analysis[0][0],Rayleigh_analysis[0][0]],[0,Rayleigh_analysis[0][1]],"r-",label="Center")
#ax.plot([Rayleigh_analysis[1][0],Rayleigh_analysis[1][0]],[0,Rayleigh_analysis[1][1]],"r-")
#ax.legend()
#
#
fig = plt.figure("Sorted Spectra over 1 FSR")
ax = fig.gca()
ax.plot(xSorted,dataCirrusSorted,"g-")
ax.plot(xSorted,dataCirrusSorted,"go",label="Cirrus")

ax.plot(xSorted,dataRayleighSorted,"b-")
ax.plot(xSorted,dataRayleighSorted,"bo",label="Rayleigh")
ax.legend()
#
#
fig = plt.figure("Fine Cirrus Spectrum")
ax = fig.gca()
ax.plot(xExtended,dataCirrusExtended,"go")
ax.plot(xFine,dataCirrusFine,"g-")
#
#
#
fig = plt.figure("Fit")
ax = fig.gca()
ax.plot(freqSorted,fitFunction(freqSorted,*popt),"r-")
ax.plot(freqSorted,dataRayleighSorted,"go")
ax.legend()
#plt.show()

h=1e-3*gates_2_altitude(2e-6,indexRayleighArray)
#print(indexRayleighArray)


h=1e-3*gates_2_altitude(2e-6,indexAverage)
fig = plt.figure("Mie Profile Average")
ax = fig.gca()
ax.errorbar(MieAverage, h, xerr=MieErrAverage, fmt='r-',label="T-Lidar")
ax.set_xlabel("Mie component")
ax.set_ylabel("Altitude to ground (km)")
ax.legend()

h=1e-3*gates_2_altitude(2e-6,indexAverage)
fig = plt.figure("Mie Profile Average")
ax = fig.gca()
ax.errorbar(MieAverage*scale_factorAverage, h, xerr=MieErrAverage, fmt='r-',label="T-Lidar")
ax.set_xlabel("Mie component")
ax.set_ylabel("Altitude to ground (km)")
ax.legend()


h=1e-3*gates_2_altitude(2e-6,indexAverage)
fig = plt.figure("Offset Profile Average")
ax = fig.gca()
ax.errorbar(1e-9*OffsetAverage, h, xerr=1e-9*OffsetErrAverage, fmt='r-',label="T-Lidar")
print(1e-9*OffsetAverage)

ax.set_xlabel("Offset (GHz)")
ax.set_ylabel("Altitude to ground (km)")
ax.legend()

plt.show()

#####转存下来照明函数
zhaomingzhengti=[]
zhaomingxishu=np.zeros([32,670])
for i in np.arange(0,670):
    zhao=0
    for j in np.arange(0,32):
        zhao=zhao+dataIllumination[j,i]
    zhaomingzhengti.append(zhao)

for i in np.arange(0,32):
    for j in np.arange(0,670):
        zhaomingxishu[i,j]=dataIllumination[i,j]/zhaomingzhengti[j]
zhaomingxishu1=zhaomingxishu[6:24,:]
sio.savemat('zhaoming.mat',{'zhaomingxishu1': zhaomingxishu1})
sio.savemat('zhaomingzheng.mat',{'zhaomingxishuzheng': dataIllumination})
123
