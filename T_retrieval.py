
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
import xlrd
import xlwt
from xlutils.copy import copy

def Retrieval(dataFringeCorrected,FileName,ROI_channels,FSR,dataCirrusCorrected,altitude_2_pressure,TRadioInterpol,dataFringe,dataRefFringe,TRadioArray,HRadioArray,k):
	#==============================================================================
	# Functions
	#==============================================================================
	def write_excel_xls(path, sheet_name, value):
		index = len(value)
		workbook = xlwt.Workbook()  # 新建一个工作簿
		sheet = workbook.add_sheet(sheet_name)  # 在工作簿中新建一个表格
		for i in range(0, index):
			for j in range(0, len(value[i])):
				sheet.write(i, j, value[i][j])  # 像表格中写入数据（对应的行和列）
		workbook.save(path)  # save 
		print("successfully!")

		
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
		print(peaks)
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
   # LOOP:
   #==============================================================================
	TArray=[]
	TErrArray=[]
	MieArray=[]
	MieErrArray=[]
	OffsetArray=[]
	OffsetErrArray=[]
	baseArray=[]
	baseErrArray=[]
	indexRayleighArray=[]
	RayleighMatrix=[]
	FitMatrix=[]


	for indexRayleigh in np.arange(0,150):
		#print (indexRayleigh)


		#==============================================================================
		# Determine complete Instrument Function over exactly 1 FSR
		#==============================================================================
		dataRayleigh=np.sum(dataFringeCorrected[ROI_channels[0]:ROI_channels[1],indexRayleigh:indexRayleigh+20],axis=1)    
		
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
		
		#plt.figure("Rayleigh spectra")
		#plt.plot(dataRayleighSorted)    
		
		
		
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
		FSR_freq=8.33e9
		
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
		wavelength= 532.e-9
		k_light=2*np.pi/wavelength
		scattering_angle=2.*np.pi/360.*(180)
		 
		
		
		
		
		
		def fitFunction(freq,offset,scale_factor,mie_intensity,T,base):
			shear_viscosity=s6.shear_viscosity_func(T)
			bulk_viscosity=s6.bulk_viscosity_func(T)
			thermal_conductivity=s6.thermal_conductivity_func(T)
			spectrum=s6.S6_spectrum(T,P,shear_viscosity,bulk_viscosity,m,thermal_conductivity,c_int,k_light,freq_array,scattering_angle,kb)
			spectrum/=max(spectrum)    
			spectrum[freq_array==0]+=abs(mie_intensity)
			
			instrumentFunction=dataCirrusFine    
			spectrumConvolved=np.convolve(spectrum,instrumentFunction,mode="same")
			
			f=interp1d(freq_array,spectrumConvolved, kind=1)    
			return scale_factor*f(freq-offset)-base
		
		xSymmetric=xSorted-xSorted[int(len(xSorted)/2)]
		freqSorted=(FSR_freq/FSR)*xSymmetric
		
		T_guess=TRadioInterpol(h)
		offset_guess=1e9
		scale_factor_guess=max(dataRayleighSorted)/fitFunction(0,0,1,0,T_guess,0)
		mie_intensity_guess=1
		base_guess=0
		
		
		popt,pcov=scipy.optimize.curve_fit(fitFunction,freqSorted,dataRayleighSorted,p0=[offset_guess,scale_factor_guess,mie_intensity_guess,T_guess,base_guess],bounds=([-1e14,0,0,200,-1000],[1e14,300,5000,400,1000]))
		[offset,scale_factor,mie_intensity,T,base]=popt
		[offset_err,scale_fator_err,mie_intensity_err,T_err,base_err]=np.sqrt(np.diag(pcov))
		'''name_a="Fit spectrum"
		name_c=str(indexRayleigh)
		name_b='-'
		name_d=str(k)
		name_e="Measurement"
		name_f="50"
		parts = [name_f,name_a,name_b,name_e,name_b,name_d,name_b,name_c,'.png']
		fig_name = ' '.join(parts)
		fig = plt.figure(fig_name)
		ax = fig.gca()
		ax.plot(freqSorted,fitFunction(freqSorted,*popt),"r-",label=("T=%.2f +/- %.2f K"%(T,T_err) ))
		ax.plot(freqSorted,dataRayleighSorted,"go")
		ax.legend()
		plt.savefig(fig_name)'''
		
		
		
		
		TArray+=[T]
		TErrArray+=[T_err]
		MieArray+=[mie_intensity]
		MieErrArray+=[mie_intensity_err]
		OffsetArray+=[offset]
		OffsetErrArray+=[offset_err]
		baseArray+=[base]
		baseErrArray+=[base_err]
		indexRayleighArray+=[indexRayleigh]
		RayleighMatrix+=[dataRayleighSorted]
		FitMatrix+=[fitFunction(freqSorted,*popt)]


	#==============================================================================
	# Sliding average:
	#==============================================================================
	TArray=np.array(TArray)
	TErrArray=np.array(TErrArray)
	MieArray=np.array(MieArray)
	MieErrArray=np.array(MieErrArray)
	OffsetArray=np.array(OffsetArray)
	OffsetErrArray=np.array(OffsetErrArray)
	baseArray=np.array(baseArray)
	baseErrArray=np.array(baseErrArray)
	indexRayleighArray=np.array(indexRayleighArray)
	window_size=10
	TAverage=[]
	TErrAverage=[]
	indexAverage=[]
	MieAverage=[]
	MieErrAverage=[]
	OffsetAverage=[]
	OffsetErrAverage=[]
	baseAverage=[]
	baseErrAverage=[]
	for i in range (0,len(TArray)-window_size+1,window_size):
		
		indexAverage+=[np.mean(indexRayleighArray[i:i+window_size])]
		TAverage+=[np.mean(TArray[i:i+window_size])]
		TErrAverage+=[np.sqrt(sum(TErrArray[i:i+window_size]**2))/window_size]
		MieAverage+=[np.mean(MieArray[i:i+window_size])]
		MieErrAverage+=[np.sqrt(sum(MieErrArray[i:i+window_size]**2))/window_size]
		OffsetAverage+=[np.mean(OffsetArray[i:i+window_size])]
		OffsetErrAverage+=[np.sqrt(sum(OffsetErrArray[i:i+window_size]**2))/window_size]
		baseAverage+=[np.mean(baseArray[i:i+window_size])]
		baseErrAverage+=[np.sqrt(sum(baseErrArray[i:i+window_size]**2))/window_size]
	


	indexAverage=np.array(indexAverage)
	TAverage=np.array(TAverage)
	TErrAverage=np.array(TErrAverage)
	MieAverage=np.array(MieAverage)
	MieErrAverage=np.array(MieErrAverage)
	OffsetAverage=np.array(OffsetAverage)
	OffsetErrAverage=np.array(OffsetErrAverage)
	baseAverage=np.array(baseAverage)
	baseErrAverage=np.array(baseErrAverage)
	h=1e-3*gates_2_altitude(2e-6,indexAverage)

	return [TAverage,TErrAverage,baseAverage,MieAverage,OffsetAverage,h]
