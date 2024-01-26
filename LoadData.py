# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 12:18:55 2015

@author: gome_pa

"""

import numpy as np 
import struct as st
    
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

#path="F:\\WM\\PMT32Channel Main Folder\\"
#
#filename="a1561811.230048"
#hex_2_month={"1":"01","2":"02","3":"03","4":"04","5":"05","6":"06","7":"07","8":"08","9":"09","A":"10","B":"11","C":"12"}
#Timestamp="20"+filename[1:3]+"/"+hex_2_month[filename[3:4]]+"/"+filename[4:6]+" "+filename[6:8]+":"+filename[9:11]+":"+filename[11:13]
#print Timestamp
#
#f=open(path+filename)
#data=f.read()
#f.close()
#start= data.index("\n\n")
#start+=2
#data=data[start:-1]
#data=data.split("\n")
#
#data_matrix=np.zeros([32,670])
#for i in range(0,32):
#     for j in range (0,670):
#         data_matrix[i,j]=st.unpack("<L",data[i][j*4:(j+1)*4])[0] # data is stored in little-endian (LSB first) format of 4 bytes
#
#
#px_array=arange(32)
#z_array=arange(670)*30#
#close("all")
#imshow(np.flipud(data_matrix.T),extent=[0,32,0,30*670],aspect="auto",interpolation="none",cmap="jet")
#colorbar()




def LoadLicelData(FullFileName):
    f=open(FullFileName,"rb")
    data=f.read()
#    return data
    f.close()
    index0= data.index(b"\r\n\r\n")
    index0+=4 #Shifts initial "\r\n\r\n"
#
#    data=data[start:-2]
#    print (len(data)-2*31)/(4.*670*32)
#    data=data.split("\r\n")
    
    
    
    data_matrix=np.zeros([32,670])
    for i in range(0,32):
        for j in range (0,670):
            index=index0+((4*670+2)*i)+4*j  # Remember: each channel has 670 times 4 byte integers. And the channels are separated by "\r\n" (2 byte)
            data_matrix[i,j]=st.unpack("<L",data[index:index+4])[0] # data is stored in little-endian (LSB first) format of 4 bytes    
             
    return data_matrix
    
    
    
    
    


