# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 18:20:21 2019

@author: Shashank
"""
# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#filename='CHOH_training'
#channel='CH3_C3H2_C3H3'
#filename='C3H2'
#Run_name='Toulene_400eV'

# fragments=['H+', 'C+', 'I+']
#fragments=['Br79+', 'CHBr279+']
#fragments=['HBr79+', 'CBr279+']  
fragments=[ 'CHBr79+', 'Br279+']
#fragments=['N+', 'N+']
Run_name='Bromoform_140eV'

# fragments=['N+', 'N+']
# Run_name='Air_160eV'

pos_filename = '_pos_000'

channel = str(np.size(fragments))
for i in range(0, np.size(fragments)):
    channel +='_'+fragments[i]

workdir = 'D:/amo/als/2021/simion_py_anbu_Chbr3/'
rep_dir=workdir+'Simion_repository/'

#basedir1='C:/Users/patha/Documents/ALS_19th_March_start/'
channeldir=workdir+Run_name+'/'+channel+'/'
#processed_dir=channeldir+'Processed/'

processed_dir = rep_dir

final_array=[]
for k in range(0,np.size(fragments)):
    filename=fragments[k]
    # filename_all=['vx','vy','tof']
    ### Simion x-> actual z, simion y-> actual x, simion z actual y

    #Load data from the file
    data = np.loadtxt(rep_dir+filename+pos_filename+'.csv',skiprows=8,delimiter=',')
    
    for j in range(0,3):
   
        #The csv file contains two rows for each ion. First row is ion start, second
        # row is for ion splash. For each KE, we iterate with x,y and z direction to 
        # get vx, vy and vz vs x, y and tof fits. 
        
        data_ion_start = data[2*j::6]
        data_ion_splash = data[2*j+1::6]

        #Separate velocity, position and tof details
        vx = data_ion_start[:,9]
        vy = data_ion_start[:,10]
        vz = data_ion_start[:,8]
        KE = data_ion_start[:,11]
        x = data_ion_splash[:,6]
        y = data_ion_splash[:,7]
        tof = data_ion_splash[:,2]*1000   # convert:ng tof from us to ns

    #    vx=np.asarray(vx)
        x=np.asarray(x)
        y=np.asarray(y)
        tof=np.asarray(tof)
    #    to
        if j==0:
            plt.figure(filename+pos_filename+'_Vx')
            plt.plot(x,vx,'o')
            plt.xlabel('X(mm)')
            plt.ylabel('Vx (mm/us)')
            m,b = np.polyfit(x,vx,1)
            if b>0:
                equation = r'$Vx = $' + str(round(m,4))+ r'X' ' + ' + str(round(b,4))
            else:
                equation = r'$Vx = $' + str(round(m,4))+ r'X' + str(round(b,4))
            plt.plot(x,m*x+b, '-',linewidth=2, color='r',label=fragments[k])
            plt.title(equation,fontsize=12)
            plt.legend(frameon=False)
            plt.savefig(rep_dir+filename+pos_filename+'_Vx.png',bbox_inches='tight')
        if j==1:
            plt.figure(filename+pos_filename+'_Vy')
            plt.plot(y,vy,'o')
            plt.xlabel('Y(mm)')
            plt.ylabel('Vy (mm/us)')
            m,b = np.polyfit(y,vy,1)
            if b>0:
                equation = r'$Vx = $' + str(round(m,4))+ r'X' ' + ' + str(round(b,4))
            else:
                equation = r'$Vx = $' + str(round(m,4))+ r'X' + str(round(b,4))
            plt.plot(y,m*y+b, '-',linewidth=2, color='r',label=fragments[k])
            plt.title(equation,fontsize=12)
            plt.legend(frameon=False)
            plt.savefig(rep_dir+filename+pos_filename+'_Vy.png',bbox_inches='tight')
        if j==2:
            plt.figure(filename+pos_filename+'_Vz')
            plt.plot(tof,vz,'o')
            plt.xlabel('TOF(ns)')
            plt.ylabel('Vz (mm/us)')
            m,b = np.polyfit(tof,vz,1)
            if b>0:
                equation = r'$Vx = $' + str(round(m,4))+ r'X' ' + ' + str(round(b,4))
            else:
                equation = r'$Vx = $' + str(round(m,4))+ r'X'  + str(round(b,4))
            plt.plot(tof,m*tof+b, '-',linewidth=2, color='r',label=fragments[k])
            plt.title(equation,fontsize=12)
            plt.legend(frameon=False)
            plt.savefig(rep_dir+filename+pos_filename+'_Vz.png',bbox_inches='tight')
      
        # if (j==0) or (j==2):
        
        #Append the fit slope and intercept for writing to velocity parameter file
        final_array.append(m)
        final_array.append(b)

np.savetxt(processed_dir+'vel_calc_fac'+pos_filename+'_'+channel+'.txt',np.asarray(final_array))
