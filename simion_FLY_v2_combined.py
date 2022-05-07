# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 17:36:53 2019

@author: Shashank Pathak (credit: Razib Obaid)

Modified: 2021 - Anbu
Changed to create single fly file for vx, vy and vz calibration.
To create fly file for ToF calibration.

"""

import numpy as np
### Simion x-> actual z, simion y-> actual x, simion z actual y
def flywrite_vel_cal(filename,directory,mass, charge, position):
    pos = '_pos_'+str(position[0])+str(position[1])+str(position[2])
    # filename_all=[filename+pos+'_vx',filename+pos+'_vy',filename+pos+'_tof']
    filename_all=filename+pos
    
    #Open File
    f = open(directory+filename_all+'.fly2', 'w')
    
    # if j==0:
    #     vec='(0,1,0)'
    # elif j==1:
    #     vec='(0,0,1)'
    # elif j==2:
    #     vec='(1,0,0)'
    #Writing the Simion FLY file from here
    foreText = "particles {coordinates = 0,"
    f.write(foreText)

    # First loop for energy range. Default 0.1 eV to 10 eV
    for i in np.arange(0.01, 10.1, 0.1):
        # For each KE step, cycle through x,y and z directions
        for vec in ('(0,1,0)', '(0,0,1)', '(1,0,0)'):
            midText = "\n standard_beam {n = 1,  tob = 0,  mass = "+str(mass)+",  charge = "+str(charge)
            midText = midText + ", ke ="+ str(np.round(i,2)) + ",  cwf = 1,  color = 0,  direction = vector"+vec
            midText = midText + ",  position = vector"+str(position)+"},"
        
            f.write(midText)

    f.write( "}")
    f.close()
 
# Write tof calibration with starting position
def flywrite_tof_cal(filename,directory,mass,charge,position):
    # filename_all=[filename+'_vx',filename+'_vy',filename+'_tof']
    filename_all = filename+'_pos_'+str(position[0])+str(position[1])+str(position[2])
    
    f = open(directory+filename_all+'.fly2', 'w')

    #Writing the Simion FLY file from here
    foreText = "particles {coordinates = 0,"
    f.write(foreText)
    
    for m in mass:    
        midText = "\n standard_beam {n = 1,  tob = 0,  mass = "+str(m)+",  charge = "+str(charge)
        midText = midText + ", ke =0,  cwf = 1,  color = 0,  direction = vector(0, 0, 0)"
        midText = midText + ",  position = vector"+str(position)+"},"
        
        f.write(midText)

    f.write( "}")
    f.close()    
    
###++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
###DO NOT CHANGE ANYTHING ABOVE THIS
####+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
workdir = 'D:/amo/als/2021/simion_py_anbu_Chbr3/'
rep_dir=workdir+'Simion_repository/'

###### Write tof calibration fly file
mass = [1.00782503223, 13.00782503223, 63.45223595, 126.90447190000, 141.927947]
position_all = [(0,0,0), (0,1,0)]#, (0,0,1), (0,1,1)]
# position_all = [(0,2,0), (0,0,2), (0,2,2)]
# position_all = [(0,5,0), (0,0,5), (0,5,5)]

# for position in position_all:
#     flywrite_tof_cal('tof_cal',rep_dir,mass,1,position) 


###### Write velocity calibration fly file
# flywrite_vel_cal('C+',rep_dir,12.0107, 1)
# flywrite_vel_cal('CH+',rep_dir,13.0186, 1)
# flywrite_vel_cal('CH2+',rep_dir,14.0266, 1)
# flywrite_vel_cal('CH3+',rep_dir,15.0345, 1)
# flywrite_vel_cal('I++',rep_dir,126.90447, 2) 
# flywrite_vel_cal('H+',rep_dir,1.00794, 1)
# flywrite_vel_cal('Ne+',rep_dir,20.1797, 1) 

for position in position_all:

    # flywrite_vel_cal('Br79+',rep_dir,78.918, 1, position)
    # flywrite_vel_cal('Br81+',rep_dir,80.916, 1, position)
    # flywrite_vel_cal('HBr79+',rep_dir,79.926, 1, position)
    # flywrite_vel_cal('HBr81+',rep_dir,81.924, 1, position)
    # flywrite_vel_cal('CHBr79+',rep_dir,91.926, 1, position)
    # flywrite_vel_cal('CHBr81+',rep_dir,93.924, 1, position)
    # flywrite_vel_cal('CHBr279+',rep_dir,170.884, 1, position)
    # flywrite_vel_cal('CHBr281+',rep_dir,174.840, 1, position)
     #flywrite_vel_cal('CBr279+',rep_dir,169.836, 1, position)
     #flywrite_vel_cal('CBr281+',rep_dir,173.832, 1, position)
   flywrite_vel_cal('Br279+',rep_dir,157.836, 1, position)

