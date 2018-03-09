# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 11:21:58 2017

@author: maria
"""
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve, convolve_fft
import os
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import pylab

from FUNCTIONS import *
from labelling_damocles_graphs import *
######observed data

os.chdir('/home/maria/Documents/SNR-DATA')
opening = open('OIII-WIFES-UNSNIPPED', 'r')
#opening = open('sn1e0-muse-oiii', 'r')
reading = opening.readlines()

#####################################   
##########################
#############MANIPULATING OBSERVED DATA

list1 = []
list2 = []

    
obsdata = appendline(reading,0,len(reading))
list1, list2 = zip(*obsdata)


obswav = np.array(makeitfloat(list1))
obsflux = np.array(makeitfloat(list2))




#scaled_obsflux = np.array([(i - 3.12e-16) for i in obsflux]) #for casa spectrum

#########converting from wavelength to vel space in observed data
obs_vel = []

c = 2.9979e5
obspeak = 5009.65  #need to figure out whether its shifted!!
labpeak = 5006.842  
 
##OII LINE
#obspeak = 3728.82
#labpeak = 3731.35

zpeak = (obspeak - labpeak)/obspeak
radvel = c * zpeak

for i in obswav:
    
    z = ((i - labpeak)/i)
    j = (c * z) - radvel
    obs_vel.append(j)
    


######BY EYE pick out where the narrow line is and remove it
def snip(array,narrow1,narrow2,size):
    a = np.empty((size*2))
    a[:] = np.nan 
    for index,line in enumerate(array):
        if line == find_nearest(array,narrow1):
            obsflux[index-size : index+size] = a
        if line == find_nearest(array,narrow2):
            obsflux[index-size : index+size] = a
        

wifes = snip(obs_vel,14.1,-2858,6)
#muse = snip(obs_vel,-33,-2864,1)
      
        
      
############################################
###################################
######PROCESSING MODELLED DATA
########################################
###########################################

#search through folders in this file

os.chdir('/home/maria/damocles-master/damocles/output')
#y = ['casagrid_nodust_profile']
z = ['integrated_line_profile.out']


#designed to run on more than one output file
#scales and convolves model data
def changemodel(mylist = []):      ###only thing i need is list of files, x
     
     twovals = []
     packedvels = []
     packedfluxes = []
     for i in mylist:
            modelflux2 = []
            modelflux3 = []
            data = open(i,'r').readlines() # check this works
            model_data = appendline(data,0,len(data))
            
            a,b,c = zip(*model_data)    
            
    
            modelvels = np.array(makeitfloat(b))
            modelflux = np.array(makeitfloat(c))
            
            
            print(np.amax(modelflux))
            scale = np.amax(modelflux)/np.nanmax(obsflux) 
             
            
            for x in modelflux:
                x =  (x/scale)
                modelflux2.append(x)
            
            
            binwidth = modelvels[1] - modelvels[2]
            resolution = 20     #172.3 for muse
            sigma = resolution/(binwidth *2.3548)
    
            g = Gaussian1DKernel(stddev=sigma)
            mod_convolve = convolve(modelflux2,g,boundary = 'extend')
            newscale =  np.amax(mod_convolve)/np.nanmax(obsflux)


            for i in mod_convolve:
                i =  (i * newscale)
                modelflux3.append(i)
                
            
            packedvels.append(modelvels)
      
            packedfluxes.append(modelflux3)
          
     return packedvels,packedfluxes
            
###need a way to search dust.in and species.in and return the dust mass and species  
        
everything = changemodel(z)  


  



####NEED TO DO A CHI SQUARED CALCN THAT DISMISSES THE WORST ONES
#
#
#
#
#




#############################################


#########PLOTTING OBS VS MODEL 

# Labels
plt.xlabel(r"Velocity (km/s)")
plt.ylabel(r"Normalised Flux")
plt.title(title)
#plt.title(r"[OIII] 5009,4959$\AA$") #for runs i want to use in poster

##formatting axis
axes = plt.gca()
axes.set_xlim([-6300,4000])
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))


plt.plot(obs_vel,obsflux,'r-',label='MUSE integrated line profile') 
for i in range(0,len(everything[1])):
    plt.plot(everything[0][i],everything[1][i],label='Damocles model')

plt.legend()
#plt.show(block=True)


os.chdir('/home/maria')

######SAVING EVERYTHING AS AN IMAGE
'''
x = os.listdir('automated-damocles-figures')
if x == []:
    out_name = "1"
    os.chdir('/home/maria/automated-damocles-figures')
    plt.savefig(out_name,dpi=500)
    plt.close()
else:
   
    os.chdir('/home/maria/automated-damocles-figures')
    
    filenumbers = []
    for i in x:
        part1,part2 = i.split('.')
        filenumbers.append(part1)
    
    
    sort_file = sorted(filenumbers, key=int)
    print(sort_file)
    new_no = int(sort_file[-1]) + 1
    out_name = str(new_no)
    plt.savefig(out_name,dpi=500)
plt.close()


###'''