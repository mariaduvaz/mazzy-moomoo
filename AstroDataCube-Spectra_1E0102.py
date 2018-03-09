# ==========================================================================================================
# AstroDataCube-Spectra_1E0102.py
# ==========================================================================================================
#
# AUTHOR:           F. Schmidt (f.schmidt.16@ucl.ac.uk)
# CREATION DATE:    25.10.2017
# PROJECT:          WiFeS
# DESCRIPTION:      Python script to extract data and plot spctra from AstroDataCube 1E0102
# NOTES:            -
#
# HISTORY:          25.10.2017: Script creation (F.Schmidt)
#
# ==========================================================================================================



# Load libraries
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
from FUNCTIONS import *
import os
from mpl_toolkits.mplot3d import Axes3D
# +--------------------------------------------------------------------------------------------------------+
# |                                            FITS File Operations                                        |
# +--------------------------------------------------------------------------------------------------------+


# Load the FITS file
os.chdir('/home/maria/Documents/SNR-DATA/')
hdulist = fits.open('mosB_1E0102_wcs_lowess_frac_0.6_it_5_removed_deblend.fits')

# Extract the data and store it in array
scidata = hdulist[0].data




# Define the dimensions
slabs = len(scidata)                    # bin no,   'no of data points in each spectrum'

yaxis = len(scidata[0])               #56, y AXIS    1st dim     

xaxis = len(scidata[0][0])            #48, x AXIS   2nd dim

print(xaxis,yaxis)

###constant for conversions
age = 1020 * (365)#*24*60*60) #age of SNR in seconds
print(age)
smc_dist = 59000 #distance to smc in parsec
pc = 3.086e13 #km in a parsec
arcsec_to_km = (np.pi * smc_dist * pc)/(180*3600)

phys_size = 12 * arcsec_to_km
print(phys_size)

age_s = phys_size/2150     #in seconds
age_years = age_s/(365*24*60*60)
print(age_years*(365))


c = 2.9979e5     #in km
obspeak = 5009.65  #need to figure out whether its shifted!!
labpeak = 5006.842 

zpeak = (obspeak - labpeak)/obspeak
radvel = c * zpeak



#########
#make x-axis

lambda_start = 4200.0                           #from vogt et al. 2017
lambda_stop = 5548.0
lambda_step = (lambda_stop - lambda_start)/slabs

obswav = np.arange(lambda_start, lambda_stop, lambda_step)

obs_vel = []

##going from wavelength on spectral axis to velocity

for i in obswav:
    
    z = ((i - labpeak)/i)
    j = (c * z) - radvel
    obs_vel.append(j)



###from obs vel observationally we want to see what velocity values we want to
#mask out. then we save an array with the indices


for index,line in enumerate(obs_vel):
    if line == find_nearest(obs_vel,-4000): #cutting off rest of spect #MAKE THESE LARGER WHEN IM PLOTTN GRID
        start = index
        
    if line == find_nearest(obs_vel,-177):
        stop = index    
        
    if line == find_nearest(obs_vel,194):  #taking out narrow lines
        start2 = index  
        
    if line == find_nearest(obs_vel,4000):
        stop2 = index

    
    
        
ind = np.arange(start,stop,1) 
ind_2 = np.arange(start2,stop2,1)

oiii_line_indices = np.concatenate((ind,ind_2)) ##all relevant information is in this array

###resetting velocity axis so its the same size as flux array
vel_axis = [obs_vel[i] for i in oiii_line_indices]




# Initialize list to store spectrum for each pixel

spectrum_pixel = []
xvals = []
yvals = []
comb = []
total_pixel = []

#Loop over all pixels in x dimension, 56 ti(mes
for i in range(0, (yaxis)):
    
    # Loop ober all pixels in y direction, 48 times
    for j in range(0, (xaxis)):
        
        comb.append([j + 1,i + 1])    #matches way pixels are read in in image (checked with QFITSview)
       
        conv_km_j = round(((j - 19)*arcsec_to_km)/age,1)
        conv_km_i = round(((i - 27)*arcsec_to_km)/age,1)
        
         #changing it so centre of expansion is point '0-0'
        
        xvals.append(conv_km_j)
        yvals.append(conv_km_i)
        # Generate temporary array
        tmp_list1 = []
        tmp_list2 = []
        
        # Loop over all slapsa
        for k in oiii_line_indices:   #for grid
        
            tmp_list1.append(scidata[k][i][j])
            
        for k in range(0,slabs):  #for integrated spectra
            
            tmp_list2.append(scidata[k][i][j])
            
        # Append tmp_list to spectrum
        spectrum_pixel.append(tmp_list1)
        total_pixel.append(tmp_list2)


###locating peak of velocity...
#creating a velocity array for every spect_pixel value, as they are all going to be different sizes
vel_axis_grid = [vel_axis] * len(spectrum_pixel)



###getting peak >0 and <0 of 1361
#for i in vel_axis_grid[1361]:
  #  if i > 0:
   #     print(i)





####

#HERE I WANT TO SELECT THE RELEVANT X,Y POINTS FOR THE SPECTRa with OIII emission
#ALSO WANT TO TRIM THESE DOWN SO WE ONLY HAVE THE LINE OF INTEREST
#SO WE HAVE ONE LOOP WHERE WE TRIM OUT SPECTRA AND THEIR ASSOCIATED POINTS IN SPACE
#AND ALSO TRIMMING THESE DOWN AROUND THE PEAK

x_neb = []
y_neb = []
x = []
y = []
oiii_spectrum = []
vel_emission = []



for i in range(0,len(xvals)): #loop over every spectra
    if sum(spectrum_pixel[i]) > 6e-16:
        
        
        
        peak = np.amax((np.absolute(spectrum_pixel[i])))
        for index,line in enumerate(spectrum_pixel[i]):  #loop over every entry in every spectra
            
            if line == peak:
              #print(i)
              x_neb.append(xvals[i])
              y_neb.append(yvals[i])
              #find peak
              oiii_spectrum.append(spectrum_pixel[i][index - 10: index + 10])
              vel_emission.append(vel_axis_grid[i][index - 10: index + 10])
###trimmed relevant spectra, only relevant ones will be read in
                
                

#print(np.shape(oiii_spectrum))




datagrid = ([[0,0,0,0]] * len(oiii_spectrum) * len(oiii_spectrum[0]))
##ALSO NEEDS TO MATCH DIMENSIONS, OF EVERY DIFFRENT SPECT AND VEL AXIS FOR EERY SPECTRUM
print(len(datagrid))



#writing out whole grid

counter = 20
for i in range(0,len(oiii_spectrum)):
    
         #DIFFERENT FOR EACH ONE
         #spect length
        
        x_index = [x_neb[i]] * 20
        y_index = [y_neb[i]] * 20
        
       
       
        spect_chunk = (list(zip(x_index,y_index,vel_emission[i],oiii_spectrum[i])))
        print(spect_chunk)
        datagrid[(counter - 20) : counter] = spect_chunk
        counter = counter + 20
        



os.chdir('/home/maria/Documents')
with open('e0102_mapped_datagrid', 'w') as fp:
     fp.write('\n'.join('{} {} {} {}'.format(x[0],x[1],x[2],x[3]) for x in datagrid))



##############checking how grid reads out - 3d plotting
list1,list2,list3,list4 = list(zip(*datagrid)) 
v_x = (makeitfloat(list1))
v_y = (makeitfloat(list2))
v_z = (makeitfloat(list3))
density = makeitfloat(list4)


fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')

ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')


plt.title("0 rot")
ax.scatter(v_x,v_y,v_z,c=density,cmap='spring',label='Wifes High density')


'''

#####   
####if there are any empty spectra full of nans, then we're taking note of these so we don't plot the
#this is because we're plotting entire spectrum so we need to remove then so they dont affect the summation
#might not actually need any of this
nan_cube = []
for i in range(0,len(total_pixel)): #looping over every spectra
    tmp_list = []

    for j in range(0,len(total_pixel[0])):    #looping over every intensity point                 
        x = np.isnan(total_pixel[i][j])   
        tmp_list.append(x)
    
    if set(tmp_list) == {True}:
        nan_cube.append(i)
   

spect_minusnan = []

for i in range(0,len(total_pixel)):
        if i not in nan_cube:
            spect_minusnan.append(total_pixel[i])




integ_spect = []
for i in zip(*spect_minusnan):

    integ_spect.append(sum(i))




# Close the FITS file
hdulist.close()

# +--------------------------------------------------------------------------------------------------------+
# |                                                   Plots                                                |
# +--------------------------------------------------------------------------------------------------------+




# Plot the individual spectra for each pixel
# ------------------------------------------
fig = plt.figure()
pl = fig.add_subplot(111)

# Axis labels
pl.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

# Labels
#pl.set_title(r"WiFeS 1E0102".format(i+1))
pl.set_xlabel(r"Wavelength ($\AA$)")
pl.set_ylabel(r"Flux")

# Plot data
plt.plot(vel_axis_grid[1361],spectrum_pixel[1361], color="darkred",linewidth=0.7)


# Plot the spectrum for the entire nebula
# ---------------------------------------
        
# Initialise figure
fig = plt.figure()
pl = fig.add_subplot(111)

# Axis labels
pl.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))

# Labels
#pl.set_title(r"WiFeS 1E0102".format(i+1))
pl.set_xlabel(r"Wavelength ($\AA$)")
pl.set_ylabel(r"Flux")

# Plot data
plt.plot(obs_vel,integ_spect, color="darkred",linewidth=0.7)
#x,y = list(zip(*spatial_pixel))
#plt.scatter(x,y)
plt.tight_layout()
plt.show()


integ_spect_data = list(zip(vel_axis,integ_spect))

os.chdir('/home/maria/damocles-master/damocles/input')
np.savetxt('integrated_oiii_e0102.in',integ_spect_data,fmt="%d %d")


# Save .png file
#out_name = "1E0102.png".format(num=i+1)
#fig.savefig(out_name, dpi=500)
#plt.close(fig)
##'''