import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from osgeo import gdal
from skimage import filters
import pandas as pd

import glob
import seaborn as sns
import sys
sys.path.append('/path/to/adjacencyanalysis.py')
import adjacencyanalysis as aa

custom_params = {'figure.constrained_layout.use' : True}
sns.set_theme(style='ticks', palette='Spectral_r', font='sans-serif', font_scale=1.2, color_codes=True, rc=custom_params)


wavelength=5
totalIF=[0.0150,0.0325]
fracs1=[35,70]
fracs01=[2,10]
fracs101=[0,3]
fracs011=[0,2]
fracs010=[0.02,0.18]

figname=scenariodirectoryname+'elephantpanel_'+scenariodirectoryname[0:-1]+'_'+str(wavelength)+'um.pdf'
totalcubefile=scenariodirectoryname+elephantnameprefix+str(wavelength)+'.colorCCD.xml'

print('image',totalcubefile)
#print('Final fig name = ',figname)

elephantcubes=['_____1.Jcube.colorCCD.xml','____01.Jcube.colorCCD.xml',
               '___101.Jcube.colorCCD.xml','___011.Jcube.colorCCD.xml','___010.Jcube.colorCCD.xml']

sns.set_style({'axes.grid' : False})

wavefig,[ax,ax1,ax01,ax101,ax011,ax010]=plt.subplots(1,6,figsize=(9,16),layout="constrained")

CUBE = aa.getcubefromPDS4(totalcubefile) 
CUBE=CUBE[2:CUBE.shape[0]-2,2:CUBE.shape[1]-2]
pp=ax.imshow(CUBE,cmap='Greys_r',vmin=totalIF[0],vmax=totalIF[1])
cbar=plt.colorbar(pp, ax=ax,  extend='both',fraction=0.046, pad=0.005,location='bottom',label='Total I/F')
ticklabs=np.linspace(totalIF[0],totalIF[1],3)
cbar.set_ticks(ticklabs)
cbar.set_ticklabels([f"{label:.3f}" for label in ticklabs])


ax.set_title(str(wavelength)+' $\mu$m')
ax.axis('off')

axs=(ax1,ax01,ax101,ax011,ax010)
axslims=(fracs1,fracs01,fracs101,fracs011,fracs010)

for c,a,l in zip(elephantcubes,axs,axslims):
    thiscube=scenariodirectoryname+str(wavelength)+'/'+elephantnameprefix+str(wavelength)+c
    #print('thiscube',thiscube)
    cube=aa.getcubefromPDS4(thiscube) 

    cube = cube[2:cube.shape[0]-2,2:cube.shape[1]-2]
    pp=a.imshow(cube/CUBE*100,cmap='Greys_r',vmin=l[0],vmax=l[1])
    cbarlab=''
    if a==ax01:
        cbarlab='Fraction of signal (%)'
    cbar=plt.colorbar(pp, ax=a,  extend='both',fraction=0.046, pad=0.005,location='bottom',label=cbarlab)
    ticklabs=np.linspace(l[0],l[1],4)
    cbar.set_ticks(ticklabs)
    if l[1]<5:
        cbar.set_ticklabels([f"{label:.1f}" for label in ticklabs])
    if l[1]<0.2:
        ticklabs=np.linspace(l[0],l[1],3)
        cbar.set_ticks(ticklabs)
        cbar.set_ticklabels([f"{label:.3f}" for label in ticklabs])
    
    else:
        cbar.set_ticklabels([f"{label:.0f}" for label in ticklabs])
        
ax1.set_title('1')
ax1.axis('off')
ax01.set_title('01')
ax01.axis('off')
ax101.set_title('101')
ax101.axis('off')
ax011.set_title('011')
ax011.axis('off')
ax010.set_title('010')
ax010.axis('off')

wavefig.savefig(figname)



wavelength=2
totalIF=[0.04,0.056]
fracs1=[20,35]
fracs01=[0.5,3]
fracs101=[0.5,3.5]
fracs011=[0.2,1]
fracs010=[0,0.015]

figname=scenariodirectoryname+'elephantpanel_'+scenariodirectoryname[0:-1]+'_'+str(wavelength)+'um.pdf'
totalcubefile=scenariodirectoryname+elephantnameprefix+str(wavelength)+'.colorCCD.xml'

print('image',totalcubefile)
#print('Final fig name = ',figname)

elephantcubes=['_____1.Jcube.colorCCD.xml','____01.Jcube.colorCCD.xml',
               '___101.Jcube.colorCCD.xml','___011.Jcube.colorCCD.xml','___010.Jcube.colorCCD.xml']

sns.set_style({'axes.grid' : False})

wavefig,[ax,ax1,ax01,ax101,ax011,ax010]=plt.subplots(1,6,figsize=(9,16),layout="constrained")

CUBE = aa.getcubefromPDS4(totalcubefile) 
CUBE=CUBE[2:CUBE.shape[0]-2,2:CUBE.shape[1]-2]
pp=ax.imshow(CUBE,cmap='Greys_r',vmin=totalIF[0],vmax=totalIF[1])
cbar=plt.colorbar(pp, ax=ax,  extend='both',fraction=0.046, pad=0.005,location='bottom',label='Total I/F')
ticklabs=np.linspace(totalIF[0],totalIF[1],3)
cbar.set_ticks(ticklabs)
cbar.set_ticklabels([f"{label:.3f}" for label in ticklabs])


ax.set_title(str(wavelength)+' $\mu$m')
ax.axis('off')

axs=(ax1,ax01,ax101,ax011,ax010)
axslims=(fracs1,fracs01,fracs101,fracs011,fracs010)

for c,a,l in zip(elephantcubes,axs,axslims):
    thiscube=scenariodirectoryname+str(wavelength)+'/'+elephantnameprefix+str(wavelength)+c
    #print('thiscube',thiscube)
    cube=aa.getcubefromPDS4(thiscube) 

    cube = cube[2:cube.shape[0]-2,2:cube.shape[1]-2]
    pp=a.imshow(cube/CUBE*100,cmap='Greys_r',vmin=l[0],vmax=l[1])
    cbarlab=''
    if a==ax01:
        cbarlab='Fraction of signal (%)'
    cbar=plt.colorbar(pp, ax=a,  extend='both',fraction=0.046, pad=0.005,location='bottom',label=cbarlab)
    ticklabs=np.linspace(l[0],l[1],4)
    cbar.set_ticks(ticklabs)
    if l[1]<5:
        cbar.set_ticklabels([f"{label:.1f}" for label in ticklabs])
    if l[1]<0.2:
        ticklabs=np.linspace(l[0],l[1],3)
        cbar.set_ticks(ticklabs)
        cbar.set_ticklabels([f"{label:.3f}" for label in ticklabs])
    
    else:
        cbar.set_ticklabels([f"{label:.0f}" for label in ticklabs])
        
ax1.set_title('1')
ax1.axis('off')
ax01.set_title('01')
ax01.axis('off')
ax101.set_title('101')
ax101.axis('off')
ax011.set_title('011')
ax011.axis('off')
ax010.set_title('010')
ax010.axis('off')

wavefig.savefig(figname)
