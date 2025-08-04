import pyvims as pv
from pyvims import titan,plot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.pyplot import cm
from skimage import filters
from lmfit.models import StepModel, ConstantModel
import seaborn as sns


sns.set_style({'axes.grid' : False})
w = 5
scenario='Kivu_elephantA=0.1_'
directory=wavedirs[0]
figname='landA=0.1_'+str(w)+'um_elephants.pdf'
totalcubefile='landA=0.1/'+scenario+str(w)+'.colorCCD.xml'


wavefig,[[ax,ax0,ax1],[ax01,ax010,ax011],[ax10,ax100,ax101]]=plt.subplots(3,3,figsize=(9,9),layout="constrained")
CUBE = aa.getcubefromPDS4(totalcubefile) 
CUBE=CUBE[2:CUBE.shape[0]-2,2:CUBE.shape[1]-2]
pp=ax.imshow(CUBE,cmap='Greys_r')
plt.colorbar(pp, ax=ax,  extend='both',fraction=0.046, pad=0.005)
ax.set_title(str(w)+' $\mu$m Total')
ax.axis('off')

axs=(ax0,ax1,ax01,ax010,ax011,ax10,ax100,ax101)
for c,a in zip(cubes,axs):
    cube=aa.getcubefromPDS4(directory+scenario+str(w)+c) 
    cube = cube[2:cube.shape[0]-2,2:cube.shape[1]-2]
    pp=a.imshow(cube/CUBE,cmap='Greys_r')
    plt.colorbar(pp, ax=a,  extend='both',fraction=0.046, pad=0.005)

    
ax0.set_title('0')
ax0.axis('off')
ax1.set_title('1')
ax1.axis('off')
ax01.set_title('01')
ax01.axis('off')
ax010.set_title('010')
ax010.axis('off')
ax011.set_title('011')
ax011.axis('off')
ax10.set_title('10')
ax10.axis('off')
ax100.set_title('100')
ax100.axis('off')
ax101.set_title('101')
ax101.axis('off')
wavefig.suptitle(suptitlestr+'\n'+str(w)+' $\mu$m')
wavefig.savefig(figname)


sns.set_style({'axes.grid' : False})
w = 2
scenario='Kivu_elephantA=0.1_'
directory=wavedirs[3]
suptitlestr='Land Albedo = 0.1'
totalcubefile='landA=0.1/Kivu_elephantA=0.1_2.colorCCD.xml'
figname='landA=0.1_'+str(w)+'um_elephants.pdf'


wavefig,[[ax,ax0,ax1],[ax01,ax010,ax011],[ax10,ax100,ax101]]=plt.subplots(3,3,figsize=(9,9),layout="constrained")
CUBE = aa.getcubefromPDS4(totalcubefile) 
CUBE=CUBE[2:CUBE.shape[0]-2,2:CUBE.shape[1]-2]
pp=ax.imshow(CUBE,cmap='Greys_r')
plt.colorbar(pp, ax=ax,  extend='both',fraction=0.046, pad=0.005)
ax.set_title(str(w)+' $\mu$m Total')
ax.axis('off')

axs=(ax0,ax1,ax01,ax010,ax011,ax10,ax100,ax101)
for c,a in zip(cubes,axs):
    cube=aa.getcubefromPDS4(directory+scenario+str(w)+c) 
    cube = cube[2:cube.shape[0]-2,2:cube.shape[1]-2]
    pp=a.imshow(cube/CUBE,cmap='Greys_r')
    plt.colorbar(pp, ax=a,  extend='both',fraction=0.046, pad=0.005)

    
ax0.set_title('0')
ax0.axis('off')
ax1.set_title('1')
ax1.axis('off')
ax01.set_title('01')
ax01.axis('off')
ax010.set_title('010')
ax010.axis('off')
ax011.set_title('011')
ax011.axis('off')
ax10.set_title('10')
ax10.axis('off')
ax100.set_title('100')
ax100.axis('off')
ax101.set_title('101')
ax101.axis('off')

wavefig.suptitle(suptitlestr+'\n'+str(w)+' $\mu$m')
wavefig.savefig(figname)