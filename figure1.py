import pyvims as pv
from pyvims import titan,plot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.pyplot import cm
from skimage import filters
from lmfit.models import StepModel, ConstantModel
import seaborn as sns



#### -- Figure 1 --
kivu=pv.VIMS('C1721854264_1_ir') #T85
context=pv.VIMS('C1753530947_1_ir') # T93
edges=np.load('T85edges.npy')

#### 1a
fig,ax1=plt.subplots(1,1,figsize=(10,10))
pv.plot.plot_polar(context,('339:351', '138:141', '121:122'),ticks=True,grid='white',lat_min=86,ax=ax1)
pv.plot.plot_polar(kivu,('339:351', '138:141', '121:122'),ticks=True,grid='white',lat_min=86,ax=ax1)
fig.savefig('maps.pdf',facecolor=fig.get_facecolor(), edgecolor='none')

#### 1b angles
fig,[ax1,ax2,ax3] = plt.subplots(1,3)
l1=ax1.imshow(kivu.phase)
plt.colorbar(l1, ax=ax1,  extend='both',fraction=0.046, pad=0.04)
ax1.set_title('Phase')

l2=ax2.imshow(kivu.inc)
plt.colorbar(l2, ax=ax2,  extend='both',fraction=0.046, pad=0.04)
ax2.set_title('Incidence')

l3=ax3.imshow(kivu.eme)
plt.colorbar(l3, ax=ax3,  extend='both',fraction=0.046, pad=0.04)
ax3.set_title('Emission')
plt.tight_layout()



#### 1b wavelengths 0.93, 1.08, 1.6, 2.7, 2.8um

sns.set_theme(context='paper', style='ticks', palette='Spectral_r', font='sans-serif', font_scale=1.5, color_codes=True, rc={'lines.markersize':8})

fig,[ax1,ax2,ax3,ax4,ax5]=plt.subplots(1,5,figsize=(6,12))

ax1.imshow(kivu@(0.93),cmap='Greys_r')
ax1.set_title('0.93 $\mu$m')
ax1.axis('off')

ax2.imshow(kivu@(1.08),cmap='Greys_r')
ax2.set_title('1.08 $\mu$m')
ax2.axis('off')

ax3.imshow(kivu@(1.6),cmap='Greys_r')
ax3.set_title('1.6 $\mu$m')
ax3.axis('off')

ax4.imshow(kivu@(2.7),cmap='Greys_r')
ax4.set_title('2.7 $\mu$m')
ax4.axis('off')

ax5.imshow(kivu@(2.8),cmap='Greys_r')
ax5.set_title('2.8 $\mu$m')
ax5.axis('off')




#### 1b wavelengths 1.28, 2, 5 um

cube5um=kivu[336:352]
cube2um=kivu[164:169]
cube13um=kivu[121:122]

cmap = plt.get_cmap('viridis')
cls=cmap(np.linspace(0.1,1,len(edge)))
cdict=dict([(key, value) for key, value in enumerate(cls)])


sns.set_theme(context='paper', style='ticks', palette='Spectral_r', font='sans-serif', font_scale=1.5, color_codes=True, rc={'lines.markersize':8})

fig,[ax1,ax2,ax3]=plt.subplots(1,3,figsize=(6,9),)

ax1.imshow(cube13um,cmap='Greys_r')
for i,p in enumerate(edge):
    if i%4 ==0:
        ax1.plot(p[1],p[0],c=cdict[i],marker='s',markersize=3,alpha=0.75)
ax1.set_title('1.28 $\mu$m')
ax1.axis('off')

ax2.imshow(cube2um,cmap='Greys_r')
for i,p in enumerate(edge):
    if i%4 ==0:
        ax2.plot(p[1],p[0],c=cdict[i],marker='s',markersize=3,alpha=0.75)
ax2.set_title('2 $\mu$m')
ax2.axis('off')


ax3.imshow(cube5um,cmap='Greys_r')
for i,p in enumerate(edge):
    if i%4 ==0:
        ax3.plot(p[1],p[0],c=cdict[i],marker='s',markersize=3,alpha=0.75)
ax3.set_title('5 $\mu$m')
ax3.axis('off')

fig.tight_layout()



#### 1c transects
plt.style.use('seaborn-white')


kbDF=pd.DataFrame(columns=['wavelength','row','edgex','edgey','transect_IF','profile_pix',
                             'fit_best','AC_fit_sigma','AC_fit_constant','AC_fit_center','AC_fit_amplitude','AC_fit_rsquared'])

fig,ax=plt.subplots(2,4,figsize=(20,10))
ax=ax.flatten()
cmap = plt.get_cmap('viridis')
cls=cmap(np.linspace(0,1,32))

for i,w in enumerate(waveranges.keys()):
    
    for l,e in zip(np.arange(0,kivub.shape[2]),edgeb):
        #print('Analyzing Row ',l,' with edges',e)
        pixdist2boundary=[x-e[1] for x in np.arange(1,kivub.shape[1]+1)]
        profile=kivub.res[l,:]*pixdist2boundary
        transect=kivub[waveranges[w][0]:waveranges[w][1]][l,0:32]
        
        pixs=[kivub@(l+1,x+1) for x in np.arange(0,32)]
        incs=[p.inc for p in pixs]
        emes=[p.eme for p in pixs]
        phases=[p.phase for p in pixs]
        
        edgepix=kivub@(e[1]+1,e[0]+1)
        edgei=edgepix.inc
        edgee=edgepix.eme
        edgep=edgepix.phase
    
        thisresult=aa.AC_fit(profile,transect)
        
        lakeinf=transect[0]
        lakeboundary=transect[e[1]-1]
        landboundary=transect[e[1]+1]
        landinf=transect[-1]
        
        thisdf=pd.DataFrame({'wavelength':w,'row':l,'edgex':e[1],'edgey':e[0],'edgei':edgei,'edgee':edgee,'edgep':edgep,
                             'rowincmean':np.mean(incs),'rowincstd':np.std(incs),
                             'rowememean':np.mean(emes),'rowemestd':np.std(emes),
                             'rowphasemean':np.mean(phases),'rowphasestd':np.std(phases),
                             'transect_IF':[transect],'profile_pix':[profile],
                             'fit_best':[thisresult.best_fit],
                             'AC0_lake':aa.AC_0(lakeinf,landinf,lakeboundary,landboundary)[0],
                             'AC_fit_sigma':thisresult.best_values['sigma'],
                             'AC_fit_constant':thisresult.best_values['c'],
                             'AC_fit_center':thisresult.best_values['center'],
                             'AC_fit_amplitude':thisresult.best_values['amplitude'],
                             'AC_fit_rsquared':thisresult.rsquared},index=np.arange(0,kivub.shape[2]))
        kbDF=pd.concat([kbDF,thisdf])
        
        
        
        if l%4 ==0:
            
            #dely = thisresult.eval_uncertainty(sigma=3)
            ax[i].scatter(profile,transect,color=cls[l])
            ax[i].plot(profile, thisresult.best_fit,c=cls[l],label=l)
            #ax[i].fill_between(profile, thisresult.best_fit-dely, thisresult.best_fit+dely, color=cls[l],alpha=0.5)
    
    ax[i].set_title(str(w)+' $\mu$m')

    ax[i].set(xlabel='Distance from boundary (km)', ylabel='I/F')

kbDF['fit_sigma_abs']=[abs(f) for f in kbDF.AC_fit_sigma]
kbDF['fit_amplitude_abs']=[abs(f) for f in kbDF.AC_fit_amplitude]

      
plt.tight_layout()

### fig 1d

#generated in 'VIMS Transects for paper' notebook
obs=pd.read_csv('observedKivusigma.csv')
mod=pd.read_csv('modeledKivusigma_swonly.csv')
data=pd.concat([obs,mod])
data=data.drop(columns=['Unnamed: 0'])


fig,ax1=plt.subplots(1,1,figsize=(6,4))

ax1.errorbar(data[data.scenario=='C1721854264_1'].wavelength,
             data[data.scenario=='C1721854264_1'].sig_avg,
             yerr=data[data.scenario=='C1721854264_1'].sig_sd,
             marker='o',linestyle='-',label='C1721854264_1',c='k')

ax1.errorbar(data[data.scenario=='Land Albedo = 0.1; Specular Lake'].wavelength,
             data[data.scenario=='Land Albedo = 0.1; Specular Lake'].sig_avg,color=scenariocolordict['Land Albedo = 0.1; Specular Lake'],
             yerr=data[data.scenario=='Land Albedo = 0.1; Specular Lake'].sig_sd,marker='o',
             linestyle='-',label='Land Albedo = 0.1; Specular Lake')

ax1.errorbar(data[data.scenario=='Land Albedo = 0.2; Specular Lake'].wavelength,
             data[data.scenario=='Land Albedo = 0.2; Specular Lake'].sig_avg,color=scenariocolordict['Land Albedo = 0.2; Specular Lake'],
             yerr=data[data.scenario=='Land Albedo = 0.2; Specular Lake'].sig_sd,marker='o',
             linestyle='-',label='Land Albedo = 0.2; Specular Lake')

ax1.errorbar(data[data.scenario=='Land Albedo = 0.1; Lambertian Lake'].wavelength,
             data[data.scenario=='Land Albedo = 0.1; Lambertian Lake'].sig_avg,color=scenariocolordict['Land Albedo = 0.1; Lambertian Lake'],
             yerr=data[data.scenario=='Land Albedo = 0.1; Lambertian Lake'].sig_sd,marker='o',
             linestyle='-',label='Land Albedo = 0.1; Lambertian Lake')

ax1.errorbar(data[data.scenario=='Land Albedo = 0.2; Lambertian Lake'].wavelength,
             data[data.scenario=='Land Albedo = 0.2; Lambertian Lake'].sig_avg,color=scenariocolordict['Land Albedo = 0.2; Lambertian Lake'],
             yerr=data[data.scenario=='Land Albedo = 0.2; Lambertian Lake'].sig_sd,marker='o',
             linestyle='-',label='Land Albedo = 0.2; Lambertian Lake')

ax1.hlines(0.8290924780979793,0.8,5.2,linestyle='--',color='k',alpha=0.5,linewidth=2) #max res of kivu cube


ax1.set_xlim(0.85,5.1)
ax1.set_ylim(0,5)
ax1.set_xlabel('Wavelength ($\mu$m)')
ax1.set_ylabel('$\sigma$ (km)')

plt.tight_layout()


## fig 1e

fig,ax1=plt.subplots(1,1,figsize=(6,4))

ax1.errorbar(data[data.scenario=='C1721854264_1'].wavelength,
             data[data.scenario=='C1721854264_1'].amp_avg,
             yerr=data[data.scenario=='C1721854264_1'].amp_sd,
             marker='o',linestyle='-',label='C1721854264_1',c='k')

ax1.errorbar(data[data.scenario=='Land Albedo = 0.1; Specular Lake'].wavelength,
             data[data.scenario=='Land Albedo = 0.1; Specular Lake'].amp_avg,color=scenariocolordict['Land Albedo = 0.1; Specular Lake'],
             yerr=data[data.scenario=='Land Albedo = 0.1; Specular Lake'].amp_sd,marker='o',
             linestyle='-',label='Land Albedo = 0.1; Specular Lake')

ax1.errorbar(data[data.scenario=='Land Albedo = 0.2; Specular Lake'].wavelength,
             data[data.scenario=='Land Albedo = 0.2; Specular Lake'].amp_avg,color=scenariocolordict['Land Albedo = 0.2; Specular Lake'],
             yerr=data[data.scenario=='Land Albedo = 0.2; Specular Lake'].amp_sd,marker='o',
             linestyle='-',label='Land Albedo = 0.2; Specular Lake')

ax1.errorbar(data[data.scenario=='Land Albedo = 0.1; Lambertian Lake'].wavelength,
             data[data.scenario=='Land Albedo = 0.1; Lambertian Lake'].amp_avg,color=scenariocolordict['Land Albedo = 0.1; Lambertian Lake'],
             yerr=data[data.scenario=='Land Albedo = 0.1; Lambertian Lake'].amp_sd,marker='o',
             linestyle='-',label='Land Albedo = 0.1; Lambertian Lake')

ax1.errorbar(data[data.scenario=='Land Albedo = 0.2; Lambertian Lake'].wavelength,
             data[data.scenario=='Land Albedo = 0.2; Lambertian Lake'].amp_avg,color=scenariocolordict['Land Albedo = 0.2; Lambertian Lake'],
             yerr=data[data.scenario=='Land Albedo = 0.2; Lambertian Lake'].amp_sd,marker='o',
             linestyle='-',label='Land Albedo = 0.2; Lambertian Lake')

ax1.set_xlim(0.85,5.1)
ax1.set_ylim(0,0.03)
ax1.set_xlabel('Wavelength ($\mu$m)')
ax1.set_ylabel('$a$ ($\Delta$I/F)')

plt.locator_params(axis='x', nbins=6)

plt.tight_layout()









