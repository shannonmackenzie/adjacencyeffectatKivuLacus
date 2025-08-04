import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import glob
import seaborn as sns
import sys
import adjacencyanalysis as aa
from matplotlib.ticker import MultipleLocator

pkls=glob.glob('totalimages_AC0calcs.pkl')
#print(pkls)
pkls=['landA=0_lambertianlake/totalimages_AC0calcs.pkl', 'landA=0.1_specularlake/totalimages_AC0calcs.pkl','landA=0.2_specularlake/totalimages_AC0calcs.pkl', 
      'landA=0.2_lambertianlake/totalimages_AC0calcs.pkl', 'landA=0.1_lambertianlake/totalimages_AC0calcs.pkl']


DF = pd.DataFrame()

for p in pkls:
    
    thispkl=pd.read_pickle(p)
    DF=pd.concat([DF,thispkl])
    
DF['fit_sigma_abs']=[abs(s) for s in DF.fit_sigma]
DF['fit_amplitude_abs']=[abs(a) for a in DF.fit_amplitude]
DF['aovera0']=[a/a0 if a0 !=0 else np.nan for a,a0 in zip(DF.fit_amplitude_abs,DF.landalbedo)]

DF ['aoversigma']=[a/s for a,s in zip(DF.fit_amplitude_abs,DF.fit_sigma_abs)]




sns.set_theme(context='paper', style='ticks', 
              palette='Spectral_r', font='sans-serif', font_scale=1.5, color_codes=True, 
              rc={'lines.markersize':8})

fig,[[ax5a,ax5b],[ax2a,ax2b],[ax13a,ax13b]]=plt.subplots(3,2,figsize=(10,15), layout="constrained")


lpf = 'Lambertian'
la = 0.2
Ns=[139,148,156,159,162] ##downselecting for a cleaner figure
cmap = plt.get_cmap('rainbow')
cls=cmap(np.linspace(0.1,1,len(Ns)))
cdict=dict([(key, value) for key, value in zip(Ns,cls)])

w=5
test=DF[(DF.wavelength==w)&(DF.lakephasefunc==lpf)&(DF.landalbedo==la)]
cube = aa.getcubefromPDS4(test.iloc[0].file)
vma=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].max()
vmi=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].min()

pp=ax5b.imshow(cube,vmin=vmi, vmax=vma,cmap='Greys_r')
plt.colorbar(pp, ax=ax5b,  extend='both',fraction=0.046, pad=0.04,label='I/F')
ax5b.axis('off')
ax5a.set_title(str(w)+' $\mu$m')
for i,r in test[test['edgeN'].isin(Ns)].iterrows():
    ax5b.scatter(r.edgeX,r.edgeY,marker='s',color=cdict[r.edgeN])
    ax5a.scatter(r.profile_pix,r.transect_IF,color=cdict[r.edgeN],
                label=str(round(r.fit_sigma,2))+'     '+str(round(r.fit_amplitude,3)))
    ax5a.plot(r.profile_pix,r.fit_best,color=cdict[r.edgeN])
    ax5a.legend(title='$\sigma$     $a$')
    ax5a.set(ylabel='I/F',xlabel='Distance (km)')


w=2
test=DF[(DF.wavelength==w)&(DF.lakephasefunc==lpf)&(DF.landalbedo==la)]
cube = aa.getcubefromPDS4(test.iloc[0].file)
vma=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].max()
vmi=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].min()
pp=ax2b.imshow(cube,vmin=vmi, vmax=vma,cmap='Greys_r')
plt.colorbar(pp, ax=ax2b,  extend='both',fraction=0.046, pad=0.04,label='I/F')
ax2b.axis('off')
ax2a.set_title(str(w)+' $\mu$m')
for i,r in test[test['edgeN'].isin(Ns)].iterrows():
    ax2b.scatter(r.edgeX,r.edgeY,marker='s',color=cdict[r.edgeN])
    ax2a.scatter(r.profile_pix,r.transect_IF,color=cdict[r.edgeN],
                label=str(round(r.fit_sigma,2))+'     '+str(round(r.fit_amplitude,3)))
    ax2a.plot(r.profile_pix,r.fit_best,color=cdict[r.edgeN])
    ax2a.legend(title='$\sigma$     $a$')
    ax2a.set(ylabel='I/F',xlabel='Distance (km)')

w=1.3
test=DF[(DF.wavelength==w)&(DF.lakephasefunc==lpf)&(DF.landalbedo==la)]
cube = aa.getcubefromPDS4(test.iloc[0].file)
vma=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].max()
vmi=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].min()
ax13a.set_title(str(w)+' $\mu$m')
pp=ax13b.imshow(cube,vmin=vmi, vmax=vma,cmap='Greys_r')
plt.colorbar(pp, ax=ax13b,  extend='both',fraction=0.046, pad=0.04,label='I/F')
ax13b.axis('off')
for i,r in test[test['edgeN'].isin(Ns)].iterrows():
    ax13b.scatter(r.edgeX,r.edgeY,marker='s',color=cdict[r.edgeN])
    ax13a.scatter(r.profile_pix,r.transect_IF,color=cdict[r.edgeN],
            label=str(round(r.fit_sigma,2))+'     '+str(round(r.fit_amplitude,3)))
    ax13a.plot(r.profile_pix,r.fit_best,color=cdict[r.edgeN])
    ax13a.legend(title='$\sigma$     $a$')
    ax13a.set(ylabel='I/F',xlabel='Distance (km)')
  
    
    
fig.suptitle('Land Albedo = 0.2, Lambertian Lake')