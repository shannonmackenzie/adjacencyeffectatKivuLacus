import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from osgeo import gdal

import pandas as pd
import glob
import seaborn as sns
import sys 
from lmfit.models import StepModel, ConstantModel
 
gdal.UseExceptions()

def sobeledges(cube,threshold=0.005):
    sobl = filters.sobel(cube)

    if len(np.argwhere(sobl >threshold)) == 0:
        print('WARNING: SOBL THRESHOLD TOO HIGH. No edges detected.')
    if len(np.argwhere(sobl > threshold)) >200:
        print('WARNING: SOBL THRESHOLD MAY BE TOO LOW. >200 edges detected.')
  
    return np.argwhere(sobl > threshold) 

def generateedges(threshold=0.005):
    dataset=gdal.OpenEx('lambertianlake_probablyA=0.1/5.0/Kivu_elephantA=0.1_5_____0.Jcube.colorCCD.xml')
    edgecube=dataset.GetRasterBand(1).ReadAsArray()
    dataset = None

    return sobeledges(edgecube[1:100,1:100]) 

def AC_fit(profile,transect):

    model = StepModel(form='logistic') + ConstantModel()
    params = model.make_params(center=0, sigma=1, amplitude=1., c=-0.5)
    # print('transect:\t',type(transect),transect)
    # print('params:\t',type(params),params)
    result = model.fit(transect, params, x=profile)
    return result


def getcubefromPDS4(file,band=1):
    #print('about to gdal:',file)
    dataset=gdal.OpenEx(file)
    cube=dataset.GetRasterBand(band).ReadAsArray()
    #thiscubeAC=pd.DataFrame(columns=['case','x','y','AC0_ext','AC0_int','inf_exterior','inf_interior',
    #                        'near_exterior','near_interior','inf_dist','near_dist'])
    
    dataset=None
    return cube

def exteriorx(x0,y0,lakeinteriorx,lakeinteriory,dist):
    d = np.sqrt((lakeinteriory-y0)**2 + (lakeinteriorx-x0)**2)
    if x0<lakeinteriorx:
        deltax = lakeinteriorx-x0
        return lakeinteriorx - deltax*(d+dist)/d 
    else :
        deltax = x0 - lakeinteriorx
        return deltax*(d+dist)/d+lakeinteriorx

def exteriory(x0,y0,lakeinteriorx,lakeinteriory,dist):
    d = np.sqrt((lakeinteriory-y0)**2 + (lakeinteriorx-x0)**2)
    if y0<lakeinteriory:
        deltay = lakeinteriory-y0
        return lakeinteriory - deltay*(d+dist)/d 
    else :
        deltay = y0 - lakeinteriory
        return deltay*(d+dist)/d+lakeinteriory    

def interiorx(x0,y0,lakeinteriorx,lakeinteriory,neardist):
    d = int(np.sqrt((lakeinteriory-y0)**2 + (lakeinteriorx-x0)**2))

    if x0<lakeinteriorx:
        deltax = int(lakeinteriorx - x0)
        return lakeinteriorx - deltax*(d-neardist)/d 
    else :
        deltax = int(x0 - lakeinteriorx)
        return deltax*(d-neardist)/d+lakeinteriorx

def interiory(x0,y0,lakeinteriorx,lakeinteriory,neardist):
    d = int(np.sqrt((lakeinteriory-y0)**2 + (lakeinteriorx-x0)**2))
    if y0<lakeinteriory:
        deltay = int(lakeinteriory - y0)
        return lakeinteriory - deltay*(d-neardist)/d
    else :
        deltay = int(y0 - lakeinteriory)
        return deltay*(d-neardist)/d + lakeinteriory

def getinteriorexteriorpoints(edge,cube,neardist,infdist,lakeinteriorx,lakeinteriory):
    y0=edge[0]
    x0=edge[1]

    infexteriorx = exteriorx(x0,y0,lakeinteriorx,lakeinteriory,infdist)
    infexteriory = exteriory(x0,y0,lakeinteriorx,lakeinteriory,infdist)

    nearinteriorx = interiorx(x0,y0,lakeinteriorx,lakeinteriory,neardist)
    nearinteriory = interiory(x0,y0,lakeinteriorx,lakeinteriory,neardist)

    nearexteriorx = exteriorx(x0,y0,lakeinteriorx,lakeinteriory,neardist)
    nearexteriory = exteriory(x0,y0,lakeinteriorx,lakeinteriory,neardist)

    ## check extents relative to cube. Remember that SRTC++ edge pixels are the "buckets" for
    ## all other photons, so we want to skip those
    if infexteriorx > cube.shape[1]-3: 
        infexteriorx = cube.shape[1]-3
    if infexteriory > cube.shape[0]-3: 
        infexteriory = cube.shape[0]-3
    if infexteriorx <= 1: 
        infexteriorx = 3
    if infexteriory <= 1: 
        infexteriory = 3

    return {'infexteriorx':infexteriorx,'infexteriory':infexteriory,
        'nearinteriorx':nearinteriorx,'nearinteriory':nearinteriory,
        'nearexteriorx':nearexteriorx,'nearexteriory':nearexteriory}

def calculateAEMagnitude(file,neardistarray,edges,lakeinteriorxy=[47,49],infdist=15,makeplots=True):
    ### WHOOPS. I broke the dataframes :/ March 4 2024
    ### 


    ## infdist = location of the "infitie" point, presumed to represent the "true" value in pixels
    lakeinteriorx,lakeinteriory = lakeinteriorxy 
        
    count = 0
    casestr=file[file.rfind('/')+1:file.rfind('.colorCCD')]
    if '.Jcube' in file:
        casestr=file[file.rfind('/')+1:file.rfind('.Jcube')]
    
    # print(file)
    cube = getcubefromPDS4(file)
    gcubex = getcubefromPDS4(file[:-13]+'_geo'+file[-13:],band=3) #right now, SRTC++ does square pixels, so x=y
    #gcubey=getcubefromPDS4(file[:-13]+'_geo'+file[-13:],band=4)
    #print('read in cube and geocube')
    wavelength=file[file.find('/')+1:file.rfind('/')]
    print(wavelength)
    
    # thiscubeAC=pd.DataFrame(columns=['case','x','y','AC0_ext','AC0_int','inf_exterior','inf_interior',
    #                         'near_exterior','near_interior','inf_dist','near_dist'])
    thiscubeAC=pd.DataFrame(columns=['case','wavelength','edgeY','edgeX','edgeN','transect_IF','profile_pix','fit_center',
                    'fit_amplitude','fit_sigma','3sigma','fit_best','AC0_ext','AC0_int',
                    'AC0_ext','AC0_int','inf_dist','near_dist','kmperpixel'])

    if cube.mean() == cube.max() == cube.min():
        print(casestr,'has no signal.')
    
    else:
        if makeplots:
            fig,[ax1,ax2]=plt.subplots(1,2,figsize=(10,5), layout="constrained")
            pp=ax2.imshow(cube,vmin=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].min(), vmax=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].max(),cmap='Greys_r')
            plt.colorbar(pp, ax=ax2,  extend='both',fraction=0.046, pad=0.04)
            ax1.set_title(casestr)
        
        for i,p in enumerate(edges):
            # calculate AC0
            
            y0=p[0]
            x0=p[1] 
            
            # generate profile xys inbetween inf interior and inf exterior
            profilexs = [exteriorx(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]
            profilexs.reverse()
            profileys = [exteriory(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]
            profileys.reverse()
            profilexs += [interiorx(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]
            profileys += [interiory(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]

            transect = [cube[int(y),int(x)] for y,x in zip(profileys,profilexs)]
            kmperpixel = gcubex[int(lakeinteriory),int(lakeinteriorx)] 
            profdists = np.arange(-infdist,infdist)
            fitresult=AC_fit(profdists,transect)
            for n in neardistarray:
                neardist=n
                
                ## get points for AC0 
                thesepoints=getinteriorexteriorpoints(p,cube,neardist,infdist,lakeinteriorx,lakeinteriory)
            
                ## probably should build in that limits check somewhere....
                IA_inf_plus=cube[int(thesepoints['infexteriory']),int(thesepoints['infexteriorx'])]
                IA_inf_minus=cube[int(lakeinteriory),int(lakeinteriorx)]
                IA_0_plus=cube[int(thesepoints['nearexteriory']),int(thesepoints['nearexteriorx'])]
                IA_0_minus=cube[int(thesepoints['nearinteriory']),int(thesepoints['nearinteriorx'])]

                ac_0=AC_0(IA_inf_plus,IA_inf_minus,IA_0_plus,IA_0_minus)
                # thisresult=pd.DataFrame({'case':casestr,'x':p[1],'y':p[0],
                #                 'AC0_ext':ac_0[0],
                #                 'AC0_int':ac_0[1],
                #                 'AC_fit_sigma':fitresult.best_values['sigma'],
                #                 'AC_fit_constant':fitresult.best_values['c'],
                #                 'AC_fit_center':fitresult.best_values['center'],
                #                 'AC_fit_amplitude':fitresult.best_values['amplitude'],
                #                 'AC_fit_rsquared':fitresult.rsquared,
                #                 'inf_exterior':IA_inf_plus,
                #                 'inf_interior':IA_inf_minus,
                #                 'near_exterior':IA_0_plus,
                #                 'near_interior':IA_0_minus,
                #                 'inf_dist':infdist,'near_dist':neardist,'kmperpixel':kmperpixel},index=[count])
                thisresult=pd.DataFrame({'case':casestr,'wavelength':float(wavelength),'edgeY':p[0],'edgeX':p[1],'edgeN':i,
                    'transect_IF':[transect],'profile_pix':[profdists],'fit_center':fitresult.best_values['c'],
                    'fit_amplitude':fitresult.best_values['amplitude'],'fit_sigma':fitresult.best_values['sigma'],
                    '3sigma':[fitresult.eval_uncertainty(sigma=3)],'fit_best':[fitresult.best_fit],
                    'AC0_ext':ac_0[0],'AC0_int':ac_0[1],'inf_dist':infdist,'near_dist':neardist,'kmperpixel':kmperpixel},index=[count])
                
                thiscubeAC=pd.concat([thiscubeAC.reset_index(),thisresult.reset_index()],ignore_index=True)
                count+=1
                
                # if (i % 12 == 0) & makeplots & (neardist ==3):
                #     ax2.plot(x0,y0,marker='s')
                #     ax2.plot(profilexs,profileys,marker='s',markersize=1,
                #          label=round(thisresult['AC0_int'].iloc[0],4))
                    
                #     ax1.plot(profdists,transect)

                # if makeplots & (np.isnan(thiscubeAC.AC0_int.max()) is False) & (neardist ==3):
                #     plt.savefig(f[:-12]+'.transectprofiles.png')
                    
                #     fig, [ax1,ax2]=plt.subplots(1,2,figsize=(16,9))
                #     a1=ax1.imshow(cube,cmap='Greys_r',vmax=cube[1:100,1:100].max(),vmin=cube[1:100,1:100].min())
                #     ax2.imshow(cube,cmap='Greys_r',vmax=cube[1:100,1:100].max(),vmin=cube[1:100,1:100].min())
                #     a2=ax2.scatter(thiscubeAC.x,thiscubeAC.y, marker='s', s=5,c=thiscubeAC.AC0_int,cmap=plt.cm.coolwarm,vmax=thiscubeAC.AC0_int.max(),vmin=thiscubeAC.AC0_int.min())
                #     fig.colorbar(a1, ax=ax1,fraction=0.04,format=lambda x, _: f"{x:.3}")
                #     fig.colorbar(a2, ax=ax2,fraction=0.04,format=lambda x, _: f"{x:.0}")
                #     ax1.set_title(casestr)
                #     #ax2.set_title('AC$_0$ \n (max='+str(thiscubeAC.AC0.max())+', min='+str(thiscubeAC.AC0.min())+', $\sigma$='+str(thiscubeAC.AC0.std()))
                #     if '.Jcube' in file:
                #         plt.savefig(file[:-12]+'.AC0calcs.png')
                #     else:
                #         plt.savefig(file[:-12]+'.AC0calcs.png')

                #     plt.close('all')


    return thiscubeAC



def calculateAEMagnitude_fitsonly(file,edges,lakeinteriorxy=[47,49],infdist=15,makeplots=True):
    
    ## infdist = location of the "infitie" point, presumed to represent the "true" value in pixels
    lakeinteriorx,lakeinteriory = lakeinteriorxy 
        
    count = 0
    casestr=file[file.rfind('/')+1:file.rfind('.colorCCD')]
    if '.Jcube' in file:
        casestr=file[file.rfind('/')+1:file.rfind('.Jcube')]
    
    # print(file)
    cube = getcubefromPDS4(file)
    gcubex = getcubefromPDS4(file[:-13]+'_geo'+file[-13:],band=3) #right now, SRTC++ does square pixels, so x=y
    #gcubey=getcubefromPDS4(file[:-13]+'_geo'+file[-13:],band=4)
    #print('read in cube and geocube')
    wavelength=file[file.find('/')+1:file.rfind('/')]
    print(wavelength)
    
    thiscubeAC=pd.DataFrame(columns=['case','wavelength','edgeY','edgeX','edgeN','transect_IF','profile_pix','fit_center',
                    'fit_amplitude','fit_sigma','3sigma','fit_best','kmperpixel'])



    if cube.mean() == cube.max() == cube.min():
        print(casestr,'has no signal.')
    
    else:
        if makeplots:
            fig,[ax1,ax2]=plt.subplots(1,2,figsize=(10,5), layout="constrained")
            pp=ax2.imshow(cube,vmin=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].min(), vmax=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].max(),cmap='Greys_r')
            plt.colorbar(pp, ax=ax2,  extend='both',fraction=0.046, pad=0.04)
            ax1.set_title(casestr)
        
        for i,p in enumerate(edges):
            # calculate AC0
            
            y0=p[0]
            x0=p[1] 
            
            # generate profile xys inbetween inf interior and inf exterior
            profilexs = [exteriorx(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]
            profilexs.reverse()
            profileys = [exteriory(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]
            profileys.reverse()
            profilexs += [interiorx(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]
            profileys += [interiory(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]

            transect = [cube[int(y),int(x)] for y,x in zip(profileys,profilexs)]
            kmperpixel = gcubex[int(lakeinteriory),int(lakeinteriorx)] 
            profdists = np.arange(-infdist,infdist)
            fitresult=AC_fit(profdists,transect)
            

            thisresult=pd.DataFrame({'case':casestr,'wavelength':float(wavelength),'edgeY':p[0],'edgeX':p[1],'edgeN':i,
                'transect_IF':[transect],'profile_pix':[profdists],'fit_center':fitresult.best_values['c'],
                'fit_amplitude':fitresult.best_values['amplitude'],'fit_sigma':fitresult.best_values['sigma'],
                '3sigma':[fitresult.eval_uncertainty(sigma=3)],'fit_best':[fitresult.best_fit],'kmperpixel':kmperpixel},index=[count])


            thiscubeAC=pd.concat([thiscubeAC,thisresult],ignore_index=True)
            if (i % 12 == 0) & makeplots :
                ax2.plot(x0,y0,marker='s')
                ax2.plot(profilexs,profileys,marker='s',markersize=1,
                     label=round(thisresult['AC0_int'].iloc[0],4))
                
                ax1.plot(profdists,transect)

            if makeplots & (np.isnan(thiscubeAC.AC0_int.max()) is False) & (neardist ==3):
                plt.savefig(f[:-12]+'.transectprofiles.png')
                
                fig, [ax1,ax2]=plt.subplots(1,2,figsize=(16,9))
                a1=ax1.imshow(cube,cmap='Greys_r',vmax=cube[1:100,1:100].max(),vmin=cube[1:100,1:100].min())
                ax2.imshow(cube,cmap='Greys_r',vmax=cube[1:100,1:100].max(),vmin=cube[1:100,1:100].min())
                a2=ax2.scatter(thiscubeAC.x,thiscubeAC.y, marker='s', s=5,c=thiscubeAC.AC0_int,cmap=plt.cm.coolwarm,vmax=thiscubeAC.AC0_int.max(),vmin=thiscubeAC.AC0_int.min())
                fig.colorbar(a1, ax=ax1,fraction=0.04,format=lambda x, _: f"{x:.3}")
                fig.colorbar(a2, ax=ax2,fraction=0.04,format=lambda x, _: f"{x:.0}")
                ax1.set_title(casestr)
                #ax2.set_title('AC$_0$ \n (max='+str(thiscubeAC.AC0.max())+', min='+str(thiscubeAC.AC0.min())+', $\sigma$='+str(thiscubeAC.AC0.std()))
                if '.Jcube' in file:
                    plt.savefig(file[:-12]+'.AC0calcs.png')
                else:
                    plt.savefig(file[:-12]+'.AC0calcs.png')

                plt.close('all')
                    
    return thiscubeAC


## version for 2024 file naming scheme. Be wary!
def justransect(file,edges,scenariostr,lakeinteriorxy=[47,49],infdist=15):
    
    infdist= 15 # location of the "infitie" point, presumed to represent the "true" value in pixels
    lakeinteriorx,lakeinteriory = lakeinteriorxy 
        
    count = 0
    casestr=file[file.rfind('/')+1:file.rfind('.colorCCD')]
    if '.Jcube' in file:
        casestr=file[file.rfind('/')+1:file.rfind('.Jcube')]
    
    print(file)
    cube = getcubefromPDS4(file)
    cube = getcubefromPDS4(file)
    gcubex = getcubefromPDS4(file[:-13]+'_geo'+file[-13:],band=3) #right now, SRTC++ does square pixels, so x=y
    
    wavelength=file[file.rfind('_')+1:file.rfind('/')]
    if wavelength == '':
        wavelength=file[file.rfind('_')+1:file.rfind('.co')]
    if '.Jcube' in file:
        wavelength=0 # do it in post; too hard to grab numbers from file naming convention
    
    print(wavelength)
    
    thiscubeAC=pd.DataFrame(columns=['case','wavelength','edgeY','edgeX','edgeN','transect_IF','profile_pix','fit_center',
                    'fit_amplitude','fit_sigma','fit_3sigma','fit_redchisqd','fit_best','kmperpixel'])


    if cube.mean() == cube.max() == cube.min():
        print(casestr,'has no signal.')
    
    else:
        fig,[ax1,ax2]=plt.subplots(1,2,figsize=(10,5), layout="constrained")
        pp=ax2.imshow(cube,vmin=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].min(),vmax=cube[2:cube.shape[0]-2,2:cube.shape[1]-2].max(),cmap='Greys_r')
        plt.colorbar(pp, ax=ax2,  extend='both',fraction=0.046, pad=0.04)
        ax1.set_title(casestr)
        
        fig2, [ax21,ax22,ax23]=plt.subplots(3,1,figsize=(9,16),layout='tight')
        a21=ax21.imshow(cube,cmap='Greys_r',vmax=cube[1:100,1:100].max(),vmin=cube[1:100,1:100].min())
        ax22.imshow(cube,cmap='Greys_r',vmax=cube[1:100,1:100].max(),vmin=cube[1:100,1:100].min())
        ax23.imshow(cube,cmap='Greys_r',vmax=cube[1:100,1:100].max(),vmin=cube[1:100,1:100].min())
        amplitudes=[]
        sigmas=[]
        for i,p in enumerate(edges):

            y0=p[0]
            x0=p[1] 

            # generate profile xys inbetween inf interior and inf exterior
            profilexs = [exteriorx(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]
            profilexs.reverse()
            profileys = [exteriory(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]
            profileys.reverse()
            profilexs += [interiorx(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]
            profileys += [interiory(x0,y0,lakeinteriorx,lakeinteriory,d) for d in np.arange(0,infdist)]

            transect = [cube[int(y),int(x)] for y,x in zip(profileys,profilexs)]
            profdists = np.arange(-infdist,infdist)
            fitresult=AC_fit(profdists,transect)
            
            ##clean up outliers; 2024-03-20
            if abs(fitresult.best_values['sigma']) < 15: 
                sigmas.append(abs(fitresult.best_values['sigma']))  
            else:
                sigmas.append(-1)
            
            if abs(fitresult.best_values['amplitude']) < 0.2:
                amplitudes.append(abs(fitresult.best_values['amplitude']))
            else:
                amplitudes.append(-1)
                
            kmperpixel = gcubex[int(lakeinteriory),int(lakeinteriorx)] 

            thisresult=pd.DataFrame({'case':casestr,'wavelength':float(wavelength),'edgeY':p[0],'edgeX':p[1],'edgeN':i,
                    'transect_IF':[transect],'profile_pix':[profdists],'fit_center':fitresult.best_values['c'],
                    'fit_amplitude':fitresult.best_values['amplitude'],'fit_sigma':fitresult.best_values['sigma'],
                    'fit_3sigma':[fitresult.eval_uncertainty(sigma=3)],'fit_redchisqd':fitresult.redchi,
                    'fit_best':[fitresult.best_fit],'kmperpixel':kmperpixel},index=[count])
                
            thiscubeAC=pd.concat([thiscubeAC,thisresult],ignore_index=True)
            count+=1    
            if (i % 12 == 0)  :
                ax2.plot(x0,y0,marker='s')
                ax2.plot(profilexs,profileys,marker='s',markersize=1)
                ax1.plot(profdists,transect)
        
        
        ax21.set_title(str(wavelength)+r' $\mu$m')


        a22=ax22.scatter([p[1] for p in edges],[p[0] for p in edges], marker='s', s=5,c=sigmas,cmap=plt.cm.coolwarm,vmin=0,vmax=np.max(sigmas))
        ax22.set_title(r'$\sigma$ (km)')
        
        a23=ax23.scatter([p[1] for p in edges],[p[0] for p in edges], marker='s', s=5,c=amplitudes,cmap=plt.cm.coolwarm,vmin=0,vmax=np.max(amplitudes))
        ax23.set_title(r'$\Delta$I/F ')
            
        fig2.colorbar(a21, ax=ax21,fraction=0.04,format=lambda x, _: f"{x:.3}")
        fig2.colorbar(a22, ax=ax22,fraction=0.04,format=lambda x, _: f"{x:.4}")
        fig2.colorbar(a23, ax=ax23,fraction=0.04,format=lambda x, _: f"{x:.4}")
        
        if '.Jcube' in file:
            fig2.suptitle(casestr)
            fig.savefig(file[:-18]+'transectprofiles.png')
            fig2.savefig(file[:-18]+'fitcalcs.png')
        else:
            fig2.suptitle(scenariostr)
            fig.savefig(file[:-12]+'transectprofiles.png')
            fig2.savefig(file[:-12]+'fitcalcs.png')
        
        plt.close('all')
    return thiscubeAC
                    
