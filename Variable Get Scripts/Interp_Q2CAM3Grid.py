#!/usr/bin/env python

## Import Packages
import numpy as np
import scipy.stats as stats
import os
import pickle
import scipy.io as scio
import pickle
import Area_Avg
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import griddata
from netCDF4 import Dataset
import warnings
from matplotlib import pyplot as plt
warnings.filterwarnings('ignore')
np.set_printoptions(suppress=True)

#Load in Dimensions
LLL = Dataset('/export/data/waccmdata/b40.rcp4_5.2deg.wcm.001/atm/b40.rcp4_5.2deg.wcm.001.2005.01.nc')
WACCM4_Lat = np.squeeze(LLL.variables['lat'])
WACCM4_Lon = np.squeeze(LLL.variables['lon'])
LLL.close()

LLL3 = Dataset('CAM3 Kernels/CAM3_AirTemperature_Kernel.nc')
CAM3_Lat = np.squeeze(LLL3.variables['lat'])
CAM3_Lon = np.squeeze(LLL3.variables['lon'])

#lists for iterating
Months = np.arange(0,12)
FP_Ensemble_Members = np.arange(0,5)
Ensemble_number = ['1','2','3','4','5']

CMIP_plevs = np.squeeze(LLL3.variables['lev'])

WACCM4_plevs_RCP45 = np.zeros([5,12,47,96,144])
WACCM4_plevs_wAVD = np.zeros([5,12,47,96,144])
for i in range(5):
    WACCM4_plevs_RCP45[i,:,:,:,:] = np.squeeze(Dataset(\
    'Midpoint Pressure Levels/Future Projection/plevs_mean_RCP45E'\
    +Ensemble_number[i]+'_WACCM4.nc').variables['pmid'])/100

    WACCM4_plevs_wAVD[i,:,:,:,:] = np.squeeze(Dataset(\
    'Midpoint Pressure Levels/Future Projection/plevs_mean_wAVDE'\
    +Ensemble_number[i]+'_WACCM4.nc').variables['pmid'])/100

#make meshgrids of all horizontal resolutions for regridding
cx,cy = np.meshgrid(CAM3_Lat,CAM3_Lon)
wx,wy = np.meshgrid(WACCM4_Lat,WACCM4_Lon)

RCP45_dQ = pickle.load(open(\
    "Future Projection CC Responses/RCP45_dlnQ.pickle","rb"))
wAVD_dQ = pickle.load(open(\
    "Future Projection CC Responses/wAVD_dlnQ.pickle","rb"))

WACCM4_plevs_RCP45_camHgrid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
WACCM4_plevs_wAVD_camHgrid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            WACCM4_plevs_RCP45_camHgrid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),WACCM4_plevs_RCP45[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')

            WACCM4_plevs_wAVD_camHgrid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),WACCM4_plevs_wAVD[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')


#water vapour
RCP45_dQ_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            RCP45_dQ_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),RCP45_dQ[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')

RCP45_dQ_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    RCP45_dQ_Tprof=np.squeeze(RCP45_dQ_CAM3Grid[e,month,:,x,y])
                    RCP45_dQ_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_RCP45_camHgrid[e,month,:,x,y]),\
                    RCP45_dQ_Tprof,CMIP_plevs,method="nearest")


wAVD_dQ_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            wAVD_dQ_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),wAVD_dQ[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')

wAVD_dQ_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    wAVD_dQ_Tprof=np.squeeze(wAVD_dQ_CAM3Grid[e,month,:,x,y])
                    wAVD_dQ_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_wAVD_camHgrid[e,month,:,x,y]),\
                    wAVD_dQ_Tprof,CMIP_plevs,method="nearest")

#save all interpolated arrays to new files

RCP45_dQ_CAM3Grid_Full_file = open("RCP45_dlnQ_CAM3Grid.pickle","wb")
pickle.dump(RCP45_dQ_CAM3Grid_Full,RCP45_dQ_CAM3Grid_Full_file)
RCP45_dQ_CAM3Grid_Full_file.close()

wAVD_dQ_CAM3Grid_Full_file = open("wAVD_dlnQ_CAM3Grid.pickle","wb")
pickle.dump(wAVD_dQ_CAM3Grid_Full,wAVD_dQ_CAM3Grid_Full_file)
wAVD_dQ_CAM3Grid_Full_file.close()
