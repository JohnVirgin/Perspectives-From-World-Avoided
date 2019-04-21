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

RCP45_Q1 = pickle.load(open(\
    "Vars for Specific Humidity Calculations/RCP45_Q_2006-2015.pickle",\
    "rb"))

RCP45_Q2 = pickle.load(open(\
    "Vars for Specific Humidity Calculations/RCP45_Q_2056-2065.pickle",\
    "rb"))

RCP45_Ta1 = pickle.load(open(\
    "Vars for Specific Humidity Calculations/RCP45_AT_2006-2015.pickle",\
    "rb"))

RCP45_Ta2 = pickle.load(open(\
    "Vars for Specific Humidity Calculations/RCP45_AT_2056-2065.pickle",\
    "rb"))

wAVD_Q1 = pickle.load(open(\
    "Vars for Specific Humidity Calculations/wAVD_Q_2006-2015.pickle",\
    "rb"))

wAVD_Q2 = pickle.load(open(\
    "Vars for Specific Humidity Calculations/wAVD_Q_2056-2065.pickle",\
    "rb"))

wAVD_Ta1 = pickle.load(open(\
    "Vars for Specific Humidity Calculations/wAVD_AT_2006-2015.pickle",\
    "rb"))

wAVD_Ta2 = pickle.load(open(\
    "Vars for Specific Humidity Calculations/wAVD_AT_2056-2065.pickle",\
    "rb"))

RCP45_Ta1_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
WACCM4_plevs_RCP45_camHgrid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            RCP45_Ta1_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),RCP45_Ta1[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')
            WACCM4_plevs_RCP45_camHgrid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),WACCM4_plevs_RCP45[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')


RCP45_Ta1_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    RCP45_Ta1_Tprof=np.squeeze(RCP45_Ta1_CAM3Grid[e,month,:,x,y])
                    RCP45_Ta1_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_RCP45_camHgrid[e,month,:,x,y]),\
                    RCP45_Ta1_Tprof,CMIP_plevs,method="nearest")


wAVD_Ta1_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
WACCM4_plevs_wAVD_camHgrid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            wAVD_Ta1_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),wAVD_Ta1[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')
            WACCM4_plevs_wAVD_camHgrid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),WACCM4_plevs_wAVD[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')


wAVD_Ta1_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    wAVD_Ta1_Tprof=np.squeeze(wAVD_Ta1_CAM3Grid[e,month,:,x,y])
                    wAVD_Ta1_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_wAVD_camHgrid[e,month,:,x,y]),\
                    wAVD_Ta1_Tprof,CMIP_plevs,method="nearest")



RCP45_Ta2_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            RCP45_Ta2_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),RCP45_Ta2[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')


RCP45_Ta2_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    RCP45_Ta2_Tprof=np.squeeze(RCP45_Ta2_CAM3Grid[e,month,:,x,y])
                    RCP45_Ta2_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_RCP45_camHgrid[e,month,:,x,y]),\
                    RCP45_Ta2_Tprof,CMIP_plevs,method="nearest")


wAVD_Ta2_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            wAVD_Ta2_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),wAVD_Ta2[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')


wAVD_Ta2_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    wAVD_Ta2_Tprof=np.squeeze(wAVD_Ta2_CAM3Grid[e,month,:,x,y])
                    wAVD_Ta2_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_wAVD_camHgrid[e,month,:,x,y]),\
                    wAVD_Ta2_Tprof,CMIP_plevs,method="nearest")





RCP45_Q1_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            RCP45_Q1_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),RCP45_Q1[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')


RCP45_Q1_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    RCP45_Q1_Tprof=np.squeeze(RCP45_Q1_CAM3Grid[e,month,:,x,y])
                    RCP45_Q1_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_RCP45_camHgrid[e,month,:,x,y]),\
                    RCP45_Q1_Tprof,CMIP_plevs,method="nearest")


wAVD_Q1_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            wAVD_Q1_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),wAVD_Q1[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')

wAVD_Q1_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    wAVD_Q1_Tprof=np.squeeze(wAVD_Q1_CAM3Grid[e,month,:,x,y])
                    wAVD_Q1_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_wAVD_camHgrid[e,month,:,x,y]),\
                    wAVD_Q1_Tprof,CMIP_plevs,method="nearest")



RCP45_Q2_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            RCP45_Q2_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),RCP45_Q2[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')


RCP45_Q2_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    RCP45_Q2_Tprof=np.squeeze(RCP45_Q2_CAM3Grid[e,month,:,x,y])
                    RCP45_Q2_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_RCP45_camHgrid[e,month,:,x,y]),\
                    RCP45_Q2_Tprof,CMIP_plevs,method="nearest")


wAVD_Q2_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            wAVD_Q2_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),wAVD_Q2[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')


wAVD_Q2_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    wAVD_Q2_Tprof=np.squeeze(wAVD_Q2_CAM3Grid[e,month,:,x,y])
                    wAVD_Q2_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_wAVD_camHgrid[e,month,:,x,y]),\
                    wAVD_Q2_Tprof,CMIP_plevs,method="nearest")


RCP45_Ta_Base_CAM3Grid_file = open('RCP45_Ta_2006-2015_CAM3Grid.pickle','wb')
pickle.dump(RCP45_Ta1_CAM3Grid_Full,RCP45_Ta_Base_CAM3Grid_file)
RCP45_Ta_Base_CAM3Grid_file.close()

RCP45_Ta_Response_CAM3Grid_file = open('RCP45_Ta_2056-2065_CAM3Grid.pickle','wb')
pickle.dump(RCP45_Ta2_CAM3Grid_Full,RCP45_Ta_Response_CAM3Grid_file)
RCP45_Ta_Response_CAM3Grid_file.close()

wAVD_Ta_Base_CAM3Grid_file = open('wAVD_Ta_2006-2015_CAM3Grid.pickle','wb')
pickle.dump(wAVD_Ta1_CAM3Grid_Full,wAVD_Ta_Base_CAM3Grid_file)
wAVD_Ta_Base_CAM3Grid_file.close()

wAVD_Ta_Response_CAM3Grid_file = open('wAVD_Ta_2056-2065_CAM3Grid.pickle','wb')
pickle.dump(wAVD_Ta2_CAM3Grid_Full,wAVD_Ta_Response_CAM3Grid_file)
wAVD_Ta_Response_CAM3Grid_file.close()

RCP45_Q_Base_CAM3Grid_file = open('RCP45_Q_2006-2015_CAM3Grid.pickle','wb')
pickle.dump(RCP45_Q1_CAM3Grid_Full,RCP45_Q_Base_CAM3Grid_file)
RCP45_Q_Base_CAM3Grid_file.close()

RCP45_Q_Response_CAM3Grid_file = open('RCP45_Q_2056-2065_CAM3Grid.pickle','wb')
pickle.dump(RCP45_Q2_CAM3Grid_Full,RCP45_Q_Response_CAM3Grid_file)
RCP45_Q_Response_CAM3Grid_file.close()

wAVD_Q_Base_CAM3Grid_file = open('wAVD_Q_2006-2015_CAM3Grid.pickle','wb')
pickle.dump(wAVD_Q1_CAM3Grid_Full,wAVD_Q_Base_CAM3Grid_file)
wAVD_Q_Base_CAM3Grid_file.close()

wAVD_Q_Response_CAM3Grid_file = open('wAVD_Q_2056-2065_CAM3Grid.pickle','wb')
pickle.dump(wAVD_Q2_CAM3Grid_Full,wAVD_Q_Response_CAM3Grid_file)
wAVD_Q_Response_CAM3Grid_file.close()
