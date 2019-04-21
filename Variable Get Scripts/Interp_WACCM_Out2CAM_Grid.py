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


#load in climate change Responses
RCP45_dST = pickle.load(open(\
    "Future Projection CC Responses/RCP45_dST.pickle","rb"))
wAVD_dST = pickle.load(open(\
    "Future Projection CC Responses/wAVD_dST.pickle","rb"))

RCP45_dALB = pickle.load(open(\
    "Future Projection CC Responses/RCP45_dALB.pickle","rb"))
wAVD_dALB = pickle.load(open(\
    "Future Projection CC Responses/wAVD_dALB.pickle","rb"))

RCP45_dFLUT = pickle.load(open(\
    "Future Projection CC Responses/RCP45_dFLUT.pickle","rb"))
wAVD_dFLUT = pickle.load(open(\
    "Future Projection CC Responses/wAVD_dFLUT.pickle","rb"))

RCP45_dFLUTC = pickle.load(open(\
    "Future Projection CC Responses/RCP45_dFLUTC.pickle","rb"))
wAVD_dFLUTC = pickle.load(open(\
    "Future Projection CC Responses/wAVD_dFLUTC.pickle","rb"))

RCP45_dFSNT = pickle.load(open(\
    "Future Projection CC Responses/RCP45_dFSNT.pickle","rb"))
wAVD_dFSNT = pickle.load(open(\
    "Future Projection CC Responses/wAVD_dFSNT.pickle","rb"))

RCP45_dFSNTC = pickle.load(open(\
    "Future Projection CC Responses/RCP45_dFSNTC.pickle","rb"))
wAVD_dFSNTC = pickle.load(open(\
    "Future Projection CC Responses/wAVD_dFSNTC.pickle","rb"))

RCP45_dAT = pickle.load(open(\
    "Future Projection CC Responses/RCP45_dAT.pickle","rb"))
wAVD_dAT = pickle.load(open(\
    "Future Projection CC Responses/wAVD_dAT.pickle","rb"))

RCP45_dQ = pickle.load(open(\
    "Future Projection CC Responses/RCP45_dlnQ.pickle","rb"))
wAVD_dQ = pickle.load(open(\
    "Future Projection CC Responses/wAVD_dlnQ.pickle","rb"))


RCP45_dAT_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
WACCM4_plevs_RCP45_camHgrid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
RCP45_dST_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            RCP45_dAT_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),RCP45_dAT[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')
            WACCM4_plevs_RCP45_camHgrid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),WACCM4_plevs_RCP45[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')
            RCP45_dST_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),RCP45_dST[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')


RCP45_dAT_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    RCP45_dAT_Tprof=np.squeeze(RCP45_dAT_CAM3Grid[e,month,:,x,y])
                    RCP45_dAT_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_RCP45_camHgrid[e,month,:,x,y]),\
                    RCP45_dAT_Tprof,CMIP_plevs,method="nearest")


wAVD_dAT_CAM3Grid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
WACCM4_plevs_wAVD_camHgrid = np.zeros((5,12,47,len(CAM3_Lon),len(CAM3_Lat)))
wAVD_dST_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for z in range(47):
        for e in range(5):
            wAVD_dAT_CAM3Grid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),wAVD_dAT[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')
            WACCM4_plevs_wAVD_camHgrid[e,month,z,:,:]=griddata(\
                                    (wx.flatten(),wy.flatten()),WACCM4_plevs_wAVD[e,month,z,:,:].T.flatten(),\
                                    (cx,cy),method='linear')
            wAVD_dST_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),wAVD_dST[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')


wAVD_dAT_CAM3Grid_Full = np.zeros((5,12,17,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for x in range(len(CAM3_Lon)):
            for y in range(len(CAM3_Lat)):
                for e in range(5):
                    wAVD_dAT_Tprof=np.squeeze(wAVD_dAT_CAM3Grid[e,month,:,x,y])
                    wAVD_dAT_CAM3Grid_Full[e,month,:,x,y]=\
                    griddata(np.squeeze(WACCM4_plevs_wAVD_camHgrid[e,month,:,x,y]),\
                    wAVD_dAT_Tprof,CMIP_plevs,method="nearest")


#save all interpolated arrays to new files

RCP45_dAT_CAM3Grid_Full_file = open("RCP45_dAT_CAM3Grid.pickle","wb")
pickle.dump(RCP45_dAT_CAM3Grid_Full,RCP45_dAT_CAM3Grid_Full_file)
RCP45_dAT_CAM3Grid_Full_file.close()

wAVD_dAT_CAM3Grid_Full_file = open("wAVD_dAT_CAM3Grid.pickle","wb")
pickle.dump(wAVD_dAT_CAM3Grid_Full,wAVD_dAT_CAM3Grid_Full_file)
wAVD_dAT_CAM3Grid_Full_file.close()

#save arrays to files

RCP45_dST_CAM3Grid_file = open("RCP45_dST_CAM3Grid.pickle","wb")
pickle.dump(RCP45_dST_CAM3Grid,RCP45_dST_CAM3Grid_file)
RCP45_dST_CAM3Grid_file.close()

wAVD_dST_CAM3Grid_file = open("wAVD_dST_CAM3Grid.pickle","wb")
pickle.dump(wAVD_dST_CAM3Grid,wAVD_dST_CAM3Grid_file)
wAVD_dST_CAM3Grid_file.close()



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


#albedo

RCP45_dALB_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for e in range(5):
        RCP45_dALB_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),RCP45_dALB[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')

wAVD_dALB_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for e in range(5):
        wAVD_dALB_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),wAVD_dALB[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')

#save arrays to files

RCP45_dALB_CAM3Grid_file = open("RCP45_dALB_CAM3Grid.pickle","wb")
pickle.dump(RCP45_dALB_CAM3Grid,RCP45_dALB_CAM3Grid_file)
RCP45_dALB_CAM3Grid_file.close()

wAVD_dALB_CAM3Grid_file = open("wAVD_dALB_CAM3Grid.pickle","wb")
pickle.dump(wAVD_dALB_CAM3Grid,wAVD_dALB_CAM3Grid_file)
wAVD_dALB_CAM3Grid_file.close()




#FLUT

RCP45_dFLUT_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for e in range(5):
        RCP45_dFLUT_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),RCP45_dFLUT[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')

wAVD_dFLUT_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for e in range(5):
        wAVD_dFLUT_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),wAVD_dFLUT[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')

#save arrays to files
RCP45_dFLUT_CAM3Grid_file = open("RCP45_dFLUT_CAM3Grid.pickle","wb")
pickle.dump(RCP45_dFLUT_CAM3Grid,RCP45_dFLUT_CAM3Grid_file)
RCP45_dFLUT_CAM3Grid_file.close()

wAVD_dFLUT_CAM3Grid_file = open("wAVD_dFLUT_CAM3Grid.pickle","wb")
pickle.dump(wAVD_dFLUT_CAM3Grid,wAVD_dFLUT_CAM3Grid_file)
wAVD_dFLUT_CAM3Grid_file.close()


#dFLUTC

RCP45_dFLUTC_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for e in range(5):
        RCP45_dFLUTC_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),RCP45_dFLUTC[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')

wAVD_dFLUTC_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for e in range(5):
        wAVD_dFLUTC_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),wAVD_dFLUTC[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')

#save arrays to files
RCP45_dFLUTC_CAM3Grid_file = open("RCP45_dFLUTC_CAM3Grid.pickle","wb")
pickle.dump(RCP45_dFLUTC_CAM3Grid,RCP45_dFLUTC_CAM3Grid_file)
RCP45_dFLUTC_CAM3Grid_file.close()

wAVD_dFLUTC_CAM3Grid_file = open("wAVD_dFLUTC_CAM3Grid.pickle","wb")
pickle.dump(wAVD_dFLUTC_CAM3Grid,wAVD_dFLUTC_CAM3Grid_file)
wAVD_dFLUTC_CAM3Grid_file.close()





#dFSNT
RCP45_dFSNT_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for e in range(5):
        RCP45_dFSNT_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),RCP45_dFSNT[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')

wAVD_dFSNT_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for e in range(5):
        wAVD_dFSNT_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),wAVD_dFSNT[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')

#save arrays to files

RCP45_dFSNT_CAM3Grid_file = open("RCP45_dFSNT_CAM3Grid.pickle","wb")
pickle.dump(RCP45_dFSNT_CAM3Grid,RCP45_dFSNT_CAM3Grid_file)
RCP45_dFSNT_CAM3Grid_file.close()

wAVD_dFSNT_CAM3Grid_file = open("wAVD_dFSNT_CAM3Grid.pickle","wb")
pickle.dump(wAVD_dFSNT_CAM3Grid,wAVD_dFSNT_CAM3Grid_file)
wAVD_dFSNT_CAM3Grid_file.close()


#dFSNTC
RCP45_dFSNTC_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for e in range(5):
        RCP45_dFSNTC_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),RCP45_dFSNTC[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')

wAVD_dFSNTC_CAM3Grid = np.zeros((5,12,len(CAM3_Lon),len(CAM3_Lat)))
for month in range(12):
    for e in range(5):
        wAVD_dFSNTC_CAM3Grid[e,month,:,:]=griddata(\
                                (wx.flatten(),wy.flatten()),wAVD_dFSNTC[e,month,:,:].T.flatten(),\
                                (cx,cy),method='linear')




RCP45_dFSNTC_CAM3Grid_file = open("RCP45_dFSNTC_CAM3Grid.pickle","wb")
pickle.dump(RCP45_dFSNTC_CAM3Grid,RCP45_dFSNTC_CAM3Grid_file)
RCP45_dFSNTC_CAM3Grid_file.close()

wAVD_dFSNTC_CAM3Grid_file = open("wAVD_dFSNTC_CAM3Grid.pickle","wb")
pickle.dump(wAVD_dFSNTC_CAM3Grid,wAVD_dFSNTC_CAM3Grid_file)
wAVD_dFSNTC_CAM3Grid_file.close()
