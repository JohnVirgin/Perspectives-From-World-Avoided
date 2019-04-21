#!/usr/bin/env python

## Import Packages
import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import pickle
import Area_Avg
from netCDF4 import Dataset
import warnings
warnings.filterwarnings('ignore')

## Read in all WACCM Data from all 9 ensemble members and calculate decadal averages, then regrid to modified WACCM resolution


Months = ['01','02','03','04','05','06','07','08','09','10','11','12']
Sims = ['001','002','003','004','005']
Sims_alt_name = ['001','002','003','001','002']
Years = list(map(str, (range(2006,2066))))

LLL = Dataset('/export/data/waccmdata/b40.rcp4_5.2deg.wcm.001/atm/b40.rcp4_5.2deg.wcm.001.2005.01.nc')
Lat = np.squeeze(LLL.variables['lat'])
Arctic_Lat = Lat[80:]
Lon = np.squeeze(LLL.variables['lon'])
LLL.close()

#establish full directories for each run

datadir_RCP45_001 = '/export/data/waccmdata/b40.rcp4_5.2deg.wcm.001/atm/b40.rcp4_5.2deg.wcm.'
datadir_RCP45_002 = '/export/data/waccmdata/b40.rcp4_5.2deg.wcm.002/atm/b40.rcp4_5.2deg.wcm.'
datadir_RCP45_003 = '/export/data/waccmdata/b40.rcp4_5.2deg.wcm.003/atm/b40.rcp4_5.2deg.wcm.'
datadir_RCP45_004 = '/export/data/waccmdata/b40.rcp4_5.2deg.wcm.KLS.001/atm/b40.rcp4_5.2deg.wcm.KLS.'
datadir_RCP45_005 = '/export/data/waccmdata/b40.rcp4_5.2deg.wcm.KLS.002/atm/b40.rcp4_5.2deg.wcm.KLS.'
datadir_wAVD_001 = '/export/data/waccmdata/b.e10.BRCP45WCN.f19_g16.waLMP.001/atm/b.e10.BRCP45WCN.f19_g16.waLMP.'
datadir_wAVD_002 = '/export/data/waccmdata/b.e10.BRCP45WCN.f19_g16.waLMP.002/atm/b.e10.BRCP45WCN.f19_g16.waLMP.'
datadir_wAVD_003 = '/export/data/waccmdata/b.e10.BRCP45WCN.f19_g16.waLMP.003/atm/b.e10.BRCP45WCN.f19_g16.waLMP.'
datadir_wAVD_004 = '/export/data/waccmdata/b40.rcp4_5.2deg.wcm.WAVD.KLS.001/atm/b40.rcp4_5.2deg.wcm.WAVD.KLS.'
datadir_wAVD_005 = '/export/data/waccmdata/b40.rcp4_5.2deg.wcm.WAVD.KLS.002/atm/b40.rcp4_5.2deg.wcm.WAVD.KLS.'

Datadir_RCP45 = [datadir_RCP45_001,datadir_RCP45_002,datadir_RCP45_003,\
                 datadir_RCP45_004,datadir_RCP45_005]
Datadir_wAVD = [datadir_wAVD_001,datadir_wAVD_002,datadir_wAVD_003,\
                datadir_wAVD_004,datadir_wAVD_005]

#SURFACE TEMPERATURE

RCP45_ST = dict()
for s in range(len(Sims)):
    RCP45_ST[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_ST[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_ST[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['TREFHT'])
            RCP45_nc.close()

wAVD_ST = dict()
for s in range(len(Sims)):
    wAVD_ST[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_ST[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_ST[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['TREFHT'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_ST_Panel = pd.Panel.from_dict(RCP45_ST, orient = 'items')
wAVD_ST_Panel = pd.Panel.from_dict(wAVD_ST, orient = 'items')

RCP45_ST_1DA = RCP45_ST_Panel['001':'005','2006':'2015','01':'12']
RCP45_ST_1DA = np.asarray(RCP45_ST_1DA)
RCP45_ST_1DA = np.mean(RCP45_ST_1DA[:,:,:], axis=(1,))

RCP45_ST_2DA = RCP45_ST_Panel['001':'005','2056':'2065','01':'12']
RCP45_ST_2DA = np.asarray(RCP45_ST_2DA)
RCP45_ST_2DA = np.mean(RCP45_ST_2DA[:,:,:], axis=(1,))

RCP45_ST_FDec = np.zeros([5,12,96,144])
RCP45_ST_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    RCP45_ST_FDec[i,:,:,:] = np.stack(RCP45_ST_1DA[i,:])
    RCP45_ST_SDec[i,:,:,:] = np.stack(RCP45_ST_2DA[i,:])

RCP45_dST = RCP45_ST_SDec-RCP45_ST_FDec

wAVD_ST_1DA = wAVD_ST_Panel['001':'005','2006':'2015','01':'12']
wAVD_ST_1DA = np.asarray(wAVD_ST_1DA)
wAVD_ST_1DA = np.mean(wAVD_ST_1DA[:,:,:], axis=(1,))

wAVD_ST_2DA = wAVD_ST_Panel['001':'005','2056':'2065','01':'12']
wAVD_ST_2DA = np.asarray(wAVD_ST_2DA)
wAVD_ST_2DA = np.mean(wAVD_ST_2DA[:,:,:], axis=(1,))

wAVD_ST_FDec = np.zeros([5,12,96,144])
wAVD_ST_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    wAVD_ST_FDec[i,:,:,:] = np.stack(wAVD_ST_1DA[i,:])
    wAVD_ST_SDec[i,:,:,:] = np.stack(wAVD_ST_2DA[i,:])

wAVD_dST = wAVD_ST_SDec-wAVD_ST_FDec

RCP45_Delta_ST_file = open("RCP45_dST.pickle",'wb')
pickle.dump(RCP45_dST,RCP45_Delta_ST_file)
RCP45_Delta_ST_file.close()

wAVD_Delta_ST_file = open("wAVD_dST.pickle",'wb')
pickle.dump(wAVD_dST,wAVD_Delta_ST_file)
wAVD_Delta_ST_file.close()

#AIR TEMPERATURE

RCP45_AT = dict()
for s in range(len(Sims)):
    RCP45_AT[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_AT[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_AT[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['T'])
            RCP45_nc.close()

wAVD_AT = dict()
for s in range(len(Sims)):
    wAVD_AT[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_AT[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_AT[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['T'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_AT_Panel = pd.Panel.from_dict(RCP45_AT, orient = 'items')
wAVD_AT_Panel = pd.Panel.from_dict(wAVD_AT, orient = 'items')

RCP45_AT_1DA = RCP45_AT_Panel['001':'005','2006':'2015','01':'12']
RCP45_AT_1DA = np.asarray(RCP45_AT_1DA)
RCP45_AT_1DA = np.mean(RCP45_AT_1DA[:,:,:], axis=(1,))

RCP45_AT_2DA = RCP45_AT_Panel['001':'005','2056':'2065','01':'12']
RCP45_AT_2DA = np.asarray(RCP45_AT_2DA)
RCP45_AT_2DA = np.mean(RCP45_AT_2DA[:,:,:], axis=(1,))

RCP45_AT_FDec = np.zeros([5,12,47,96,144])
RCP45_AT_SDec = np.zeros([5,12,47,96,144])

for i in range(len(Sims)):
    RCP45_AT_FDec[i,:,:,:,:] = np.stack(RCP45_AT_1DA[i,:])
    RCP45_AT_SDec[i,:,:,:,:] = np.stack(RCP45_AT_2DA[i,:])

RCP45_dAT = RCP45_AT_SDec-RCP45_AT_FDec

wAVD_AT_1DA = wAVD_AT_Panel['001':'005','2006':'2015','01':'12']
wAVD_AT_1DA = np.asarray(wAVD_AT_1DA)
wAVD_AT_1DA = np.mean(wAVD_AT_1DA[:,:,:], axis=(1,))

wAVD_AT_2DA = wAVD_AT_Panel['001':'005','2056':'2065','01':'12']
wAVD_AT_2DA = np.asarray(wAVD_AT_2DA)
wAVD_AT_2DA = np.mean(wAVD_AT_2DA[:,:,:], axis=(1,))

wAVD_AT_FDec = np.zeros([5,12,47,96,144])
wAVD_AT_SDec = np.zeros([5,12,47,96,144])

for i in range(len(Sims)):
    wAVD_AT_FDec[i,:,:,:,:] = np.stack(wAVD_AT_1DA[i,:])
    wAVD_AT_SDec[i,:,:,:,:] = np.stack(wAVD_AT_2DA[i,:])

wAVD_dAT = wAVD_AT_SDec-wAVD_AT_FDec

RCP45_Delta_AT_file = open("RCP45_dAT.pickle",'wb')
pickle.dump(RCP45_dAT,RCP45_Delta_AT_file)
RCP45_Delta_AT_file.close()

wAVD_Delta_AT_file = open("wAVD_dAT.pickle",'wb')
pickle.dump(wAVD_dAT,wAVD_Delta_AT_file)
wAVD_Delta_AT_file.close()

#WATER VAPOUR

RCP45_Q = dict()
for s in range(len(Sims)):
    RCP45_Q[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_Q[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_Q[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['Q'])
            RCP45_nc.close()

wAVD_Q = dict()
for s in range(len(Sims)):
    wAVD_Q[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_Q[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_Q[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['Q'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_Q_Panel = pd.Panel.from_dict(RCP45_Q, orient = 'items')
wAVD_Q_Panel = pd.Panel.from_dict(wAVD_Q, orient = 'items')

RCP45_Q_1DA = RCP45_Q_Panel['001':'005','2006':'2015','01':'12']
RCP45_Q_1DA = np.asarray(RCP45_Q_1DA)
RCP45_Q_1DA = np.mean(RCP45_Q_1DA[:,:,:], axis=(1,))

RCP45_Q_2DA = RCP45_Q_Panel['001':'005','2056':'2065','01':'12']
RCP45_Q_2DA = np.asarray(RCP45_Q_2DA)
RCP45_Q_2DA = np.mean(RCP45_Q_2DA[:,:,:], axis=(1,))

RCP45_Q_FDec = np.zeros([5,12,47,96,144])
RCP45_Q_SDec = np.zeros([5,12,47,96,144])

for i in range(len(Sims)):
    RCP45_Q_FDec[i,:,:,:,:] = np.stack(RCP45_Q_1DA[i,:])
    RCP45_Q_SDec[i,:,:,:,:] = np.stack(RCP45_Q_2DA[i,:])

RCP45_dQ = RCP45_Q_SDec-RCP45_Q_FDec

wAVD_Q_1DA = wAVD_Q_Panel['001':'005','2006':'2015','01':'12']
wAVD_Q_1DA = np.asarray(wAVD_Q_1DA)
wAVD_Q_1DA = np.mean(wAVD_Q_1DA[:,:,:], axis=(1,))

wAVD_Q_2DA = wAVD_Q_Panel['001':'005','2056':'2065','01':'12']
wAVD_Q_2DA = np.asarray(wAVD_Q_2DA)
wAVD_Q_2DA = np.mean(wAVD_Q_2DA[:,:,:], axis=(1,))

wAVD_Q_FDec = np.zeros([5,12,47,96,144])
wAVD_Q_SDec = np.zeros([5,12,47,96,144])

for i in range(len(Sims)):
    wAVD_Q_FDec[i,:,:,:,:] = np.stack(wAVD_Q_1DA[i,:])
    wAVD_Q_SDec[i,:,:,:,:] = np.stack(wAVD_Q_2DA[i,:])

wAVD_dQ = wAVD_Q_SDec-wAVD_Q_FDec

RCP45_Delta_Q_file = open("RCP45_dQ.pickle",'wb')
pickle.dump(RCP45_dQ,RCP45_Delta_Q_file)
RCP45_Delta_Q_file.close()

wAVD_Delta_Q_file = open("wAVD_dQ.pickle",'wb')
pickle.dump(wAVD_dQ,wAVD_Delta_Q_file)
wAVD_Delta_Q_file.close()

#NET UPWELLING FLUX AT THE TOA

RCP45_FLUT = dict()
for s in range(len(Sims)):
    RCP45_FLUT[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_FLUT[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_FLUT[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['FLUT'])
            RCP45_nc.close()

wAVD_FLUT = dict()
for s in range(len(Sims)):
    wAVD_FLUT[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_FLUT[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_FLUT[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['FLUT'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_FLUT_Panel = pd.Panel.from_dict(RCP45_FLUT, orient = 'items')
wAVD_FLUT_Panel = pd.Panel.from_dict(wAVD_FLUT, orient = 'items')

RCP45_FLUT_1DA = RCP45_FLUT_Panel['001':'005','2006':'2015','01':'12']
RCP45_FLUT_1DA = np.asarray(RCP45_FLUT_1DA)
RCP45_FLUT_1DA = np.mean(RCP45_FLUT_1DA[:,:,:], axis=(1,))

RCP45_FLUT_2DA = RCP45_FLUT_Panel['001':'005','2056':'2065','01':'12']
RCP45_FLUT_2DA = np.asarray(RCP45_FLUT_2DA)
RCP45_FLUT_2DA = np.mean(RCP45_FLUT_2DA[:,:,:], axis=(1,))

RCP45_FLUT_FDec = np.zeros([5,12,96,144])
RCP45_FLUT_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    RCP45_FLUT_FDec[i,:,:,:] = np.stack(RCP45_FLUT_1DA[i,:])
    RCP45_FLUT_SDec[i,:,:,:] = np.stack(RCP45_FLUT_2DA[i,:])

RCP45_dFLUT = RCP45_FLUT_SDec-RCP45_FLUT_FDec

wAVD_FLUT_1DA = wAVD_FLUT_Panel['001':'005','2006':'2015','01':'12']
wAVD_FLUT_1DA = np.asarray(wAVD_FLUT_1DA)
wAVD_FLUT_1DA = np.mean(wAVD_FLUT_1DA[:,:,:], axis=(1,))

wAVD_FLUT_2DA = wAVD_FLUT_Panel['001':'005','2056':'2065','01':'12']
wAVD_FLUT_2DA = np.asarray(wAVD_FLUT_2DA)
wAVD_FLUT_2DA = np.mean(wAVD_FLUT_2DA[:,:,:], axis=(1,))

wAVD_FLUT_FDec = np.zeros([5,12,96,144])
wAVD_FLUT_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    wAVD_FLUT_FDec[i,:,:,:] = np.stack(wAVD_FLUT_1DA[i,:])
    wAVD_FLUT_SDec[i,:,:,:] = np.stack(wAVD_FLUT_2DA[i,:])

wAVD_dFLUT = wAVD_FLUT_SDec-wAVD_FLUT_FDec

RCP45_Delta_FLUT_file = open("RCP45_dFLUT.pickle",'wb')
pickle.dump(RCP45_dFLUT,RCP45_Delta_FLUT_file)
RCP45_Delta_FLUT_file.close()

wAVD_Delta_FLUT_file = open("wAVD_dFLUT.pickle",'wb')
pickle.dump(wAVD_dFLUT,wAVD_Delta_FLUT_file)
wAVD_Delta_FLUT_file.close()

#CLEAR SKY SAME THING

RCP45_FLUTC = dict()
for s in range(len(Sims)):
    RCP45_FLUTC[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_FLUTC[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_FLUTC[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['FLUTC'])
            RCP45_nc.close()

wAVD_FLUTC = dict()
for s in range(len(Sims)):
    wAVD_FLUTC[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_FLUTC[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_FLUTC[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['FLUTC'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_FLUTC_Panel = pd.Panel.from_dict(RCP45_FLUTC, orient = 'items')
wAVD_FLUTC_Panel = pd.Panel.from_dict(wAVD_FLUTC, orient = 'items')

RCP45_FLUTC_1DA = RCP45_FLUTC_Panel['001':'005','2006':'2015','01':'12']
RCP45_FLUTC_1DA = np.asarray(RCP45_FLUTC_1DA)
RCP45_FLUTC_1DA = np.mean(RCP45_FLUTC_1DA[:,:,:], axis=(1,))

RCP45_FLUTC_2DA = RCP45_FLUTC_Panel['001':'005','2056':'2065','01':'12']
RCP45_FLUTC_2DA = np.asarray(RCP45_FLUTC_2DA)
RCP45_FLUTC_2DA = np.mean(RCP45_FLUTC_2DA[:,:,:], axis=(1,))

RCP45_FLUTC_FDec = np.zeros([5,12,96,144])
RCP45_FLUTC_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    RCP45_FLUTC_FDec[i,:,:,:] = np.stack(RCP45_FLUTC_1DA[i,:])
    RCP45_FLUTC_SDec[i,:,:,:] = np.stack(RCP45_FLUTC_2DA[i,:])

RCP45_dFLUTC = RCP45_FLUTC_SDec-RCP45_FLUTC_FDec

wAVD_FLUTC_1DA = wAVD_FLUTC_Panel['001':'005','2006':'2015','01':'12']
wAVD_FLUTC_1DA = np.asarray(wAVD_FLUTC_1DA)
wAVD_FLUTC_1DA = np.mean(wAVD_FLUTC_1DA[:,:,:], axis=(1,))

wAVD_FLUTC_2DA = wAVD_FLUTC_Panel['001':'005','2056':'2065','01':'12']
wAVD_FLUTC_2DA = np.asarray(wAVD_FLUTC_2DA)
wAVD_FLUTC_2DA = np.mean(wAVD_FLUTC_2DA[:,:,:], axis=(1,))

wAVD_FLUTC_FDec = np.zeros([5,12,96,144])
wAVD_FLUTC_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    wAVD_FLUTC_FDec[i,:,:,:] = np.stack(wAVD_FLUTC_1DA[i,:])
    wAVD_FLUTC_SDec[i,:,:,:] = np.stack(wAVD_FLUTC_2DA[i,:])

wAVD_dFLUTC = wAVD_FLUTC_SDec-wAVD_FLUTC_FDec

RCP45_Delta_FLUTC_file = open("RCP45_dFLUTC.pickle",'wb')
pickle.dump(RCP45_dFLUTC,RCP45_Delta_FLUTC_file)
RCP45_Delta_FLUTC_file.close()

wAVD_Delta_FLUTC_file = open("wAVD_dFLUTC.pickle",'wb')
pickle.dump(wAVD_dFLUTC,wAVD_Delta_FLUTC_file)
wAVD_Delta_FLUTC_file.close()

#NET SOLAR FLUX AT THE TOA

RCP45_FSNT = dict()
for s in range(len(Sims)):
    RCP45_FSNT[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_FSNT[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_FSNT[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['FSNT'])
            RCP45_nc.close()

wAVD_FSNT = dict()
for s in range(len(Sims)):
    wAVD_FSNT[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_FSNT[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_FSNT[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['FSNT'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_FSNT_Panel = pd.Panel.from_dict(RCP45_FSNT, orient = 'items')
wAVD_FSNT_Panel = pd.Panel.from_dict(wAVD_FSNT, orient = 'items')

RCP45_FSNT_1DA = RCP45_FSNT_Panel['001':'005','2006':'2015','01':'12']
RCP45_FSNT_1DA = np.asarray(RCP45_FSNT_1DA)
RCP45_FSNT_1DA = np.mean(RCP45_FSNT_1DA[:,:,:], axis=(1,))

RCP45_FSNT_2DA = RCP45_FSNT_Panel['001':'005','2056':'2065','01':'12']
RCP45_FSNT_2DA = np.asarray(RCP45_FSNT_2DA)
RCP45_FSNT_2DA = np.mean(RCP45_FSNT_2DA[:,:,:], axis=(1,))

RCP45_FSNT_FDec = np.zeros([5,12,96,144])
RCP45_FSNT_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    RCP45_FSNT_FDec[i,:,:,:] = np.stack(RCP45_FSNT_1DA[i,:])
    RCP45_FSNT_SDec[i,:,:,:] = np.stack(RCP45_FSNT_2DA[i,:])

RCP45_dFSNT = RCP45_FSNT_SDec-RCP45_FSNT_FDec

wAVD_FSNT_1DA = wAVD_FSNT_Panel['001':'005','2006':'2015','01':'12']
wAVD_FSNT_1DA = np.asarray(wAVD_FSNT_1DA)
wAVD_FSNT_1DA = np.mean(wAVD_FSNT_1DA[:,:,:], axis=(1,))

wAVD_FSNT_2DA = wAVD_FSNT_Panel['001':'005','2056':'2065','01':'12']
wAVD_FSNT_2DA = np.asarray(wAVD_FSNT_2DA)
wAVD_FSNT_2DA = np.mean(wAVD_FSNT_2DA[:,:,:], axis=(1,))

wAVD_FSNT_FDec = np.zeros([5,12,96,144])
wAVD_FSNT_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    wAVD_FSNT_FDec[i,:,:,:] = np.stack(wAVD_FSNT_1DA[i,:])
    wAVD_FSNT_SDec[i,:,:,:] = np.stack(wAVD_FSNT_2DA[i,:])

wAVD_dFSNT = wAVD_FSNT_SDec-wAVD_FSNT_FDec

RCP45_Delta_FSNT_file = open("RCP45_dFSNT.pickle",'wb')
pickle.dump(RCP45_dFSNT,RCP45_Delta_FSNT_file)
RCP45_Delta_FSNT_file.close()

wAVD_Delta_FSNT_file = open("wAVD_dFSNT.pickle",'wb')
pickle.dump(wAVD_dFSNT,wAVD_Delta_FSNT_file)
wAVD_Delta_FSNT_file.close()

#CLEAR SKY

RCP45_FSNTC = dict()
for s in range(len(Sims)):
    RCP45_FSNTC[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_FSNTC[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_FSNTC[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['FSNTC'])
            RCP45_nc.close()

wAVD_FSNTC = dict()
for s in range(len(Sims)):
    wAVD_FSNTC[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_FSNTC[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_FSNTC[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['FSNTC'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_FSNTC_Panel = pd.Panel.from_dict(RCP45_FSNTC, orient = 'items')
wAVD_FSNTC_Panel = pd.Panel.from_dict(wAVD_FSNTC, orient = 'items')

RCP45_FSNTC_1DA = RCP45_FSNTC_Panel['001':'005','2006':'2015','01':'12']
RCP45_FSNTC_1DA = np.asarray(RCP45_FSNTC_1DA)
RCP45_FSNTC_1DA = np.mean(RCP45_FSNTC_1DA[:,:,:], axis=(1,))

RCP45_FSNTC_2DA = RCP45_FSNTC_Panel['001':'005','2056':'2065','01':'12']
RCP45_FSNTC_2DA = np.asarray(RCP45_FSNTC_2DA)
RCP45_FSNTC_2DA = np.mean(RCP45_FSNTC_2DA[:,:,:], axis=(1,))

RCP45_FSNTC_FDec = np.zeros([5,12,96,144])
RCP45_FSNTC_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    RCP45_FSNTC_FDec[i,:,:,:] = np.stack(RCP45_FSNTC_1DA[i,:])
    RCP45_FSNTC_SDec[i,:,:,:] = np.stack(RCP45_FSNTC_2DA[i,:])

RCP45_dFSNTC = RCP45_FSNTC_SDec-RCP45_FSNTC_FDec

wAVD_FSNTC_1DA = wAVD_FSNTC_Panel['001':'005','2006':'2015','01':'12']
wAVD_FSNTC_1DA = np.asarray(wAVD_FSNTC_1DA)
wAVD_FSNTC_1DA = np.mean(wAVD_FSNTC_1DA[:,:,:], axis=(1,))

wAVD_FSNTC_2DA = wAVD_FSNTC_Panel['001':'005','2056':'2065','01':'12']
wAVD_FSNTC_2DA = np.asarray(wAVD_FSNTC_2DA)
wAVD_FSNTC_2DA = np.mean(wAVD_FSNTC_2DA[:,:,:], axis=(1,))

wAVD_FSNTC_FDec = np.zeros([5,12,96,144])
wAVD_FSNTC_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    wAVD_FSNTC_FDec[i,:,:,:] = np.stack(wAVD_FSNTC_1DA[i,:])
    wAVD_FSNTC_SDec[i,:,:,:] = np.stack(wAVD_FSNTC_2DA[i,:])

wAVD_dFSNTC = wAVD_FSNTC_SDec-wAVD_FSNTC_FDec

RCP45_Delta_FSNTC_file = open("RCP45_dFSNTC.pickle",'wb')
pickle.dump(RCP45_dFSNTC,RCP45_Delta_FSNTC_file)
RCP45_Delta_FSNTC_file.close()

wAVD_Delta_FSNTC_file = open("wAVD_dFSNTC.pickle",'wb')
pickle.dump(wAVD_dFSNTC,wAVD_Delta_FSNTC_file)
wAVD_Delta_FSNTC_file.close()



#FSNS

RCP45_FSNS = dict()
for s in range(len(Sims)):
    RCP45_FSNS[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_FSNS[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_FSNS[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['FSNS'])
            RCP45_nc.close()

wAVD_FSNS = dict()
for s in range(len(Sims)):
    wAVD_FSNS[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_FSNS[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_FSNS[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['FSNS'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_FSNS_Panel = pd.Panel.from_dict(RCP45_FSNS, orient = 'items')
wAVD_FSNS_Panel = pd.Panel.from_dict(wAVD_FSNS, orient = 'items')

RCP45_FSNS_1DA = RCP45_FSNS_Panel['001':'005','2006':'2015','01':'12']
RCP45_FSNS_1DA = np.asarray(RCP45_FSNS_1DA)
RCP45_FSNS_1DA = np.mean(RCP45_FSNS_1DA[:,:,:], axis=(1,))

RCP45_FSNS_2DA = RCP45_FSNS_Panel['001':'005','2056':'2065','01':'12']
RCP45_FSNS_2DA = np.asarray(RCP45_FSNS_2DA)
RCP45_FSNS_2DA = np.mean(RCP45_FSNS_2DA[:,:,:], axis=(1,))

RCP45_FSNS_FDec = np.zeros([5,12,96,144])
RCP45_FSNS_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    RCP45_FSNS_FDec[i,:,:,:] = np.stack(RCP45_FSNS_1DA[i,:])
    RCP45_FSNS_SDec[i,:,:,:] = np.stack(RCP45_FSNS_2DA[i,:])

wAVD_FSNS_1DA = wAVD_FSNS_Panel['001':'005','2006':'2015','01':'12']
wAVD_FSNS_1DA = np.asarray(wAVD_FSNS_1DA)
wAVD_FSNS_1DA = np.mean(wAVD_FSNS_1DA[:,:,:], axis=(1,))

wAVD_FSNS_2DA = wAVD_FSNS_Panel['001':'005','2056':'2065','01':'12']
wAVD_FSNS_2DA = np.asarray(wAVD_FSNS_2DA)
wAVD_FSNS_2DA = np.mean(wAVD_FSNS_2DA[:,:,:], axis=(1,))

wAVD_FSNS_FDec = np.zeros([5,12,96,144])
wAVD_FSNS_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    wAVD_FSNS_FDec[i,:,:,:] = np.stack(wAVD_FSNS_1DA[i,:])
    wAVD_FSNS_SDec[i,:,:,:] = np.stack(wAVD_FSNS_2DA[i,:])


#FSDS

RCP45_FSDS = dict()
for s in range(len(Sims)):
    RCP45_FSDS[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_FSDS[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_FSDS[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['FSDS'])
            RCP45_nc.close()

wAVD_FSDS = dict()
for s in range(len(Sims)):
    wAVD_FSDS[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_FSDS[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_FSDS[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['FSDS'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_FSDS_Panel = pd.Panel.from_dict(RCP45_FSDS, orient = 'items')
wAVD_FSDS_Panel = pd.Panel.from_dict(wAVD_FSDS, orient = 'items')

RCP45_FSDS_1DA = RCP45_FSDS_Panel['001':'005','2006':'2015','01':'12']
RCP45_FSDS_1DA = np.asarray(RCP45_FSDS_1DA)
RCP45_FSDS_1DA = np.mean(RCP45_FSDS_1DA[:,:,:], axis=(1,))

RCP45_FSDS_2DA = RCP45_FSDS_Panel['001':'005','2056':'2065','01':'12']
RCP45_FSDS_2DA = np.asarray(RCP45_FSDS_2DA)
RCP45_FSDS_2DA = np.mean(RCP45_FSDS_2DA[:,:,:], axis=(1,))

RCP45_FSDS_FDec = np.zeros([5,12,96,144])
RCP45_FSDS_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    RCP45_FSDS_FDec[i,:,:,:] = np.stack(RCP45_FSDS_1DA[i,:])
    RCP45_FSDS_SDec[i,:,:,:] = np.stack(RCP45_FSDS_2DA[i,:])

wAVD_FSDS_1DA = wAVD_FSDS_Panel['001':'005','2006':'2015','01':'12']
wAVD_FSDS_1DA = np.asarray(wAVD_FSDS_1DA)
wAVD_FSDS_1DA = np.mean(wAVD_FSDS_1DA[:,:,:], axis=(1,))

wAVD_FSDS_2DA = wAVD_FSDS_Panel['001':'005','2056':'2065','01':'12']
wAVD_FSDS_2DA = np.asarray(wAVD_FSDS_2DA)
wAVD_FSDS_2DA = np.mean(wAVD_FSDS_2DA[:,:,:], axis=(1,))

wAVD_FSDS_FDec = np.zeros([5,12,96,144])
wAVD_FSDS_SDec = np.zeros([5,12,96,144])

for i in range(len(Sims)):
    wAVD_FSDS_FDec[i,:,:,:] = np.stack(wAVD_FSDS_1DA[i,:])
    wAVD_FSDS_SDec[i,:,:,:] = np.stack(wAVD_FSDS_2DA[i,:])


#Albedo Calculation

RCP45_Albedo_Decade1 = 1-(np.divide(RCP45_FSNS_FDec,RCP45_FSDS_FDec))
RCP45_Albedo_Decade2 = 1-(np.divide(RCP45_FSNS_SDec,RCP45_FSDS_SDec))
RCP45_dALB = RCP45_Albedo_Decade2-RCP45_Albedo_Decade1

wAVD_Albedo_Decade1 = 1-(np.divide(wAVD_FSNS_FDec,wAVD_FSDS_FDec))
wAVD_Albedo_Decade2 = 1-(np.divide(wAVD_FSNS_SDec,wAVD_FSDS_SDec))
wAVD_dALB = wAVD_Albedo_Decade2-wAVD_Albedo_Decade1

RCP45_dALB_file = open("RCP45_dALB.pickle","wb")
pickle.dump(RCP45_dALB,RCP45_dALB_file)
RCP45_dALB_file.close()

wAVD_dALB_file = open("wAVD_dALB.pickle","wb")
pickle.dump(wAVD_dALB,wAVD_dALB_file)
wAVD_dALB_file.close()
