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


RCP45_Q_FDec_file = open("RCP45_Q_2006-2015.pickle",'wb')
pickle.dump(RCP45_Q_FDec,RCP45_Q_FDec_file)
RCP45_Q_FDec_file.close()

RCP45_Q_SDec_file = open("RCP45_Q_2056-2065.pickle",'wb')
pickle.dump(RCP45_Q_SDec,RCP45_Q_SDec_file)
RCP45_Q_SDec_file.close()

wAVD_Q_FDec_file = open("wAVD_Q_2006-2015.pickle",'wb')
pickle.dump(wAVD_Q_FDec,wAVD_Q_FDec_file)
wAVD_Q_FDec_file.close()

wAVD_Q_SDec_file = open("wAVD_Q_2056-2065.pickle",'wb')
pickle.dump(wAVD_Q_SDec,wAVD_Q_SDec_file)
wAVD_Q_SDec_file.close()
