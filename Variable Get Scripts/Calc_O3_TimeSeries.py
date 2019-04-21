#!/usr/bin/env python

## Import Packages
import numpy as np
import pandas as pd
import scipy.stats as stats
import os
import pickle
from netCDF4 import Dataset
import Area_Avg as AA
import Jacks_Functions as JF
import warnings
warnings.filterwarnings('ignore')


Months = ['01','02','03','04','05','06','07','08','09','10','11','12']
Sims = ['001','002','003','004','005']
Sims_alt_name = ['001','002','003','001','002']
Years = list(map(str, (range(2006,2066))))

LLL = Dataset('/export/data/waccmdata/b40.rcp4_5.2deg.wcm.001/atm/b40.rcp4_5.2deg.wcm.001.2005.01.nc')
WACCM4_Lat = np.squeeze(LLL.variables['lat'])
WACCM4_Lon = np.squeeze(LLL.variables['lon'])
hyam = np.squeeze(LLL.variables['hyam']) #coef A for midpoints
hybm = np.squeeze(LLL.variables['hybm']) #coef B for midpoints
hyai = np.squeeze(LLL.variables['hyai']) #coef A for interfaces
hybi = np.squeeze(LLL.variables['hybi']) #coef B for interfaces
p0 = 1000 #reference height pressure
LLL.close()

#establish full directories for each run
print('Set directories for each run')

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

print('load ozone data')

#Ozone
RCP45_O3 = dict()
for s in range(len(Sims)):
    RCP45_O3[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_O3[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_O3[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['O3'])
            RCP45_nc.close()

wAVD_O3 = dict()
for s in range(len(Sims)):
    wAVD_O3[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_O3[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_O3[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['O3'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_O3_Panel = pd.Panel.from_dict(RCP45_O3, orient = 'items')
wAVD_O3_Panel = pd.Panel.from_dict(wAVD_O3, orient = 'items')

RCP45_O3_TmSrs_Stacked = RCP45_O3_Panel['001','2006':'2065','01':'12']
RCP45_O3_TmSrs_Stacked = np.asarray(RCP45_O3_TmSrs_Stacked)

wAVD_O3_TmSrs_Stacked = wAVD_O3_Panel['001','2006':'2065','01':'12']
wAVD_O3_TmSrs_Stacked = np.asarray(wAVD_O3_TmSrs_Stacked)

print('unpack and store as a single time series array')

RCP45_Ozone_TmSrs = np.stack((np.hstack(RCP45_O3_TmSrs_Stacked[0:60,0:12])[0:720]),axis=0)
wAVD_Ozone_TmSrs = np.stack((np.hstack(wAVD_O3_TmSrs_Stacked[0:60,0:12])[0:720]),axis=0)


print('Surface Pressure now!')

RCP45_PS = dict()
for s in range(len(Sims)):
    RCP45_PS[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45_PS[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45_PS[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['PS'])
            RCP45_nc.close()

wAVD_PS = dict()
for s in range(len(Sims)):
    wAVD_PS[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD_PS[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD_PS[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['PS'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_PS_Panel = pd.Panel.from_dict(RCP45_PS, orient = 'items')
wAVD_PS_Panel = pd.Panel.from_dict(wAVD_PS, orient = 'items')

RCP45_PS_TmSrs_Stacked = RCP45_PS_Panel['001','2006':'2065','01':'12']
RCP45_PS_TmSrs_Stacked = np.asarray(RCP45_PS_TmSrs_Stacked)

wAVD_PS_TmSrs_Stacked = wAVD_PS_Panel['001','2006':'2065','01':'12']
wAVD_PS_TmSrs_Stacked = np.asarray(wAVD_PS_TmSrs_Stacked)

print('unpack and store as a single time series array')

RCP45_PS_TmSrs = np.zeros([720,96,144])
wAVD_PS_TmSrs = np.zeros([720,96,144])

RCP45_PS_TmSrs = np.stack((np.hstack(RCP45_PS_TmSrs_Stacked[0:60,0:12])[0:720]),axis=0)
wAVD_PS_TmSrs = np.stack((np.hstack(wAVD_PS_TmSrs_Stacked[0:60,0:12])[0:720]),axis=0)

print('Calculate Pressure levels')

#calculate pressure levels
plevs_RCP45 = np.zeros([720,47,96,144])
plevs_wAVD = np.zeros([720,47,96,144])

for i in range(720):
    plevs_RCP45[i,:,:,:] = JF.HybSig2Plev(P0=p0,PS=RCP45_PS_TmSrs[i,:,:],\
        coef_A=hyam,coef_B=hybm)/100 #in hPa

    plevs_wAVD[i,:,:,:] = JF.HybSig2Plev(P0=p0,PS=wAVD_PS_TmSrs[i,:,:],\
        coef_A=hyam,coef_B=hybm)/100 #in hPa

print('Tropopause estimation')

p_tropopause_zonalmean_linear = np.zeros([96])

for y in range(96):
    if y <= len(WACCM4_Lat)/2:
        if y == 0:
            p_tropopause_zonalmean_linear[y] = 300
        elif y == len(WACCM4_Lat)/2:
            p_tropopause_zonalmean_linear[y] = 100
        else:
            p_tropopause_zonalmean_linear[y] = \
            p_tropopause_zonalmean_linear[y-1]\
            -(200/(len(WACCM4_Lat)/2-1))
    else:
        if y == len(WACCM4_Lat)/2:
            p_tropopause_zonalmean_linear[y] = 100
        elif y == len(WACCM4_Lat):
            p_tropopause_zonalmean_linear[y] = 300
        else:
            p_tropopause_zonalmean_linear[y] = \
            p_tropopause_zonalmean_linear[y-1]\
            +(200/(len(WACCM4_Lat)/2-1))

Q = np.tile(p_tropopause_zonalmean_linear.T,(WACCM4_Lon.size,1))
R = np.transpose(Q[:,:,None],(0,1,2))
S = np.tile(R,(1,1,47))
T = np.transpose(S[:,:,:,None],(0,1,2,3))
U = np.tile(T,(1,1,1,720))
V = np.swapaxes(U,3,0)
p_tropopause_WACCM4 = np.swapaxes(V,2,1)

print('Mask Values below the stratosphere')

RCP45_Ozone_TmSrs_Strato = np.zeros([720,47,96,144])
wAVD_Ozone_TmSrs_Strato = np.zeros([720,47,96,144])

    RCP45_Ozone_TmSrs_Strato = RCP45_Ozone_TmSrs*(plevs_RCP45<=p_tropopause_WACCM4)
    wAVD_Ozone_TmSrs_Strato[ = wAVD_Ozone_TmSrs*(plevs_wAVD<=p_tropopause_WACCM4)

RCP45_Ozone_TmSrs_Strato = np.ma.masked_equal(RCP45_Ozone_TmSrs_Strato,0)
wAVD_Ozone_TmSrs_Strato = np.ma.masked_equal(wAVD_Ozone_TmSrs_Strato,0)

RCP45_Ozone_TmSrs_Strato_AA = np.ma.zeros([720,47])
RCP45_Ozone_TmSrs_Strato_GA = np.ma.zeros([720,47])

wAVD_Ozone_TmSrs_Strato_AA = np.ma.zeros([720,47])
wAVD_Ozone_TmSrs_Strato_GA = np.ma.zeros([720,47])

print('area average variables')

for i in range(720):
    RCP45_Ozone_TmSrs_Strato_AA = AA.LatLonavg_Time(\
    RCP45_Ozone_TmSrs_Strato[i,:,80:,:],WACCM4_Lat[80:],WACCM4_Lon)

    RCP45_Ozone_TmSrs_Strato_GA = AA.LatLonavg_Time(\
    RCP45_Ozone_TmSrs_Strato[i,:,:,:],WACCM4_Lat,WACCM4_Lon)

    wAVD_Ozone_TmSrs_Strato_AA = AA.LatLonavg_Time(\
    wAVD_Ozone_TmSrs_Strato[i,:,80:,:],WACCM4_Lat[80:],WACCM4_Lon)

    wAVD_Ozone_TmSrs_Strato_GA = AA.LatLonavg_Time(\
    wAVD_Ozone_TmSrs_Strato[i,:,:,:],WACCM4_Lat,WACCM4_Lon)

RCP45_Ozone_TmSrs_Strato_AA_file = open('RCP45_Ozone_TmSrs_ArcStrato_01.pickle','wb')
pickle.dump(RCP45_Ozone_TmSrs_Strato_AA,RCP45_Ozone_TmSrs_Strato_AA_file)
RCP45_Ozone_TmSrs_Strato_AA_file.close()

RCP45_Ozone_TmSrs_Strato_GA_file = open('RCP45_Ozone_TmSrs_GlobStrato_01.pickle','wb')
pickle.dump(RCP45_Ozone_TmSrs_Strato_GA,RCP45_Ozone_TmSrs_Strato_GA_file)
RCP45_Ozone_TmSrs_Strato_GA_file.close()

wAVD_Ozone_TmSrs_Strato_AA_file = open('wAVD_Ozone_TmSrs_ArcStrato_01.pickle','wb')
pickle.dump(wAVD_Ozone_TmSrs_Strato_AA,wAVD_Ozone_TmSrs_Strato_AA_file)
wAVD_Ozone_TmSrs_Strato_AA_file.close()

wAVD_Ozone_TmSrs_Strato_GA_file = open('wAVD_Ozone_TmSrs_GlobStrato_01.pickle','wb')
pickle.dump(wAVD_Ozone_TmSrs_Strato_GA,wAVD_Ozone_TmSrs_Strato_GA_file)
wAVD_Ozone_TmSrs_Strato_GA_file.close()

print('Fin')
