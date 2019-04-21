#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset
import pandas as pd
import scipy.io as scio
import pickle
import os

Months = ['01','02','03','04','05','06','07','08','09','10','11','12']
Sims = ['001','002','003','004','005']
Sims_alt_name = ['001','002','003','001','002']

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

Years = ['2006','2007','2008','2009','2010','2011','2012','2013','2014','2015']

RCP45 = dict()
for s in range(len(Sims)):
    RCP45[Sims[s]] = dict()
    for m in range(len(Months)):
        RCP45[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_RCP45 = Datadir_RCP45[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            RCP45_nc = Dataset(files_RCP45)
            RCP45[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(RCP45_nc.variables['PS'])
            RCP45_nc.close()

wAVD = dict()
for s in range(len(Sims)):
    wAVD[Sims[s]] = dict()
    for m in range(len(Months)):
        wAVD[Sims[s]][Months[m]] = dict()
        for y in range(len(Years)):
            files_wAVD = Datadir_wAVD[s]+\
            Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
            wAVD_nc = Dataset(files_wAVD)
            wAVD[Sims[s]][Months[m]][Years[y]] = \
            np.squeeze(wAVD_nc.variables['PS'])
            wAVD_nc.close()

#turn all temp data into panels
RCP45_Panel = pd.Panel.from_dict(RCP45, orient = 'items')
wAVD_Panel = pd.Panel.from_dict(wAVD, orient = 'items')


RCP45_Spressure = RCP45_Panel.loc['001':'005','2006':'2015','01':'12']
RCP45_Spressure = np.asarray(RCP45_Spressure)
RCP45_Spressure = np.mean(RCP45_Spressure,axis=(0,1))
RCP45_Spressure_001 = np.stack(RCP45_Spressure[:],axis=0)


wAVD_Spressure = wAVD_Panel.loc['001':'005','2006':'2015','01':'12']
wAVD_Spressure = np.asarray(wAVD_Spressure)
wAVD_Spressure = np.mean(wAVD_Spressure,axis=(0,1))
wAVD_Spressure = np.stack(wAVD_Spressure[:],axis=0)

#pickle arrays
wAVD_Spressure_file = open('wAVD_Base_PS.pickle','wb')
pickle.dump(wAVD_Spressure,wAVD_Spressure_file)
wAVD_Spressure_file.close()

RCP45_Spressure_file = open('RCP45_Base_PS.pickle','wb')
pickle.dump(RCP45_Spressure,RCP45_Spressure_file)
RCP45_Spressure_file.close()
