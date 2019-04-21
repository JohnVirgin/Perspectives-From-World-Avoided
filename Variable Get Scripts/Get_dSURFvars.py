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

Years = list(map(str, (range(2006,2066))))
Months = ['01','02','03','04','05','06','07','08','09','10','11','12']
Vars = ['FLNS','FLNSC','FSNS','FSNSC','FLUT','FLUTC','FSNT','FSNTC','LHFLX','SHFLX','CLDLOW','CLDHGH','CLDMED']
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
Datadir_WAVD = [datadir_wAVD_001,datadir_wAVD_002,datadir_wAVD_003,\
                datadir_wAVD_004,datadir_wAVD_005]

print('Reading in data')

RCP45_Vars = dict()
WAVD_Vars = dict()
for v in range(len(Vars)):
    RCP45_Vars[Vars[v]] = dict()
    WAVD_Vars[Vars[v]] = dict()
    for s in range(len(Sims)):
        RCP45_Vars[Vars[v]][Sims[s]] = dict()
        WAVD_Vars[Vars[v]][Sims[s]] = dict()
        for m in range(len(Months)):
            RCP45_Vars[Vars[v]][Sims[s]][Months[m]] = dict()
            WAVD_Vars[Vars[v]][Sims[s]][Months[m]] = dict()
            for y in range(len(Years)):

                files_RCP45 = Datadir_RCP45[s]+Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
                RCP45_nc = Dataset(files_RCP45)
                RCP45_Vars[Vars[v]][Sims[s]][Months[m]][Years[y]] = np.squeeze(RCP45_nc.variables[Vars[v]])
                RCP45_nc.close()

                files_WAVD = Datadir_WAVD[s]+Sims_alt_name[s]+'.'+Years[y]+'.'+Months[m]+'.nc'
                WAVD_nc = Dataset(files_WAVD)
                WAVD_Vars[Vars[v]][Sims[s]][Months[m]][Years[y]] = np.squeeze(WAVD_nc.variables[Vars[v]])
                WAVD_nc.close()

print('Done')

print('Modifying')

RCP45_Vars_Panel = dict()
WAVD_Vars_Panel = dict()

RCP45_Vars_1D = dict()
WAVD_Vars_1D = dict()

RCP45_Vars_2D = dict()
WAVD_Vars_2D = dict()

RCP45_Vars_1Dstack = dict()
WAVD_Vars_1Dstack = dict()

RCP45_Vars_2Dstack = dict()
WAVD_Vars_2Dstack = dict()

RCP45_dVars = dict()
WAVD_dVars = dict()

for keys,values in RCP45_Vars.iteritems():
    RCP45_Vars_Panel[keys] = pd.Panel.from_dict(RCP45_Vars[keys], orient = 'items')
    WAVD_Vars_Panel[keys] = pd.Panel.from_dict(WAVD_Vars[keys], orient = 'items')

    RCP45_Vars_1D[keys] = np.asarray(RCP45_Vars_Panel[keys]['001':'005','2006':'2015','01':'12'])
    RCP45_Vars_1D[keys] = np.mean(RCP45_Vars_1D[keys],axis=1)
    RCP45_Vars_2D[keys] = np.asarray(RCP45_Vars_Panel[keys]['001':'005','2056':'2065','01':'12'])
    RCP45_Vars_2D[keys] = np.mean(RCP45_Vars_2D[keys],axis=1)

    WAVD_Vars_1D[keys] = np.asarray(WAVD_Vars_Panel[keys]['001':'005','2006':'2015','01':'12'])
    WAVD_Vars_1D[keys] = np.mean(WAVD_Vars_1D[keys],axis=1)
    WAVD_Vars_2D[keys] = np.asarray(WAVD_Vars_Panel[keys]['001':'005','2056':'2065','01':'12'])
    WAVD_Vars_2D[keys] = np.mean(WAVD_Vars_2D[keys],axis=1)

    RCP45_Vars_1Dstack[keys] = np.zeros([5,12,96,144])
    RCP45_Vars_2Dstack[keys] = np.zeros([5,12,96,144])

    WAVD_Vars_1Dstack[keys] = np.zeros([5,12,96,144])
    WAVD_Vars_2Dstack[keys] = np.zeros([5,12,96,144])

    for i in range(len(Sims)):
        RCP45_Vars_1Dstack[keys][i,:,:,:] = np.stack(RCP45_Vars_1D[keys][i,:])
        RCP45_Vars_2Dstack[keys][i,:,:,:] = np.stack(RCP45_Vars_2D[keys][i,:])

        WAVD_Vars_1Dstack[keys][i,:,:,:] = np.stack(WAVD_Vars_1D[keys][i,:])
        WAVD_Vars_2Dstack[keys][i,:,:,:] = np.stack(WAVD_Vars_2D[keys][i,:])


    RCP45_dVars[keys] = RCP45_Vars_2Dstack[keys]-RCP45_Vars_1Dstack[keys]
    WAVD_dVars[keys] = WAVD_Vars_2Dstack[keys]-WAVD_Vars_1Dstack[keys]

print('Done')

print('Saving')

RCP45_dVars_file = open("RCP45_dFlux_Vars.pickle",'wb')
pickle.dump(RCP45_dVars,RCP45_dVars_file)
RCP45_dVars_file.close()

WAVD_dVars_file = open("WAVD_dFlux_Vars.pickle",'wb')
pickle.dump(WAVD_dVars,WAVD_dVars_file)
WAVD_dVars_file.close()

print('Done')
