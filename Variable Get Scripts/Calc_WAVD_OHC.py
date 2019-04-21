#!/usr/bin/env python

## Import Packages
import numpy as np
import pandas as pd
import matplotlib as mpl
import scipy.stats as stats
import scipy.io as scio
import CalcOceanHeatContent as coh
import os
import pickle
from netCDF4 import Dataset

# Read in all necessary dimensions and external variables
LLL = Dataset('CESM-WACCM Ocean Grid Cell Thickness/thkcello_CESM1-WACCM_interp.nc')
WACCM_Lat = np.squeeze(LLL.variables['lat'])
WACCM_Arctic_Lat = WACCM_Lat[150:]
WACCM_Lon = np.squeeze(LLL.variables['lon'])
WACCM_Depth = np.squeeze(LLL.variables['lev']) #in metres
dz = np.squeeze(LLL.variables['thkcello']) #in metres
dz = np.ma.masked_equal(dz,1e+20)
LLL2 = Dataset('CESM-WACCM Ocean Grid Cell Thickness/areacello_CESM1-WACCM_interp.nc')
Ocn_GridCell_Area = np.squeeze(LLL2.variables['areacello']) #in metres squared
Ocn_GridCell_Area = np.ma.masked_equal(Ocn_GridCell_Area,1e+20)
Ocn_GridCell_Area_Arctic = Ocn_GridCell_Area[150:,:]

Years = list(map(str, (range(2006,2066))))
Months = ['01','02','03','04','05','06','07','08','09','10','11','12']

Datadir_WAVD = '/export/data/waccmdata/b40.rcp4_5.2deg.wcm.WAVD.KLS.001/ocn/'

WAVD_SST = dict()
for m in range(len(Months)):
    WAVD_SST[Months[m]] = dict()
    for y in range(len(Years)):
        files_WAVD = Datadir_WAVD+'b40.rcp4_5.2deg.wcm.WAVD.KLS.001.'\
        +Years[y]+'.'+Months[m]+'.nc'
        WAVD_nc = Dataset(files_WAVD)
        WAVD_SST[Months[m]][Years[y]] = \
        np.squeeze(WAVD_nc.variables['TEMP'])
        WAVD_nc.close()


WAVD_SST_DF = pd.DataFrame.from_dict(WAVD_SST, orient = 'columns')

WAVD_SST_Array = WAVD_SST_DF.loc['2006':'2065','01':'12']
WAVD_SST_Array = np.asarray(WAVD_SST_Array)
WAVD_SST_TmSrs = np.hstack((WAVD_SST_Array[0:61,0:12]))
WAVD_SST_TmSrs = np.stack((WAVD_SST_TmSrs[0:732]), axis=0)
WAVD_SST_TmSrs = np.ma.masked_greater(WAVD_SST_TmSrs,100)

WAVD_dOHC_TmSrs_Wm = coh.dOHC_Wm(SST = WAVD_SST_TmSrs,dz = dz,Timestep = "Monthly")

WAVD_dOHC_TmSrs_Wm_file = open("WAVD_OHC_TimeSeries_Wm_EM4.pickle",'wb')
pickle.dump(WAVD_dOHC_TmSrs_Wm,WAVD_dOHC_TmSrs_Wm_file)
WAVD_dOHC_TmSrs_Wm_file.close()
