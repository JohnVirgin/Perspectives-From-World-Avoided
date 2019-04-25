## This file contains a suite of python functions designed to modify, analyze,
  ## and manipulate GCM-oriented data on Lat/Lon/Lev grids.
#Written by John Virgin
#Master File Creation Date: December 23rd, 2018
#Most Recent Revision/Update: February 6th, 2019

import numpy as np
import scipy.stats as stats

#Creation Date: December 13th,2018
def PlevThck(PS='PS', plevs='plevs',p_top='p_top'):

    ''' This function uses GCM surface pressure and a set of prescribed pressures
    throughout an atmosphere to calculute the layer thickness for all levels -
    all values must be in the same units (hpa or pa). '''

    if min(plevs) == plevs[0]:
        plevs = plevs
    else:
        plevs = plevs[::-1]

    Upper = 0
    Mid = 1
    Lower = 2

    dp = np.empty([plevs.size])
    dp[:] = np.nan

    def lower_boundary():
        for i,value in enumerate(plevs):
            if value > PS:
                    return i

        return plevs.size

    for i in range(lower_boundary()):

        if i == 0:
            dp[i] = -(p_top-(((plevs[1]+plevs[2])/2)-plevs[0]))

        elif Lower < lower_boundary():
            dp[i] = ((plevs[Mid]+plevs[Lower])/2)-((plevs[Upper]+plevs[Mid])/2)

            Upper += 1
            Mid += 1
            Lower += 1

        else:
            dp[i] = PS-(((plevs[lower_boundary()-1]+plevs[lower_boundary()-2])/2))


    return dp

#Creation Date: December 12th, 2018
def HybSig2Plev(P0 = 'P0', PS = 'PS', coef_A = 'coef_A', coef_B = 'coef_B'):

    ''' This function uses a reference pressure, the surface pressure at each lat/lon point,
    and either mid or interface points of hybrid level coefficients to calculate the
    pressure levels of the variables corresponding to your Surface pressure variable.
    Note that all of these variables be the same the units (either pascals/millibars
    or hectopascals).'''

    hyam = coef_A
    hybm = coef_B

    Plevs = np.zeros([hyam.size,PS[:,0].size,PS[0,:].size])

    for z in range(hyam.size): #loop over levels
        for x in range(PS[:,0].size): #loop over latitudes
            for y in range(PS[0,:].size): #loop over longitudes
                Plevs[z,x,y] = hyam[z]*P0+(hybm[z]*PS[x,y])

    return Plevs  #return pressure levels with the four dimensions

#Creation Date: December 23rd, 2018
def HybSig2PlevThck(P0 = 'P0', PS = 'PS', coef_A = 'coef_A', coef_B = 'coef_B'):

    ''' This function uses a reference pressure, the surface pressure at each lat/lon point,
    and either mid or interface points of hybrid level coefficients to calculate the
    pressure levels of the variables corresponding to your Surface pressure variable.
    Note that all of these variables be the same the units (either pascals/millibars
    or hectopascals). The extra addition to this function is the thickness of
    each level that's calculate, which is useful for pressure weighted variables.'''

    hyam = coef_A
    hybm = coef_B

    Dim1 = hyam.size
    Dim2 = PS[:,0].size
    Dim3 = PS[0,:].size

    Plevs = np.zeros([Dim1,Dim2,Dim3])

    for z in range(Dim1): #loop over levels
        for x in range(Dim2): #loop over latitudes
            for y in range(Dim3): #loop over longitudes
                Plevs[z,x,y] = hyam[z]*P0+(hybm[z]*PS[x,y])


    plevs_thck = np.zeros([Dim1-1,Dim2,Dim3])

    Upper=0
    Lower=1

    for z in range(Dim1-1):
        for x in range(Dim2):
            for y in range(Dim3):
                plevs_thck[z,x,y] = abs(Plevs[Upper,x,y]-Plevs[Lower,x,y])
        Upper+=1
        Lower+=1

    return plevs_thck


#Creation Date: June 22nd, 2018
def dOceanHeat_Wm2(SST = 'SST', dz = 'dz', Timestep = str):


            ''' A few notes before we begin. The SST's must be in K, not degrees
            Celcius. They must also be a function of all possible parametres
            (time,depth,lat,lon). Model depth must be in m^2, and must also be
            masked to discount the fill value. Lastly, SST's must be masked
            where land is on the grid.'''

            if Timestep == 'Monthly':
                #number of seconds per month
                dt = 2592000
            elif Timestep == 'Daily':
                #number of seconds per day
                dt = 86400
            elif Timestep == 'Yearly':
                #number of seconds per year
                dt = 3.154e+7
            elif Timestep == 'Hourly':
                #number of seconds per hour
                dt = 3600.45662112
            else:
                raise ValueError(\
                "Time step must be Monthly, Daily, Yearly, or Hourly")

            #specific heat capacity of sea water
            cP = 3985#J kg^-1 K^-1

            #reference density of sea water
            rHO = 1025#kg m^-3

            SST_K = SST+273.15 #conversion to Kelvin

            #calculate ocean heat content in J m^-2
            ohc_Jm = np.zeros(\
            [SST_K[:,0,0,0].size,SST_K[0,0,:,0].size,SST_K[0,0,0,:].size])

            for i in range(SST_K[:,0,0,0].size):
                ohc_Jm[i,:,:] = cP*rHO*\
                (np.ma.sum(np.ma.multiply(SST_K[i,:,:,:],dz),axis=0))

            #Calculate time rate of change of vertically integrated ohc in
            #W m^-2
            ohc_Wm = np.zeros([ohc_Jm[:,0,0].size,ohc_Jm[0,:,0].size,\
            ohc_Jm[0,0,:].size])

            for i in range(ohc_Wm[:,0,0].size):
                for j in range(ohc_Wm[0,:,0].size):
                    for k in range(ohc_Wm[0,0,:].size):
                        if i == 0:
                            ohc_Wm[i,j,k] = \
                            (ohc_Jm[i+1,j,k]-ohc_Jm[i,j,k])/dt
                        elif i == ohc_Jm[:,0,0].size-1:
                            ohc_Wm[i,j,k] = \
                            (ohc_Jm[i,j,k]-ohc_Jm[i-1,j,k])/dt
                        else:
                            ohc_Wm[i,j,k] = \
                            (((ohc_Jm[i+1,j,k]+ohc_Jm[i,j,k])/2)\
                            -((ohc_Jm[i,j,k]+ohc_Jm[i-1,j,k])/2))/dt

            #Return ohc in W m^-2
            return ohc_Wm

#Creation Date: June 22nd, 2018
def dOceanHeat_Jm2(SST = 'SST',dz = 'dz',OCN_GridCell_Area = 'OCN_GridCell_Area'):


        ''' A few notes before we begin. The SST's must be in K, not degrees
        Celcius. They must also be a function of all possible parametres
        (time,depth,lat,lon). Model depth must be in m^2, and must also be
        masked to discount the fill value. Lastly, SST's must be masked
        where land is on the grid.'''

        #specific heat capacity of sea water
        cP = 3985#J kg^-1 K^-1

        #reference density of sea water
        rHO = 1025#kg m^-3

        SST_K = SST+273.15 #conversion to Kelvin

        #calculate ocean heat content in J m^-2
        ohc_Jm = np.zeros(\
        [SST_K[:,0,0,0].size,SST_K[0,0,:,0].size,SST_K[0,0,0,:].size])

        for i in range(SST_K[:,0,0,0].size):
            ohc_Jm[i,:,:] = cP*rHO*\
            (np.ma.sum(np.ma.multiply(SST_K[i,:,:,:],dz),axis=0))

        OHC_J = np.zeros(\
        [ohc_Jm[:,0,0].size,ohc_Jm[0,:,0].size,ohc_Jm[0,0,:].size])
        for i in range(ohc_Jm[:,0,0].size):
            OHC_J[i,:,:] = np.ma.multiply(ohc_Jm[i,:,:],OCN_GridCell_Area)

        #Return OHC in Jules
        return OHC_J

#Creation Date: January 2nd, 2019
def Moving_Avg(var='var', window = int):

    '''This function applies a moving average to smooth a time series of data with a particular window size
    Important to note that the edges of the particular time series will be cut off using this function,
    and the series will be shorter in length'''

    Moving_Avg = np.convolve(var,np.ones(window,)/window, mode = 'valid')

    return Moving_Avg

#Creation Date: January 3rd, 2019
def Calc_SatSpec_Hum(Ta = 'Ta',P = 'P'):

    '''This function calculates saturation specific humidity at each grid point
    using air temperatures and corresponding pressure (The arrays must have the same dimensions).
    The equations are sourced from Buck (1981), and this specific function is sourced from
    Dr. Angeline Pendergrass' Github page, which was originally written for matlab'''

    Ta = Ta-273.15 #convert temperature to degress celsius
    MWratio = 0.622 #molecular weight ratio of water vapour to dry air

    #equation 1: saturation vapour pressure with respect to water
    #the hard coded values here are coefficients for general usage
    Sat_P_water = (1.0007+(3.46e-6*P))*6.1121*(np.exp((17.502*Ta)/(240.97+Ta)))

    # saturation mixing ratio with respect to liquid water (g of water/kg of dry air)
    ws_liquid = MWratio*Sat_P_water/(P-Sat_P_water)

    #equation 2: saturation vapour pressure with respect to ice
    #the hard coded values here are coefficients for general usage
    Sat_P_ice = (1.0003+(4.18e-6*P))*6.1115*(np.exp((22.452*Ta)/(272.55+Ta)))

    # saturation mixing ratio with respect to ice (g of water/kg of dry air)
    ws_ice = MWratio*Sat_P_ice/(P-Sat_P_ice)

    ws = ws_liquid
    ws[Ta<0] = ws_ice[Ta<0]

    #saturation specific humidity (g/kg)
    Qs = ws/(1+ws)

    return Qs

#Creation date: January 8th, 2019
def Mean_Confidence_Interval(Data='Data', Confidence = 'Confidence'):

    ''' Uses scipy stats library to calculate a confidence interval for the mean of
    a set of data. Returns four values in the order of Mean, lower bound interval,
    upper bound interval, and interval value itself.'''

    a = 1.0 * np.array(Data)
    n = a.size
    Mean, Std_Er = np.mean(a), stats.sem(a)
    h = Std_Er*stats.t.ppf((1 + Confidence) / 2., n-1)
    return Mean, Mean-h, Mean+h, h
