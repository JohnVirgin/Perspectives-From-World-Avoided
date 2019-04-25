import numpy as np

def LatLonavg_Time(var='var', lat='lat', lon='lon'):
    """Returns the area-averged variable for a certain lat-lon region"""
    y = lat*np.pi/180
    coslat = np.cos(y)
    coslat = np.tile(coslat,(lon.size,1)).T

    #area weighting (weight by cos(lat) because area of grid boxes get smaller closer to the pole)
    lat_tmp_cos = np.ma.zeros([var[:,0,0].size, lat.size, lon.size])
    coslat_tmp = np.ma.zeros([var[:,0,0].size, lat.size, lon.size])
    for i in range(lat.size):
        for j in range(lon.size):
            for k in range(var[:,0,0].size):
                try:
                    if var.mask[k,i,j] == False:
                        lat_tmp_cos[k,i,j] = var[k,i,j] * coslat[i,j]
                        coslat_tmp[k,i,j] = coslat[i,j]
                    else:
                        lat_tmp_cos[k,i,j] = float('nan')
                        coslat_tmp[k,i,j] = float('nan')
                except:
                    lat_tmp_cos[k,i,j] = var[k,i,j] * coslat[i,j]
                    coslat_tmp[k,i,j] = coslat[i,j]
                else:
                    pass

    #sum over all gridpoints
    cos_sum = np.nansum(coslat_tmp, axis=(1,2))

    lat_sum = np.nansum(lat_tmp_cos, axis=(1,2))

    #normalize by area
    LatLonavg = lat_sum / cos_sum

    return LatLonavg
