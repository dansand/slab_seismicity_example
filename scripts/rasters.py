#Credits: 

#To do. Would be good to be able to read geotiffs using gdal/rasterio

import sys                                                                      
sys.path.insert(1, './scripts/')
from geodesy import *


from scipy import ndimage as nd

try:
    from netCDF4 import Dataset as netcdf
except ImportError:
    from scipy.io import netcdf_file as netcdf
    print('Warning: NetCDF4 grids not supported ("netCDF4" Python module not found). '
          'Falling back to NetCDF3, rasters may fail to load.')

    
    
import scipy.interpolate as spi
import numpy as np
import gdal

def sample_grid_using_scipy(x,y,grdfile, zvar= 'z'):
    
    data=netcdf(grdfile,'r')
    try:
        lon = np.copy(data.variables['x'][:])
        lat = np.copy(data.variables['y'][:])
    except:
        lon = np.copy(data.variables['lon'][:])
        lat = np.copy(data.variables['lat'][:])
    
    Zg = data.variables[zvar][:]
    
    test = fill_ndimage(Zg)
    
    lut=spi.RectBivariateSpline(lon,lat,test.T)
    result = []
    for xi,yi in zip(x,y):
        result.append(lut(xi, yi)[0][0])
            
    return np.array(result)



def fill_ndimage(data,invalid=None):
    """Replace the value of invalid 'data' cells (indicated by 'invalid')
    by the value of the nearest valid data cell
    Parameters
    ----------
    data: numpy array of any dimension
    invalid: a binary array of same shape as 'data'. True cells set where data
    value should be replaced.
    If None (default), use: invalid = np.isnan(data)
    Returns
    -------
    Return a filled array.
    Credits
    -------
    http://stackoverflow.com/a/9262129
    """
    if invalid is None: invalid = np.isnan(data)
    ind = nd.distance_transform_edt(invalid, return_distances=False, return_indices=True)
    return data[tuple(ind)]



def preprocess_grd_file(grdfile, zvar= 'z'):
    """
    
    """
    
    data_dict = {}
    
    data=netcdf(grdfile,'r')
    
    try:
        lon = np.copy(data.variables['x'][:])
        lat = np.copy(data.variables['y'][:])
    except:
        lon = np.copy(data.variables['lon'][:])
        lat = np.copy(data.variables['lat'][:])
        
    lon = wrap_to_360(lon)
    
    zarray = data.variables[zvar][:]
    
    #next bit creates a binary array that we can use to query where the initial mask was
    #it would be simpler to use the mask directly (np.invert(zarray.mask).astype('float')), 
    #but this strategy failed in cases where there were no values to mask 
    
    mask_array = np.zeros(zarray.shape)
    mask_array[zarray.mask] = 0
    mask_array[~zarray.mask] = 1
    
    #This fills NaN values using nearest neighbour interpolation
    filled_array = fill_ndimage(zarray)
    
    #fill out the dictionary
    data_dict['lon'] = lon
    data_dict['lat'] = lat
    data_dict['filled_array'] = filled_array
    
    #data_dict['mask_array'] = np.invert(zarray.mask).astype('float')
    data_dict['mask_array'] = mask_array
    
    data_dict['lut_data'] = spi.RectBivariateSpline(lon,lat, filled_array.T)
    data_dict['lut_mask'] = spi.RectBivariateSpline(lon,lat, mask_array.T)
    
    return data_dict 


def interp_raster_at_points(data_dict, points, thresh = 0.5, pygplates_order = False):
    
    """
    threshhold us used to determine where valid/NaN data were stored in the inital array
    """
    if pygplates_order is False:
        x = wrap_to_360(points[:,0])
        y = points[:,1]
    else:
        x = wrap_to_360(points[:,1])
        y = points[:,0]
    
    lon = wrap_to_360(data_dict['lon'])
    lat = data_dict['lat']

    
    #filled_array = data_dict['filled_array']
    #mask_array = data_dict['mask_array']

    #lut=spi.RectBivariateSpline(lon,lat, filled_array.T)
    result = data_dict['lut_data'](x, y, grid=False)
    
    #lut_=spi.RectBivariateSpline(lon,lat, mask_array.T)
    values_to_mask = data_dict['lut_mask'](x, y, grid=False)
    
    mask = values_to_mask > thresh
    
    return result, mask



def raster_extent_gdal(filename):
    
    slab_raster = gdal.Open(filename, gdal.GA_ReadOnly)

    slab_data = -1.*slab_raster.ReadAsArray() #positive depths kms
    geo_transform = slab_raster.GetGeoTransform()

    origin_x = geo_transform[0]
    origin_y = geo_transform[3]
    pixel_width = geo_transform[1]
    pixel_height = geo_transform[5]

    extent_lonlat = (
        origin_x, 
        origin_x + (pixel_width * slab_raster.RasterXSize),
        origin_y + (pixel_height * slab_raster.RasterYSize),
        origin_y
    )
    
    return extent_lonlat