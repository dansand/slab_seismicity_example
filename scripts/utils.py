"""
A collection of functions to help with seismicity analysis in slabs


:license: GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)
"""

  

import sys                                                                      
sys.path.insert(1, '/usr/lib/pygplates/revision18/')
sys.path.insert(1, './scripts/')
from geodesy import *


import os
import shutil
import numpy as np
import math
from scipy.ndimage.filters import gaussian_filter
from scipy.interpolate import CubicSpline
import pygplates
import time
import obspy

import cartopy.crs as ccrs






#i don't know of a more efficeint way to get all the catalog lons / lats, 
#These need more testing. The try/except loop is for different structure in CMT and ISC-focal mech xmls.
def get_lon(cat):
    lons = []
    for ev in cat:
        try:
            lon = ev.preferred_origin()['longitude']
        except:
            lon = ev.origins[0]['longitude']
        lons.append(lon)
    return np.array(lons)
            
        
def get_lat(cat):
    lats = []
    for ev in cat:
        try:
            lat = ev.preferred_origin()['latitude']
        except:
            lat = ev.origins[0]['latitude']
            
        lats.append(lat)
    return np.array(lats)
            
        
def get_depth(cat):
    depths = []
    for ev in cat:
        
        try:
            depth = ev.preferred_origin()['depth']
        except:
            depth = ev.origins[0]['depth']
        
        if depth is None:
            depth = 0.
        depths.append(depth)    
    
    return np.array(depths)
            
        
def get_mag(cat):
    mags = []
    for ev in cat:
        try:
            m = np.max([m['mag'] for m in ev.magnitudes]) #average of provided mags
        except:
            m = 5.0
        mags.append(m)
    return np.array(mags)


def clean_trench_points(subzone, gauss_filter_trench = 0, 
                        north=None,south=None, east=None, west=None, profile_spacing_km=10., earth_radius=6371.009):

      
    '''
    Takes a  pygplates feature collection, containing a subduction trench or multiple trench segments, 
    returns tesselated data
    '''
    
      
    profile_rads = profile_spacing_km/earth_radius

    #There may be multiple features (segments) in the shapefile
    points_list= []
    for feature in subzone:
        subzone = feature.get_geometry()
        subzone  = subzone.to_tessellated(profile_rads) #resample to get hight precision on domain
        points_list.append(subzone.to_lat_lon_array())

   
    #trench_points = np.row_stack((points_list[0], points_list[1]))
    
    #for varying length lists. This does raise an error. Must be a better way. 
    trench_points = np.row_stack((p for p in points_list))


    #apply the gaussian smoothing before choosing the sub region of the analyis
    trench_points = np.column_stack((gaussian_filter(trench_points[:,0], gauss_filter_trench), 
                            gaussian_filter(trench_points[:,1], gauss_filter_trench)))

    bool_array = np.bool_(np.ones(len(trench_points)))
    maskN, maskS, maskE, maskW = bool_array, bool_array, bool_array, bool_array
    if north:
        maskN = trench_points[:,0] <= north
    if south:
        maskS = trench_points[:,0] >= south
    if east:
        maskE = wrap_to_360(trench_points[:,1]) <= wrap_to_360(east)
    if west:
        maskW = wrap_to_360(trench_points[:,1]) >= wrap_to_360(west)
        
    mask = np.all(np.column_stack((maskN, maskS, maskE, maskW)), axis = 1)
        

    subzone = pygplates.PolylineOnSphere(trench_points[mask,:])
    trench_points = subzone.to_lat_lon_array()
    
    return subzone, trench_points




def get_trench_azimuths(trench_points, plate_vel_azimuths= None):
    
    '''
    Returns azimuth of trench based on the convention that subductuion direction is on the RHS of the trench, 
    when facing along the points. Also return the trench_points is the correct order.
    
    '''
    
    trench_azimuths = []
    trench_point_dists = []

    for i in range(1, len(trench_points) - 1):
        a = calculate_initial_compass_bearing(trench_points[::1][i-1], trench_points[::1][i + 1])
        trench_azimuths.append(a)

    trench_azimuths = np.array(trench_azimuths)
    trench_azimuths = np.pad(trench_azimuths, 1, 'edge')


    #Make sure the subductuion direction is on the RHS of the trench, facing along the points    
    #if we have no information about plate velcoity, skip.
    if plate_vel_azimuths is not None:
        
        #this is perhaps too crude: assumung mean is good enough
        if (plate_vel_azimuths.mean()  - trench_azimuths.mean()) > 180.:

            trench_points = trench_points[::-1]
            trench_azimuths = []
            trench_point_dists = []

            for i in range(1, len(trench_points) - 1):
                a = calculate_initial_compass_bearing(trench_points[::1][i-1], trench_points[::1][i + 1])
                trench_azimuths.append(a)

            trench_azimuths = np.array(trench_azimuths)
            trench_azimuths = np.pad(trench_azimuths, 1, 'edge')
            
            trench_points = trench_points[::1]
            
            
    return trench_azimuths, trench_points



def create_profile_lines(trench_points, azimuths,forearc_dist_km, ocean_distance_fac,
                  profile_inc_km=2, earth_radius=6371.009):
    
    '''
    Create a set of great circles, tesselated to profile_inc_km, based in provided paramters
    '''

    slab_normal_pts = []
    ocean_normal_pts = []

    for i in range(len(trench_points)):
        
        if isinstance(azimuths, int) or isinstance(azimuths, float):
            line_azimuth = azimuths

        elif len(azimuths) == len(trench_points):
            line_azimuth = azimuths[i]
        
        else:
            line_azimuth = azimuths[0]

        #create a set of points in the subdcution direction
        ptslab = great_circle(distance=forearc_dist_km*1e3, 
                     azimuth = line_azimuth ,
                     latitude=trench_points[i][0], 
                     longitude=trench_points[i][1], 
                     rmajor=earth_radius*1e3, 
                     rminor=earth_radius*1e3)

        slab_normal_pts.append((ptslab['latitude'], ptslab['longitude']))


        #create a set of points in the seaward direction
        ptocean = great_circle(distance=forearc_dist_km*ocean_distance_fac*1e3, 
                     azimuth = line_azimuth -180.,
                     latitude=trench_points[i][0], 
                     longitude=trench_points[i][1], 
                     rmajor=earth_radius*1e3, 
                     rminor=earth_radius*1e3)

        ocean_normal_pts.append((ptocean['latitude'], ptocean['longitude']))
    
    
    
    slab_normal_pts = np.array(slab_normal_pts)
    ocean_normal_pts = np.array(ocean_normal_pts)
        
    print(slab_normal_pts.shape)
    
    min_line_length = 1e6
    
    lineList = []
    for i in range(len(trench_points)):
        gca = pygplates.GreatCircleArc(ocean_normal_pts[i], slab_normal_pts[i])
        line =  gca.to_tessellated(profile_inc_km/earth_radius)
        lineList.append(line)
        
        #Enforce llines to be the same number of points
        if len(line) < min_line_length:
            min_line_length = len(line)
    
    

    profile_lines = [pygplates.PolylineOnSphere(l[:min_line_length]) for l in lineList]
        
    #Generate the profile distances
    profile_distances = []
    for i in range(0, len(profile_lines)):

        p0  = pygplates.PointOnSphere(ocean_normal_pts[i])
        p1  = pygplates.PointOnSphere(trench_points[i])
        start_to_trench_kms = pygplates.GeometryOnSphere.distance(p0 , p1)*earth_radius
        indv_dists = []
        for j in range(0,len(profile_lines[i])):
            point = pygplates.PointOnSphere( profile_lines[i][j])
            indv_dists.append(pygplates.GeometryOnSphere.distance(p0 , point)*earth_radius -  start_to_trench_kms )
        profile_distances.append(np.array(indv_dists ))
            
    return slab_normal_pts, ocean_normal_pts, profile_lines, profile_distances


def indexes_points_in_polygons(lons, lats, polygon):
    
    '''
    Get the indexes of the points (lons, lats) that lie in the provided pologon (pygplates.PolygonOnSphere)
    '''
    
    
    #Here I'm adding a small random component, this is a hack that will help with mapping from.. 
    #pygplates points back to the original arrays
    lons = lons + np.sign(lons)*(np.random.rand(len(lons))*1e-6) 
    lats = lats + np.sign(lons)*(np.random.rand(len(lons))*1e-6)

    gpoints = pygplates.MultiPointOnSphere(np.column_stack((lats, lons)))

    points_in_poly = []
    polygon.partition(gpoints, partitioned_geometries_inside=points_in_poly)
    if len(points_in_poly) > 0:
        #The way to do this is to define a mapping back to the original arrays
        indexes = []
        lats_ = points_in_poly[0].to_lat_lon_array()[:,0]
        lons_ = points_in_poly[0].to_lat_lon_array()[:,1]
        for k in range(len(lats_)):
            test = np.abs(lats - lats_[k]).argmin()
            indexes.append(test) 
    else:
        indexes = []
        
    return indexes  



def map_points_to_profile_line(lons, lats, profile_lines):
    
    '''
    Find the closest point on a set of PolyLines (profile_lines), given a set of points, provides as lons and lats
    Returned two lists of indexes, corresponding to the closest line in profile_lines, and the closes pointr in that line.
    
    '''
    
    gpoints = pygplates.MultiPointOnSphere(np.column_stack((lats, lons)))

    point_profile_line_map = []
    profile_line_index_map = []

    for pt in gpoints:
        mindist = 1e6
        minindx1 = -1
        minindx2 = -1
        for i in range(len(profile_lines)):


            dist_, ip, il= pygplates.GeometryOnSphere.distance(pt , profile_lines[i],  return_closest_indices=True)
            if dist_ < mindist:
                mindist = dist_
                minindx1 = i
                minindx2 = il
        point_profile_line_map.append(minindx1)
        profile_line_index_map.append(minindx2)
        
    return  point_profile_line_map, profile_line_index_map



def point_distances_by_profile_line(profile_distances, point_profile_line_map, profile_line_index_map):
    
    '''
    Determinine the distance from a point, by using the known distance in the nearest line
    Point coordinates are encoded with the mapping to the profile lines:
    point_profile_line_map, profile_line_index_map, which are generated with map_points_to_profile_line()
    '''

    eq_distances = []
    for i in range(len(point_profile_line_map)):
        m1 = point_profile_line_map[i]
        m2 = profile_line_index_map[i]
        trench_dist = profile_distances[m1][m2]
        eq_distances.append(trench_dist)

    point_distances = np.array(eq_distances)
    return point_distances



def catalog_time_comparison(times1, times2, tol_seconds = 10):
    
    """
    times1/2 can be either obspy.core.utcdatetime.UTCDateTime
    or an array of seconds
    """
    
    #using the the datetime objects were really slow. 
    #Converting to seconds to run the comparison
    #Speedup is something like (1e3)!!!

    if isinstance(times1[0], obspy.core.utcdatetime.UTCDateTime):
        times1 = np.array([time.mktime(t.timetuple()) for t in np.array(times1)])
        
    if isinstance(times2[0], obspy.core.utcdatetime.UTCDateTime):
        times2 = np.array([time.mktime(t.timetuple()) for t in np.array(times2)])
    


    mask1_ = []
    mask2_ = []
    ix = 0
    for t in times1:
        dt = times2- t                           #Produces an array

        if np.any(np.abs(dt) < tol_seconds):
            mask1_.append(ix)                    #index of a shared event in first set
            mask2_.append(np.argmin(np.abs(dt))) #index of a shared event in first set

        ix+=1
    return mask1_, mask2_


#distance from trench
def nearest(points, pigPlatesObject = None):
    """
    Use the pygplates functionality to get distance between points and object
    """
    
    for g  in points:
        yield pygplates.GeometryOnSphere.distance(pigPlatesObject, g)
        
        
        
def gca_radialdepth_to_cartesian(distance, depth, earth_radius_km= 6371.009):
    radians_ = (distance)/earth_radius_km
    r_ = earth_radius_km - depth 
    cart_distance = r_*np.sin(radians_)
    cart_depth = earth_radius_km - (r_*np.cos(radians_))
    
    return cart_distance, cart_depth


def curve_distance(pts):
    dxs = np.diff(pts[:,0])
    dys = np.diff(pts[:,1])
    ds = np.sqrt(dxs**2 + dys**2)
    
    #pad the first point
    return np.cumsum(np.append([0.], ds)) 


def interp_shift_output(distance, depth, shift=0., coord_gauss = 2, dip_angle_gauss = 2):
    
    
    nanmask = np.invert(np.logical_or(np.isnan(distance), np.isnan(depth)))
    
    #generate a spline that mainly serves to span the missing data between the 
    #slab model and the bathymetry data
    
    fy = CubicSpline(gaussian_filter(distance[nanmask], coord_gauss), 
                 gaussian_filter(depth[nanmask] , coord_gauss))

    #Now resample the points at even intervals in the horizontal/tangential coordinate
    dx_mean = np.diff(distance[nanmask]).mean()
    even_distance = np.arange(distance[nanmask][0], distance[nanmask][-1], dx_mean)
    even_depth = fy(even_distance)
    
    
    #Recreate the spline to estimate derivatives
    fy = CubicSpline(even_distance,even_depth)
    y1_ = fy.derivative(1)(even_distance)
    
    #determine the slab dip. 
    #Note that a lot of filtering may be requited here to see a smooth profile
    dip_angles = gaussian_filter(np.arctan(y1_), dip_angle_gauss )
    
    #if a shift is presribed the point will be shifted in ther normal direction
    #with positive up
    if shift!= 0.:
        normals_rad_ = np.arctan(y1_) + np.pi/2.
        normals_rad = gaussian_filter(normals_rad_, dip_angle_gauss )

        even_distance = (even_distance - np.cos(normals_rad)*shift)
        even_depth = (even_depth - np.sin(normals_rad)*shift)
    
    #Deterimine the total distance along the curve
    arc_distance = curve_distance(np.column_stack((even_distance, even_depth))) + even_distance[0]
    
    
    return even_distance, even_depth, arc_distance, dip_angles
    
    


def scale_bar(ax, length=None, location=(0.5, 0.05), linewidth=3,  ):
    """
    ax is the axes to draw the scalebar on.
    length is the length of the scalebar in km.
    location is center of the scalebar in axis coordinates.
    (ie. 0.5 is the middle of the plot)
    linewidth is the thickness of the scalebar.
    """
    #Get the limits of the axis in lat long
    llx0, llx1, lly0, lly1 = ax.get_extent(ccrs.PlateCarree())
    #Make tmc horizontally centred on the middle of the map,
    #vertically at scale bar location
    sbllx = (llx1 + llx0) / 2
    sblly = lly0 + (lly1 - lly0) * location[1]
    tmc = ccrs.TransverseMercator(sbllx, sblly)
    #Get the extent of the plotted area in coordinates in metres
    x0, x1, y0, y1 = ax.get_extent(tmc)
    #Turn the specified scalebar location into coordinates in metres
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    #Calculate a scale bar length if none has been given
    #(Theres probably a more pythonic way of rounding the number but this works)
    if not length: 
        length = (x1 - x0) / 5000 #in km
        ndim = int(np.floor(np.log10(length))) #number of digits in number
        length = round(length, -ndim) #round to 1sf
        #Returns numbers starting with the list
        def scale_number(x):
            if str(x)[0] in ['1', '2', '5']: return int(x)        
            else: return scale_number(x - 10 ** ndim)
        length = scale_number(length) 

    #Generate the x coordinate for the ends of the scalebar
    bar_xs = [sbx - length * 500, sbx + length * 500]
    #Plot the scalebar
    ax.plot(bar_xs, [sby, sby], transform=tmc, color='k', linewidth=linewidth)
    #Plot the scalebar label
    ax.text(sbx, sby, str(length) + ' km', transform=tmc,
            horizontalalignment='center', verticalalignment='bottom',fontsize = 18)
