# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:09:56 2023

@author: tasmi
"""

import sys
import numpy as np

import xarray as xr
import rioxarray #Note, import not needed but install is needed to run code!

import geopandas as gpd, pandas as pd
from shapely.geometry import Point, LineString
import math

import simplekml

#station, observer_latitude, observer_longitude, observer_elevation, dem_path = sys.argv[1:]

def build_occluded_sightlines(station, observer_latitude, observer_longitude, observer_elevation, dem_path=None, \
                              azimuths=range(0, 360), RH=None, angles=range(5,31), d=25, savedir=None):
    '''
    For a given lat/lon/elevation location (station), check for blocked sight lines.
    Uses fixed angles/elevations for simplicity/speed of processing
    Only checks up to a 'd' distance away in km
    '''
    
    def extract_along_line(dem, line, n_samples=100):
        #Sample along a geographic line in n_samples points. 
        profile, sample_pts = [], []
        
        for i in range(1, n_samples + 1):
            point = line.interpolate(i / n_samples, normalized=True) #Get normalized distance along the line [0-1], excluding self start
            value = dem.sel(x=point.x, y=point.y, method="nearest").data[0] #Sample DEM at that point
            profile.append(value)
            sample_pts.append(Point(point.x, point.y)) #Return point as well
            
        return sample_pts, profile
    
    def get_point_at_distance(lat1, lon1, d, bearing, R=6371):
        from math import asin, atan2, cos, degrees, radians, sin

        #Via: https://stackoverflow.com/questions/7222382/get-lat-long-given-current-point-distance-and-bearing
        """
        lat: initial latitude, in degrees
        lon: initial longitude, in degrees
        d: target distance from initial
        bearing: (true) heading in degrees
        R: optional radius of sphere, defaults to mean radius of earth

        Returns new lat/lon coordinate {d}km from initial, in degrees
        """
        lat1 = radians(lat1)
        lon1 = radians(lon1)
        a = radians(bearing)
        lat2 = asin(sin(lat1) * cos(d/R) + cos(lat1) * sin(d/R) * cos(a))
        lon2 = lon1 + atan2(
            sin(a) * sin(d/R) * cos(lat1),
            cos(d/R) - sin(lat1) * sin(lat2))
        
        return degrees(lat2), degrees(lon2)
    
    #First get a local DEM via Copernicus, or use DEM
    if dem_path:
        dem = xr.open_dataset(dem_path)
    else:
        print('No DEM! Please provide a path to a local or global DEM')
        sys.exit()

    obs_ll = np.array([observer_longitude, observer_latitude]) #Set the observer (station) location
    
    #Store data in geodataframes
    gdf_list_point = []
    gdf_list_line = []
    
    #Now loop through all angles chosen
    for a in azimuths:
        if a % 30 == 0:
            print('Starting Azimuth', a)
        #Get a new point some distance away from observer
        target_lat, target_lon = get_point_at_distance(observer_latitude, observer_longitude, d, a)
        target_ll = np.array([target_lon, target_lat])
        
        #Build a line between the two
        line_ll = LineString([obs_ll, target_ll])
        
        #Get the DEM values along that line
        sample_pts, profile = extract_along_line(dem.band_data, line_ll, n_samples=50)
        
        #Find the minimum angle where we get topo blocking
        for i, pt in enumerate(sample_pts):
            topo_elev = profile[i]
            if np.isnan(topo_elev):
                topo_elev = 0 #Assume NaNs in DEM are water/sea level
            #Get distance between station/location to see how high up we are at a given angle
            dist = pt.distance(Point(obs_ll))
            if dist > 0: #Pass over self sample
                dist = dist * 111000 #Convert that to m (roughly), via 1dd = 111km
                
                #Now check various slope angles
                for ang in angles:
                    view_height = observer_elevation + float(RH) + dist * math.sin(ang * np.pi/180)
                    if view_height <= topo_elev:
                        print(a, ang, dist, view_height, topo_elev)
                        #Add to a list of problematic view points, coded with angle
                        loc = {'ElevAngle':ang, 'Azimuth': a, 'ViewHeight':view_height,'Topo':topo_elev}
                        gdf = gpd.GeoDataFrame(loc, geometry=[pt], crs='epsg:4326', index=[0])
                        gdf_list_point.append(gdf)

                        ln = LineString([pt, obs_ll])
                        gdf = gpd.GeoDataFrame(loc, geometry=[ln], crs='epsg:4326', index=[0])
                        gdf_list_line.append(gdf)

    #Stack everything together
    print('Collecting Blocking Points..')
    full_gdf = pd.concat(gdf_list_point) 
    out_gdf_point = gpd.GeoDataFrame(full_gdf, geometry=full_gdf.geometry, crs='epsg:4326')

    full_gdf = pd.concat(gdf_list_line)
    out_gdf_line = gpd.GeoDataFrame(full_gdf, geometry=full_gdf.geometry, crs='epsg:4326')
    
    if savedir:
        print('Saving KML output...')
        #Save out the data as geopackages if needed
        #out_gdf_point.to_file(savedir + station + '_blockedPoints.gpkg', driver='GPKG')
        #out_gdf_line.to_file(savedir + station + '_blockedLines.gpkg', driver='GPKG')
        
        #Export as KML
        def build_colors(n_colors=30, cmap='gist_rainbow'):
            from matplotlib import pyplot as plt
            from matplotlib.colors import rgb2hex
            cm = plt.get_cmap(cmap)
            return [rgb2hex(cm(1.*i/n_colors)) + 'FF' for i in range(n_colors)] #FF makes them opaque
            
        kml = simplekml.Kml()
        
        colormap = [simplekml.Color.yellow, simplekml.Color.blue, simplekml.Color.red,simplekml.Color.green,simplekml.Color.cyan,simplekml.Color.white] #build_colors(len(angles))
        a_arr = np.array(list(angles))
        #Split into two folders, one for lines one for points
#        fol = kml.newfolder(name='Points')
#        for ang in angles:
#            subset = out_gdf_point[out_gdf_point.ElevAngle == ang]
#            #Define one style for all points of a given angle to save space in the output KML
#            sharedstyle = simplekml.Style()
#            color = colormap[np.where(a_arr == ang)[0][0]]
#            sharedstyle.iconstyle.color = color
#            sharedstyle.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'
#            sharedstyle.altitudemode = simplekml.AltitudeMode.relativetoground
            
            for _, row in subset.iterrows():
                #pnt = fol.newpoint(name=ang, coords = row.geometry.coords[:]) #There can be a LOT of names here...
                pnt = fol.newpoint(coords = row.geometry.coords[:])
                pnt.style = sharedstyle 
            del sharedstyle
            
        fol = kml.newfolder(name='Lines')
        for ang in angles:
            subset = out_gdf_line[out_gdf_line.ElevAngle == ang]
            #Define one style for all lines of a given angle
            sharedstyle = simplekml.Style()
            color = colormap[np.where(a_arr == ang)[0][0]]
            sharedstyle.linestyle.color = color
            sharedstyle.iconstyle.color = color
            sharedstyle.linestyle.width = 3
            sharedstyle.altitudemode = simplekml.AltitudeMode.relativetoground
            
            for _, row in subset.iterrows():
                ln = fol.newlinestring(coords = row.geometry.coords[:])
                ln.style = sharedstyle 
            del sharedstyle

        kml.save(savedir + station + '_TopoBlocking.kml')
    
    #Return points if want to use them for something futher
    return out_gdf_point, out_gdf_line

#Test run based on Global DEM
#https://www.hydrosheds.org/hydrosheds-core-downloads
#dem_path = 'hyd_glo_dem_15s.tif'

out_gdf_point, out_gdf_line = build_occluded_sightlines(station, observer_latitude, observer_longitude, \
                                                        observer_elevation, dem_path)


    
