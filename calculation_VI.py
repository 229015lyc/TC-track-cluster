#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 15:20:18 2021

@author: cwuyichen
 calculation of ventilation index (Tang and Emanual 2012) 
 some def follow tcpyPI.utilities (Gilford, D. M. 2020: pyPI: Potential Intensity Calculations in Python, pyPI v1.3. Zenodo)
"""
import xarray as xr

import numpy as np
import numba as nb
from tcpyPI.utilities import *

syear = 1979
eyear = 2004

#%%
for year in np.arange(syear,eyear+1):
    print(year)

    #  read --------------------------------------------------------
    print('        reading........')
    
    # pi and entropy dicifit -------------------
    diri = './'
    subdiri = 'clim19792004'
    
    # only JJASON
    fn = diri+'{}/output_{}JJASON.nc'.format(subdiri,year)
    
    ds = xr.open_dataset(fn)
    
    pi = ds['vmax'].data
    
    diseq = ds['diseq']
    sst = ds['sst']
    sst = T_Ctok(sst.data)
    s_denominator = diseq.data/sst.data
    
    
    
    q_600 = ds['q'].sel(p=600)
    t_600 = ds['t'].sel(p=600)
    q_600 = q_600.data/1000.
    t_600 = T_Ctok(t_600.data)
    p = np.full_like(q_600,600.)
    
    s_600 = entropy_S(t_600,q_600,p)
    es_600_star = es_cc(T_ktoC(t_600))
    q_600_star = rv(es_600_star,p)
    s_600_star = entropy_S(t_600,q_600_star,np.full_like(q_600,600.))
    s_numerator = s_600_star-s_600
    
    #del sst,diseq,q_600,t_600,q_600_star,s_600,s_600_star
    ds.close()

    # VWS---------------------------
    diri = './VWS200850/'

    mon = np.arange(6,12)
    lat_new = np.linspace(-89.5, 89.5, 180)
    lon_new = np.linspace(0.5, 359.5, 360)
    

    for fn_m in mon:
        fn = diri+'{}/{}_vws200850_{}{}.nc'.format(subdiri,subdiri[:-13],year,str(fn_m).zfill(2))
        #print(fn)
        ds_input = xr.open_dataset(fn)
        
        ds_new = ds_input.interp(grid_xt=lon_new,grid_yt=lat_new)
        ds_new = ds_new.rename({'grid_xt':'lon','grid_yt':'lat'})
        ds_new['month'] = fn_m
        ds_new = ds_new.set_coords('month')
        ds_new = ds_new.expand_dims({'month': 1})    
        ds_input.close()
        
        if fn_m == mon[0]:
            ds = ds_new
        else:
            ds = xr.combine_by_coords([ds,ds_new]) 
    
        ds_new.close()
    del ds_input,ds_new
      
    ds = ds.roll(lon=180)
    ds['lon'] = np.linspace(-179.5,179.5,360)   
    vws = ds['vws'].data
    ds.close()  

    #  calculation vi ----------------------------------------------------
    print('        calculating........')
    vi = vws*s_numerator/s_denominator/pi
    
    ds_out = xr.Dataset(
            data_vars=dict(
                vi=(['month','lat','lon'],vi,{'standard_name':'Ventilation Index','units':'NaN'}),
                vws=(['month','lat','lon'],vws,{'standard_name':'Veritcal Wind Shear 200hPa-850hPa','units':'m/s'}),
                s_nu=(['month','lat','lon'],s_numerator,\
                      {'standard_name':'Entropy decifit (numerator)','units':'J/K/kg',\
                       'description':'computed follow (Tang and Emanual 2012) (3) at 600hPa'}),
                s_de=(['month','lat','lon'],s_denominator,\
                      {'standard_name':'Entropy decifit (denominator)','units':'J/K/kg',\
                       'description':'disequilibrium/SST at 600hPa'}),
                pi=(['month','lat','lon'],pi,{'standard_name':'Potential Intensity','units':'m/s'}),
            ),
            coords=dict(
                month=mon,
                lat=lat_new,
                lon=np.linspace(-179.5,179.5,360)
            )
    )
        
    ds_out['lon'].attrs = {'standard_name':'longitude','units':'degree_east'}
    ds_out['lat'].attrs = {'standard_name':'latitude','units':'degree_north'}
    ds_out['month'].attrs = {'standard_name':'Month','units':'Month Number'}
    
    #  save -----------------------------------------------
    diriout = './output_vi/{}/'.format(subdiri)
    print('save to  ',diriout)
    ds_out.to_netcdf(diriout+"vi_{}JJASON.nc".format(str(year)))

