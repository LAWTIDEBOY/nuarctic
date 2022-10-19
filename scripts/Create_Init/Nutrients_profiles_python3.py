#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!jupyter nbconvert --to=python Nutrients_profile_python3.ipynb


# In[ ]:


class Nut_profile:   
    '''
    class Nut_depth(runname,resultpath,savepath,meshpath,ncfileDSi,ncfileDIN,ncfileDFe,first_year,last_year,
                 savefig=False,regional=True)
                 
    c.f. Danilov et al. (2017):
    "in the vertical direction, the horizontal velocities and scalars are located at mid-levels" 
                 
    if regional = True, profiles will plotted for each main basins + Global Ocean. 
    Otherwise, just the Global Ocean.
    '''
    def __init__(self,runname,mesh,ncfileDSi,ncfileDIN,ncfileDFe):

        self.runname = runname
        self.mesh = mesh
        self.ncfileDSi = ncfileDSi
        self.ncfileDIN = ncfileDIN
        self.ncfileDFe = ncfileDFe
        
        import matplotlib.pyplot as plt
        import numpy as np
        from netCDF4 import Dataset
        from scipy.interpolate import griddata
        import skill_metrics as sm
        import cartopy.crs as ccrs
        import pyfesom2 as pf
        
        from load_interp_WOA_python3 import WOAdata
        from load_interp_PISCES_python3 import PISCESdata
        
        
        meshdiag = pf.get_meshdiag(mesh)
        runid      =  self.runname
        
        unitsDIN = 'DIN [mmol m$^{-3}$]'
        unitsDSi = 'DSi [mmol m$^{-3}$]'
        unitsDFe = 'DFe [mmol m$^{-3}$]'
        
        
        # load data -------------------------------------------------------------------------------------
        DINwoa_input = WOAdata(runid,mesh,ncfileDIN,'n_an', get_overview=False)
        DSiwoa_input = WOAdata(runid,mesh,ncfileDSi,'i_an', get_overview=False)
        DFepisces_input = PISCESdata(runid,mesh,ncfileDFe,'Fe', get_overview=False)
        
        DFepisces = DFepisces_input.pisces_int
        DFepisces[DFepisces>1000]= 0
        #DFepisces[DFefesom == 0] = 0
        
        DINwoa = DINwoa_input.woa_int
        #DINwoa[DINfesom == 0] = 0

        DSiwoa = DSiwoa_input.woa_int
        #DSiwoa[DSifesom == 0] = 0
        
        # Load and derive profiles

        nod_area = np.ma.masked_equal(meshdiag.nod_area.values, 0)
        mask = pf.get_mask(mesh, "Arctic_Basin")

        DINwoa_by_area = ((np.ma.masked_equal(DINwoa[mask,:],0) * nod_area[:-1,:].T[mask]).mean(axis=0))
        DINwoa_weighted_Arctic = DINwoa_by_area/nod_area[:-1,:].T[mask].mean(axis=0)

        DSiwoa_by_area = ((np.ma.masked_equal(DSiwoa[mask,:],0) * nod_area[:-1,:].T[mask]).mean(axis=0))
        DSiwoa_weighted_Arctic = DSiwoa_by_area/nod_area[:-1,:].T[mask].mean(axis=0)
            
        DFepisces_by_area = ((np.ma.masked_equal(DFepisces[mask,:],0) * nod_area[:-1,:].T[mask]).mean(axis=0))
        DFepisces_weighted_Arctic = DFepisces_by_area/nod_area[:-1,:].T[mask].mean(axis=0)
            
        self.DFe = DFepisces_weighted_Arctic
        self.DSi = DSiwoa_weighted_Arctic
        self.DIN = DINwoa_weighted_Arctic
        self.depth = DINwoa_input.layer_depths

