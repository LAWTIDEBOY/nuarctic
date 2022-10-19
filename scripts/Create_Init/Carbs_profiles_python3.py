#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!jupyter nbconvert --to=python Carbs_profile_python3.ipynb


# In[ ]:


class Carbs_profile:   
    '''
    class Nut_depth(runname,resultpath,savepath,meshpath,ncfileAlk,ncfileDIC,ncfileDFe,first_year,last_year,
                 savefig=False,regional=True)
                 
    c.f. Danilov et al. (2017):
    "in the vertical direction, the horizontal velocities and scalars are located at mid-levels" 
                 
    if regional = True, profiles will plotted for each main basins + Global Ocean. 
    Otherwise, just the Global Ocean.
    '''
    def __init__(self,runname,mesh,ncfileAlk,ncfileDIC):

        self.runname = runname
        self.mesh = mesh
        self.ncfileAlk = ncfileAlk
        self.ncfileDIC = ncfileDIC
        
        import matplotlib.pyplot as plt
        import numpy as np
        from netCDF4 import Dataset
        from scipy.interpolate import griddata
        import skill_metrics as sm
        import cartopy.crs as ccrs
        import pyfesom2 as pf
        
        from load_interp_GLODAP_python3 import GLODAPdata
        
        meshdiag = pf.get_meshdiag(mesh)
        runid      =  self.runname
        
        unitsDIC = 'DIC [mmol m$^{-3}$]'
        unitsAlk = 'Alk [mmol m$^{-3}$]'
        
        
        # load data -------------------------------------------------------------------------------------
        DICglodap_input = GLODAPdata(runid,mesh,ncfileDIC,'TCO2_mmol', get_overview=False)
        Alkglodap_input = GLODAPdata(runid,mesh,ncfileAlk,'TAlk_mmol', get_overview=False)
        
        DICglodap = DICglodap_input.glodap_int
        #DICglodap[DICfesom == 0] = 0

        Alkglodap = Alkglodap_input.glodap_int
        #Alkglodap[Alkfesom == 0] = 0
        
        # Load and derive profiles
        nod_area = np.ma.masked_equal(meshdiag.nod_area.values, 0)
        mask = pf.get_mask(mesh, "Arctic_Basin")

        DICglodap_by_area = ((np.ma.masked_equal(DICglodap[mask,:],0) * nod_area[:-1,:].T[mask]).mean(axis=0))
        DICglodap_weighted_Arctic = DICglodap_by_area/nod_area[:-1,:].T[mask].mean(axis=0)

        Alkglodap_by_area = ((np.ma.masked_equal(Alkglodap[mask,:],0) * nod_area[:-1,:].T[mask]).mean(axis=0))
        Alkglodap_weighted_Arctic = Alkglodap_by_area/nod_area[:-1,:].T[mask].mean(axis=0)
            
        self.Alk = Alkglodap_weighted_Arctic
        self.DIC = DICglodap_weighted_Arctic
        self.depth = DICglodap_input.layer_depths

