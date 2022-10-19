#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!jupyter nbconvert --to=python Nutrients_profile_python3.ipynb


# In[ ]:


class DO2_profile:   
    '''
    class Nut_depth(runname,resultpath,savepath,meshpath,ncfileDSi,ncfileDO2,ncfileDFe,first_year,last_year,
                 savefig=False,regional=True)
                 
    c.f. Danilov et al. (2017):
    "in the vertical direction, the horizontal velocities and scalars are located at mid-levels" 
                 
    if regional = True, profiles will plotted for each main basins + Global Ocean. 
    Otherwise, just the Global Ocean.
    '''
    def __init__(self,runname,mesh,ncfileDO2):

        self.runname = runname
        self.mesh = mesh
        self.ncfileDO2 = ncfileDO2
        
        import matplotlib.pyplot as plt
        import numpy as np
        from netCDF4 import Dataset
        from scipy.interpolate import griddata
        import skill_metrics as sm
        import cartopy.crs as ccrs
        import pyfesom2 as pf
        
        from load_interp_WOA_python3 import WOAdata
        from load_interp_PISCES_python3 import PISCESdata
        
        #mesh       = pf.load_mesh(self.meshpath)
        meshdiag = pf.get_meshdiag(mesh)
        runid      =  self.runname
        
        unitsDO2 = 'DO2 [mmol m$^{-3}$]'
        
        
        # load data -------------------------------------------------------------------------------------
        DO2woa_input = WOAdata(runid,mesh,ncfileDO2,'oxygen_mmol', get_overview=False)

        DO2woa = DO2woa_input.woa_int
        #DO2woa[DO2fesom == 0] = 0
        
        # Load and derive profiles

        nod_area = np.ma.masked_equal(meshdiag.nod_area.values, 0)
        mask = pf.get_mask(mesh, "Global Ocean")

        DO2woa_by_area = ((np.ma.masked_equal(DO2woa[mask,:],0) * nod_area[:-1,:].T[mask]).mean(axis=0))

        DO2woa_weighted_Global = DO2woa_by_area/nod_area[:-1,:].T[mask].mean(axis=0)
            
        mask = pf.get_mask(mesh, "Arctic_Basin")

        DO2woa_by_area = ((np.ma.masked_equal(DO2woa[mask,:],0) * nod_area[:-1,:].T[mask]).mean(axis=0))

        DO2woa_weighted_Arctic = DO2woa_by_area/nod_area[:-1,:].T[mask].mean(axis=0)
        
        self.DO2 = DO2woa_weighted_Arctic
        self.depth = DO2woa_input.layer_depths
        
            
#         fig, axs = plt.subplots(1,1, figsize=(10, 5), facecolor='w', edgecolor='k', constrained_layout=True, sharey=True)

#         print(axs)
#         axs[0].plot(DO2woa_weighted_Arctic, mesh.zlev[:-1]/1000,label = 'WOA', color = 'k', lw=3, linestyle = '--')
#         axs[0].set_ylabel('Depth [km]',fontsize=14)
#         axs[0].set_xlabel(unitsDO2,fontsize=14)
#         axs[1].set_title('Arctic Ocean',size=16, weight='bold')
#         axs[0].tick_params(labelsize=14)
#         axs[0].grid()
#         axs[0].legend(loc='best', borderaxespad=0., fontsize=14)

#         plt.show(block=False)

        

