#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!jupyter nbconvert --to=python load_interp_PISCES_python3.ipynb


# In[ ]:


class PISCESdata:
    '''
    Load PISCES Atlas NetCDF4 file and interpolate to FESOM mesh.
    
    Careful: Output still needs to be masked with an ocean bathymetry mask, e.g. with FESOM output!
    
    In: ncfile,ncvariable,meshpath
    
    Out: class PISCESdata with self.pisces_int containing interpolated PISCES ncvariable
    '''
    
    def __init__(self,runname,mesh,ncfile,ncvariable,get_overview=False):
    
        self.runid = runname
        self.mesh = mesh
        self.ncfile = ncfile
        self.ncvariable = ncvariable
        self.get_overview = get_overview

        import matplotlib.pyplot as plt
        import numpy as np
        from netCDF4 import Dataset
        from scipy.interpolate import griddata
        #import skill_metrics as sm
        #import cartopy.crs as ccrs
        #import pickle

        import pyfesom2 as pf
        
        # load FESOM mesh ------------------------------------------------------------------------------------
        #mesh = pf.load_mesh(self.meshpath)
        
        # load FESOM mesh diag -------------------------------------------------------------------------------
        meshdiag= self.mesh.path+'/'+self.runid+'.mesh.diag.nc'
        #!ncdump -h $meshdiag

        diag = pf.get_meshdiag(mesh,meshdiag=meshdiag, runid=self.runid)
        if 'nl' in diag.dims.mapping.keys():
            #nod_area = diag.rename_dims({"nl": "nz1", "nod_n": "nod2"}).nod_area
            #nod_area.load()
            mesh_depths = diag['Z'].values
        elif 'nzl' in diag.dims.mapping.keys():
            mesh_depths = diag['nz1'].values
        
        # load raw PISCES data -------------------------------------------------------------------------------------

        f          = Dataset(self.ncfile, 'r')
        DepthRaw   = -f.variables['depth'][:]                                # Depth is negative
        DFe        =  f.variables[self.ncvariable][:]                         # Unit [mol/L]
        DFe        = 1.e9 * DFe                                              # [mol/L] => [umol/m3] 
        DFe        = np.ma.filled(DFe, np.nan)                               # From masked array to numpy array

        DFeRaw     = np.zeros(shape=(np.shape(DFe)))                         # Change longitude from 0:360 to -180:180 
        for i in range(0,len(DepthRaw)):
          DFeRaw[i,:,:] = np.hstack((DFe[i,:,181:360],DFe[i,:,0:181]))

        x360       = np.arange(-179.5,180.,1.)
        y180       = np.arange(-89.5,89.6,1.)
        X360, Y180 = np.meshgrid(x360, y180)
        
        # interpolate PISCES data -------------------------------------------------------------------------------------

        # check maximum depth in PISCES compared to FESOM
        dmin_woa = np.min(DepthRaw)
        dmin_fesom = np.min(mesh_depths)

        if(dmin_woa <= dmin_fesom):
            print('***\nDepth greater in PISCES ({0}) than in FESOM ({1})'.format(dmin_woa, dmin_fesom))
            ilev = len(mesh_depths)
            max_zlev = mesh_depths[ilev-1]
        else:
            print('***\nDepth greater in FESOM ({1}) than in PISCES ({0})'.format(dmin_woa, dmin_fesom))
            ilev = np.where(mesh_depths >= dmin_woa)
            ilev = ilev[0][-1]
            max_zlev = mesh_depths[ilev]

        # storage container
        pisces_int = np.zeros((len(mesh_depths),len(mesh.x2)))
        #print(np.shape(din_int))

        for k in range(0,len(mesh_depths)): # layer depth as in meshi diag.Z
            lev = mesh_depths[k] # current FESOM depth
            ind1 = np.where(DepthRaw >= lev)
            ind1 = ind1[0][-1]
            ind2 = np.where(DepthRaw < lev)[0]

            if ind2.size > 0:                            # If we have not yet reached the bottom
                ind2 = ind2[0]                           # The index of the depth level below the current fesom level
                c    = DepthRaw[ind1]-DepthRaw[ind2]     # Difference in depth between the data value above and below the fesom depth
                c1   = DepthRaw[ind1]-lev                # Difference between fesom depth and data depth above
                c2   = -(DepthRaw[ind2]-lev)             # Difference between fesom depth and data depth below
                c1   = (c-c1)/c                          # Scaling coefficient for the depth above
                c2   = (c-c2)/c                          # Scaling coefficient for the depth below
            else:                                        # We have reached the bottom
                c1   = 1.
                c2   = 0.
                ind2 = ind1

            indZ  = np.where(mesh_depths == lev)                               
            # original code:
            # indZ  = np.where(-mesh.z3 == lev)          # Find the mesh index of the current fesom depth
            indZ = indZ[0] 
            if(False):
                print('\nFESOM depth = {0}, WOA depths = {1}, {2} \nDepth indices: {3} {4},  FESOM index: {5} \nScaling c1 = {6}, c2 = {7}'.format(lev,DepthRaw[ind1],DepthRaw[ind2],ind1, ind2,indZ,c1,c2))


            aux1  = DFeRaw[ind1,:,:]                     # Find the data above the current fesom depth
            aux2  = DFeRaw[ind2,:,:]                     # Find the data below the current fesom depth
            aux   = np.squeeze(c1*aux1+c2*aux2)          # Scaling the data according to vertical distribution as found above
            ind   = np.squeeze(~np.isnan(aux)) 
            #print(np.shape(aux), np.shape(ind))

            # first interpolation to original grid to close empty parts
            aux = griddata((X360[ind], Y180[ind]), aux[ind], (X360, Y180), method='nearest')                             
            # 2D field without nans                           

            # second interpolation to FESOM grid
            pisces_int[indZ,:] = griddata((X360.ravel(), Y180.ravel()), aux.ravel(), 
                                   (mesh.x2, mesh.y2), method='nearest')  
            # Final interpolated field

            if np.isnan(np.min(pisces_int)): print('WARNING: The interpolated field contains NaNs at depth',lev)                 # Testing if results contain NaNs. If yes, the routine needs adjustments

            if(False):
                print('Depth: {0} min = {1} max = {2} mean = {3}'.format(lev,np.min(pisces_int), np.max(pisces_int), np.mean(pisces_int)))

        pisces_int = np.swapaxes(pisces_int,0,1) # adjust axes layout to FESOM output   
        
        self.layer_depths = mesh_depths
        self.pisces_int = pisces_int

