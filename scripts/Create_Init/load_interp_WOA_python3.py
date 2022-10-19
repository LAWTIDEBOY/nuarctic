#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!jupyter nbconvert --to=python load_interp_WOA_python3.ipynb


# In[ ]:


class WOAdata:
    '''
    Load World Ocean Atlas NetCDF4 file and interpolate to FESOM mesh.
    
    Careful: Output still needs to be masked with an ocean bathymetry mask, e.g. with FESOM output!
    
    In: ncfile,ncvariable,meshpath
    
    Out: class WOAdata with self.woa_int containing interpolated WOA ncvariable
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

        # load NetCDF ------------------------------------------------------------------------------------
        print('***\nLoading WOA file: {0}\n***'.format(self.ncfile))
        f          = Dataset(self.ncfile, 'r')
        DepthRaw   = -f.variables['depth'][:]                                
        lonwoa     =  f.variables['lon'][:]
        lonwoa[lonwoa>180]=lonwoa[lonwoa>180]-360
        latwoa     =  f.variables['lat'][:]
        Timewoa     =  f.variables['time'][:]
        VARwoa_temp     =  f.variables[ncvariable][:]
        VARwoa_temp     = np.squeeze(VARwoa_temp)
        VARwoa_temp     = np.ma.filled(VARwoa_temp, np.nan)
        
        VARwoa     = np.zeros(shape=(np.shape(VARwoa_temp)))                         # Change longitude from 0:360 to -180:180 
        for i in range(0,len(DepthRaw)):
          VARwoa[i,:,:] = np.hstack((VARwoa_temp[i,:,181:360],VARwoa_temp[i,:,0:181]))
        
        x360       = np.arange(-179.5,180.,1.)
        y180       = np.arange(-89.5,89.6,1.)
        X360, Y180 = np.meshgrid(x360, y180)
        #X360, Y180 = np.meshgrid(lonwoa, latwoa)

        if(self.get_overview==True):
            #!ncdump -h $self.ncfile
            
            fig = plt.figure(figsize= (7,7))
            ax = plt.subplot()
            im = ax.pcolor(X360, Y180, VARwoa[0,:,:])
            cbar = fig.colorbar(im, orientation = 'horizontal')
            cbar.set_label(ncvariable) 
            plt.title('WOA var "{0} before interpolation"'.format(self.ncvariable))

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
        
        # check maximum depth in WOA compared to FESOM
        dmin_woa = np.min(DepthRaw)
        dmin_fesom = np.min(mesh_depths)#mesh.zlev)

        if(dmin_woa <= dmin_fesom):
            print('***\nDepth greater in WOA ({0}) than in FESOM ({1})'.format(dmin_woa, dmin_fesom))
            ilev = len(mesh_depths)
            max_zlev = mesh_depths[ilev-1]
        else:
            print('***\nDepth greater in FESOM ({1}) than in WOA ({0})'.format(dmin_woa, dmin_fesom))
            ilev = np.where(mesh_depths >= dmin_woa)
            ilev = ilev[0][-1]
            max_zlev = mesh_depths[ilev]

        # storage container
        woa_int = np.zeros((len(mesh.zlev)-1,len(mesh.x2)))
        #print(np.shape(din_int))

        for k in range(0,len(mesh_depths)): # depth are layer depth --> between levels, c.f. FESOM documentation
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
            if False: #(self.get_overview == True):
                print('\nFESOM depth = {0}, WOA depths = {1}, {2} \nDepth indices: {3} {4},  FESOM index: {5} \nScaling c1 = {6}, c2 = {7}'.format(lev,DepthRaw[ind1],DepthRaw[ind2],ind1, ind2,indZ,c1,c2))


            aux1  = VARwoa[ind1,:,:]                     # Find the data above the current fesom depth
            aux2  = VARwoa[ind2,:,:]                     # Find the data below the current fesom depth
            aux   = np.squeeze(c1*aux1+c2*aux2)          # Scaling the data according to vertical distribution as found above
            ind   = np.squeeze(~np.isnan(aux)) 
            #print(np.shape(aux), np.shape(ind))

            # first interpolation to original grid to close empty parts
            aux = griddata((X360[ind], Y180[ind]), aux[ind], (X360, Y180), method='nearest')                             
            # 2D field without nans                           

            # second interpolation to FESOM grid
            woa_int[indZ,:] = griddata((X360.ravel(), Y180.ravel()), aux.ravel(), 
                                   (mesh.x2, mesh.y2), method='nearest')  
            # Final interpolated field

            if np.isnan(np.min(woa_int)): print('WARNING: The interpolated field contains NaNs at depth',lev)                 # Testing if results contain NaNs. If yes, the routine needs adjustments

            if(self.get_overview ==True):
                print('Depth: {0} min = {1} max = {2} mean = {3}'.format(lev,np.min(woa_int), np.max(woa_int), np.mean(woa_int)))

        woa_int = np.swapaxes(woa_int,0,1) # adjust axes layout to FESOM output
        #print(np.shape(woa_int))    
        
        self.layer_depths = mesh_depths
        self.woa_int = woa_int

