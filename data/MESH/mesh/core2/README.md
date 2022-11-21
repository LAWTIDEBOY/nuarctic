# FESOM2 meshes

This repository contains one of standard FESOM2 meshes. The latest versions can be found on AWI GitLab, in public FESOM group (https://gitlab.awi.de/fesom/)

All meshes are not rotated (provided in geographical coordinates). Use of meshes in rotated coordinates is no longer recomended.

Please cite original publications where meshes are introduced if you use them in your work.


## core2

**download:** `git lfs clone https://gitlab.awi.de/fesom/core2.git` 

**source mesh:** core2_meanz

**topography:** A blend of several bottom topography data sets (see Wang et al., 2014).

**Reference (please cite when you use this mesh):**
 
* Wang, Q., Danilov, S., Sidorenko, D., Timmermann, R., Wekerle, C., Wang, X., Jung, T., and Schröter, J.: The Finite Element Sea Ice-Ocean Model (FESOM) v.1.4: formulation of an ocean general circulation model, Geosci. Model Dev., 7, 663–693, https://doi.org/10.5194/gmd-7-663-2014, 2014. 


## farc

**download:** `git lfs clone https://gitlab.awi.de/fesom/farc.git`

**source mesh:** farc_sorted

**contact:** W. Wang

**topography:** IBCAO and GEBCO (see details in Wang et al., 2018)

**Reference (please cite when you use this mesh):**

* Wang, Q., Wekerle, C., Danilov, S., Wang, X., and Jung, T.: A 4.5 km resolution Arctic Ocean simulation with the global multi-resolution model FESOM 1.4, Geosci. Model Dev., 11, 1229–1255, https://doi.org/10.5194/gmd-11-1229-2018, 2018. 

**! NOTE**: The nodes in this mesh was sorted, which is different from the version of Wang et al., 2018


## mr

**download:** `git lfs clone https://gitlab.awi.de/fesom/mr.git`

**source mesh:** mr

**contact:** D. Sein

**Reference (please cite when you use this mesh):**

* Sein, D. V., Danilov, S., Biastoch, A., Durgadoo, J. V., Sidorenko, D., Harig, S., and Wang, Q. (2016), Designing variable ocean model resolution based on the observed ocean variability, J. Adv. Model. Earth Syst., 8, 904– 916, doi:10.1002/2016MS000650. 

* Sein, D. V., Koldunov, N. V., Danilov, S., Wang, Q., Sidorenko, D., Fast, I., Rackow, T., Cabos, W. and Jung, T., Ocean Modeling on A Mesh with Resolution Following the Local Rossby Radius. Journal of Advances in Modeling Earth Systems, 9, 2601–2614. doi:10.1002/2017MS001099 , 2017

**For AWI-CM1 setup also cite:**

* Semmler, T., Danilov, S., Gierz, P., Goessling, H. F., Hegewald, J., Hinrichs, C., et al. (2020). Simulations for CMIP6 with the AWI climate model AWI‐CM‐1‐1. Journal of Advances in Modeling Earth Systems, 12, e2019MS002009. https://doi.org/10.1029/2019MS002009 


## orca25

**download:** `git lfs clone https://gitlab.awi.de/fesom/orca25.git`

**source mesh:** NEMO_RT_FIXED

**contact:** N. Koldunov

**topography:** RTopo-2.0.1_30sec_global_topo_2016-12-13.nc (https://doi.org/10.5194/essd-8-543-2016) with some modifications.

**Reference (please cite when you use this mesh):** Juricke, S., Danilov, S., Koldunov, N. V., Oliver, M., & Sidorenko, D. (2020). Ocean kinetic energy backscatter parametrization on unstructured grids: Impact on global eddy‐permitting simulations. Journal of Advances in Modeling Earth Systems, 12, e2019MS001855. https://doi.org/10.1029/2019MS001855



## hr

**download:** `git lfs clone https://gitlab.awi.de/fesom/hr.git`

**source mesh:** BOLD_RT_FIXED

**contact:** D. Sein

**topography:** RTopo-2.0.1_30sec_global_topo_2016-12-13.nc (https://doi.org/10.5194/essd-8-543-2016) with some modifications.

**Reference (please cite when you use this mesh):**
* Sein, D. V., Danilov, S., Biastoch, A., Durgadoo, J. V., Sidorenko, D., Harig, S., and Wang, Q. (2016), Designing variable ocean model resolution based on the observed ocean variability, J. Adv. Model. Earth Syst., 8, 904– 916, doi:10.1002/2016MS000650. 

* Sein, D. V., Koldunov, N. V., Danilov, S., Wang, Q., Sidorenko, D., Fast, I., Rackow, T., Cabos, W. and Jung, T., Ocean Modeling on A Mesh with Resolution Following the Local Rossby Radius. Journal of Advances in Modeling Earth Systems, 9, 2601–2614. doi:10.1002/2017MS001099 , 2017

**For AWI-CM2 setup also cite:**
* Sidorenko, D., Goessling, H. F., Koldunov, N. V., Scholz, P., Danilov, S., Barbi, D., et al (2019). Evaluation of FESOM2.0 coupled to ECHAM6.3: Pre‐industrial and HighResMIP simulations. Journal of Advances in Modeling Earth Systems, 11. https://doi.org/10.1029/2019MS001696

**! NOTE:** Topography used for the current mesh is different from the ones used in Sein et al., 2016, 2017. This exact version of the mesh was used in Sidorenko et al., 2019.



## pi

**download:** `git lfs clone https://gitlab.awi.de/fesom/pi.git`

**source mesh:** pi_mesh

**contact:** T. Rakow

**topography:**

**Reference:**



## core2_old

**download:** `git lfs clone https://gitlab.awi.de/fesom/core2_old.git`

**source mesh:** mesh_CORE2_finaltopo_mean

**topography:** A blend of several bottom topography data sets (see Wang et al., 2014).

**Reference (please cite when you use this mesh):**

* Wang, Q., Danilov, S., Sidorenko, D., Timmermann, R., Wekerle, C., Wang, X., Jung, T., and Schröter, J.: The Finite Element Sea Ice-Ocean Model (FESOM) v.1.4: formulation of an ocean general circulation model, Geosci. Model Dev., 7, 663–693, https://doi.org/10.5194/gmd-7-663-2014, 2014.

**! NOTE:** The only reason to use this mesh is to reproduce AWI-CM2 setup from Sidorenko et al., 2019. In other cases **DONT USE IT!**

