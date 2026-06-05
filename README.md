# LMDZ Patches

LMDZ versions need minimal modifications to create a module dedicated to communication with coupled Python scripts. This repository contains patches of modified LMDZ sources for versions: 

- [LMDZ6](http://svn.lmd.jussieu.fr/LMDZ/LMDZ6/) (trunk 5891)

## Use a patch 
- Install LMDZ6 with parallel MPI and OMP features
- Copy patch sources in the LMDZ `libf` directory, overwrite already existing sources
- Import and use Python communication module API (more details [here](https://morays-doc.readthedocs.io/en/latest/nemo.api_4.html#user-guide), template is for NEMO but works the same with LMDZ)

**NB:** Fields to be passed through the API must be defined on the dynamics MPI grid, which might imply to use ```cpl2gath()``` and ```gath2cpl()``` functions to swap the fields from physics to dynamics grid, and conversely.
- Compile with ```-cpp key_eophis -c true``` flags and OASIS_v5.0 (see this [guide](https://morays-doc.readthedocs.io/en/latest/nemo.getting_started.html#morays-environment))

## Patch Modifications
  * Architecture: OASIS coupling module `oasis.F90` was initially managed within ocean coupling dedicated routines
      - OASIS configuration is now totally managed by ```phyetat0.F90``` main routine
      - Coupling module `oasis.F90` is independent and can be called by any other module to define, send or receive coupling variables
      - Ocean coupling definition is moved back up to `cpl_mod.f90` and leverages new module-independent `oasis.F90` interface
      - Possible to perform exchange of 3D fields (OASIS_v5.0 or later required)

  * New modules:        
      - `eophis_def.F90` : reads Eophis Fortran namelist to give NEMO access to Eophis script attributes
      - `pycpl.F90` : configure NEMO coupler layer (OASIS) from `eophis_def` - exposes Python communication API
      - `pyfld.F90` : used to define and store fields returned by the Python script
