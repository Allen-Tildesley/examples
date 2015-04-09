# SCons build file for the various programs
# Run this typing 'scons'
# Clean the program by 'scons -c' (i.e., remove object, module and executable)

import os, sys

env_normal=Environment(ENV=os.environ)
env_mpi=Environment(ENV=os.environ,F90='mpif90',F90LINKER='mpif90')

# Assume that gfortran will be used. 
#Tool('gfortran')(env_normal)
#Tool('mpif90')(env_mpi)
#MY_F90FLAGS='-O2 -finline-functions -funswitch-loops -fwhole-file'
MY_F90FLAGS='-fdefault-real-8 -std=f2008 -Wall'
MY_LINKFLAGS=''

env_normal.Append(F90FLAGS=MY_F90FLAGS,LINKFLAGS=MY_LINKFLAGS,FORTRANMODDIRPREFIX='-J',FORTRANMODDIR = '${TARGET.dir}',F90PATH='${TARGET.dir}')
env_mpi.Append(F90FLAGS=MY_F90FLAGS,LINKFLAGS=MY_LINKFLAGS,FORTRANMODDIRPREFIX='-J',FORTRANMODDIR = '${TARGET.dir}',F90PATH='${TARGET.dir}')

variants={}
variants['build_mc_nvt_lj']    = (['mc_nvt_lj.f90','mc_lj_module.f90','utility_module.f90'],env_normal)
variants['build_mc_npt_lj']    = (['mc_npt_lj.f90','mc_lj_module.f90','utility_module.f90'],env_normal)
variants['build_mc_zvt_lj']    = (['mc_zvt_lj.f90','mc_lj_module.f90','utility_module.f90'],env_normal)
variants['build_md_nve_hs']    = (['md_nve_hs.f90','md_nve_hs_module.f90','utility_module.f90'],env_normal)
variants['build_mc_nvt_hs']    = (['mc_nvt_hs.f90','mc_hs_module.f90','utility_module.f90'],env_normal)
variants['build_mc_nvt_sc']    = (['mc_nvt_sc.f90','mc_sc_module.f90','utility_module.f90'],env_normal)
variants['build_md_nve_lj']    = (['md_nve_lj.f90','md_lj_module.f90','utility_module.f90'],env_normal)

# Build each variant in appropriate variant directory
for variant_dir,(sources,env) in variants.iteritems():
    SConscript('SConscript', variant_dir=variant_dir,exports={'env':env,'sources':sources},duplicate=1)
