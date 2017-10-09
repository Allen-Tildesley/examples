# SCons build file for the various programs
# Run this by typing 'scons'
# Clean the programs by 'scons -c' (i.e., remove object, module and executable files)

import os, sys

# This has been tested using SCons v2.5.1, gfortran v6.3,
# using MacOS Sierra (10.12) with compilers and libraries installed through MacPorts.
# It may not work on your system. It is possible that you can get it to work by
# changing the flags and library/include paths defined in the following few statements.
# The most likely trouble spots are the programs that use the non-standard
# environments: env_lapack, env_fftw, env_mpi, and env_omp. You may be able to make
# some progress by commenting out the corresponding "variants" lines below, and
# compiling only the programs that use env_normal.

# If you don't like using SCons, or can't get it to work,
# it is not difficult to compile the programs using other methods.
# Bear in mind that, with Fortran, it is usually essential to compile any
# modules that are used by the main program, before compiling the main program itself.
# Take a look at this file in any case, as it shows the file dependencies for each example.

# If you see an error in this file, or a specific improvement that would improve robustness
# or portability, please feel free to raise an issue on the GitHub examples repository.
# Unfortunately, due to the enormous variety of computing platforms and compilers,
# we cannot offer more specific advice on the build process.

# Assume that gfortran will be used as the default compiler
# NB by default we do not invoke any optimization
MY_FLAGS='-fdefault-real-8 -fall-intrinsics -std=f2008 -Wall'
LAPACK_LIBPATH='/opt/local/lib/lapack'
LAPACK_LIBS='lapack'
FFTW_LIBPATH='/opt/local/lib'
FFTW_LIBS='fftw3'
FFTW_INCLUDE='/opt/local/include'
OMP_FLAGS='-fopenmp'
OMP_LINKFLAGS='-fopenmp'

env_normal=Environment(ENV=os.environ)
env_normal.Append(F90FLAGS=MY_FLAGS,FORTRANMODDIRPREFIX='-J',FORTRANMODDIR='${TARGET.dir}',F90PATH=['${TARGET.dir}'])

env_lapack=env_normal.Clone()
env_lapack.Append(LIBPATH=[LAPACK_LIBPATH],LIBS=LAPACK_LIBS)

env_fftw=env_normal.Clone()
env_fftw.Append(F90PATH=[FFTW_INCLUDE])
env_fftw.Append(LIBPATH=[FFTW_LIBPATH],LIBS=FFTW_LIBS)

env_mpi=env_normal.Clone(F90='mpif90',LINK='mpif90')

env_omp=env_normal.Clone()
env_omp.Append(F90FLAGS=OMP_FLAGS,LINKFLAGS=OMP_LINKFLAGS)

utilities=['config_io_module.f90','averages_module.f90','maths_module.f90']
utnomaths=['config_io_module.f90','averages_module.f90']
utnoavrgs=['config_io_module.f90','maths_module.f90']
variants={}
variants['build_adjust']               = (['adjust.f90']+utnoavrgs,env_normal)
variants['build_bd_nvt_lj']            = (['bd_nvt_lj.f90','md_lj_module.f90','lrc_lj_module.f90']+utilities,env_normal)
variants['build_dpd']                  = (['dpd.f90','dpd_module.f90']+utilities,env_normal)
variants['build_cluster']              = (['cluster.f90','config_io_module.f90'],env_normal)
variants['build_corfun']               = (['corfun.f90','maths_module.f90'],env_fftw)
variants['build_diffusion']            = (['diffusion.f90','config_io_module.f90'],env_normal)
variants['build_diffusion_test']       = (['diffusion_test.f90']+utnoavrgs,env_normal)
variants['build_eos_lj']               = (['eos_lj.f90','eos_lj_module.f90','lrc_lj_module.f90'],env_normal)
variants['build_eos_hs']               = (['eos_hs.f90'],env_normal)
variants['build_error_calc']           = (['error_calc.f90','maths_module.f90'],env_normal)
variants['build_ewald']                = (['ewald.f90','ewald_module.f90','mesh_module.f90','config_io_module.f90'],env_fftw)
variants['build_fft3dwrap']            = (['fft3dwrap.f90'],env_fftw)
variants['build_grint']                = (['grint.f90','grint_module.f90','config_io_module.f90'],env_normal)
variants['build_hit_and_miss']         = (['hit_and_miss.f90'],env_normal)
variants['build_initialize']           = (['initialize.f90','initialize_module.f90']+utnoavrgs,env_normal)
variants['build_mc_chain_nvt_cbmc_lj'] = (['mc_chain_nvt_cbmc_lj.f90','mc_chain_lj_module.f90']+utilities,env_normal)
variants['build_mc_chain_nvt_sw']      = (['mc_chain_nvt_sw.f90','mc_chain_sw_module.f90']+utilities,env_normal)
variants['build_mc_chain_wl_sw']       = (['mc_chain_wl_sw.f90','mc_chain_sw_module.f90']+utilities,env_normal)
variants['build_mc_gibbs_lj']          = (['mc_gibbs_lj.f90','mc_gibbs_lj_module.f90','lrc_lj_module.f90']+utilities,env_normal)
variants['build_mc_npt_hs']            = (['mc_npt_hs.f90','mc_hs_module.f90']+utilities,env_normal)
variants['build_mc_npt_lj']            = (['mc_npt_lj.f90','mc_lj_module.f90','lrc_lj_module.f90']+utilities,env_normal)
variants['build_mc_npt_lj_ll']         = (['mc_npt_lj.f90','mc_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+utilities,env_normal)
variants['build_mc_npt_sc']            = (['mc_npt_sc.f90','mc_sc_module.f90']+utilities,env_normal)
variants['build_mc_nvt_hs']            = (['mc_nvt_hs.f90','mc_hs_module.f90']+utilities,env_normal)
variants['build_mc_nvt_lj']            = (['mc_nvt_lj.f90','mc_lj_module.f90','lrc_lj_module.f90']+utilities,env_normal)
variants['build_mc_nvt_lj_ll']         = (['mc_nvt_lj.f90','mc_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+utilities,env_normal)
variants['build_mc_nvt_lj_re']         = (['mc_nvt_lj_re.f90','mc_lj_module.f90','lrc_lj_module.f90']+utilities,env_mpi)
variants['build_mc_nvt_poly_lj']       = (['mc_nvt_poly_lj.f90','mc_poly_lj_module.f90']+utilities,env_normal)
variants['build_mc_nvt_sc']            = (['mc_nvt_sc.f90','mc_sc_module.f90']+utilities,env_normal)
variants['build_mc_zvt_lj']            = (['mc_zvt_lj.f90','mc_lj_module.f90','lrc_lj_module.f90']+utilities,env_normal)
variants['build_mc_zvt_lj_ll']         = (['mc_zvt_lj.f90','mc_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+utilities,env_normal)
variants['build_md_chain_mts_lj']      = (['md_chain_mts_lj.f90','md_chain_lj_module.f90']+utnomaths,env_lapack)
variants['build_md_chain_nve_lj']      = (['md_chain_nve_lj.f90','md_chain_lj_module.f90']+utilities,env_lapack)
variants['build_md_lj_mts']            = (['md_lj_mts.f90','md_lj_mts_module.f90','lrc_lj_module.f90']+utnomaths,env_normal)
variants['build_md_npt_lj']            = (['md_npt_lj.f90','md_lj_module.f90','lrc_lj_module.f90']+utilities,env_normal)
variants['build_md_nve_hs']            = (['md_nve_hs.f90','md_nve_hs_module.f90']+utnomaths,env_normal)
variants['build_md_nve_lj']            = (['md_nve_lj.f90','md_lj_module.f90','lrc_lj_module.f90']+utnomaths,env_normal)
variants['build_md_nve_lj_vl']         = (['md_nve_lj.f90','md_lj_vl_module.f90','lrc_lj_module.f90','verlet_list_module.f90']+utnomaths,env_normal)
variants['build_md_nve_lj_ll']         = (['md_nve_lj.f90','md_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+utnomaths,env_normal)
variants['build_md_nve_lj_omp']        = (['md_nve_lj.f90','md_lj_omp_module.f90','lrc_lj_module.f90']+utnomaths,env_omp)
variants['build_md_nvt_lj']            = (['md_nvt_lj.f90','md_lj_module.f90','lrc_lj_module.f90']+utilities,env_normal)
variants['build_md_nvt_lj_ll']         = (['md_nvt_lj.f90','md_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+utilities,env_normal)
variants['build_md_nvt_lj_le']         = (['md_nvt_lj_le.f90','md_lj_le_module.f90']+utnomaths,env_normal)
variants['build_md_nvt_lj_llle']       = (['md_nvt_lj_le.f90','md_lj_llle_module.f90','link_list_module.f90']+utnomaths,env_normal)
variants['build_md_nvt_poly_lj']       = (['md_nvt_poly_lj.f90','md_poly_lj_module.f90']+utilities,env_normal)
variants['build_mesh']                 = (['mesh.f90','mesh_module.f90'],env_normal)
variants['build_pair_distribution']    = (['pair_distribution.f90','config_io_module.f90'],env_normal)
variants['build_qmc_pi_sho']           = (['qmc_pi_sho.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_qmc_pi_lj']            = (['qmc_pi_lj.f90','qmc_pi_lj_module.f90','lrc_lj_module.f90']+utilities,env_normal)
variants['build_qmc_walk_sho']         = (['qmc_walk_sho.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_sample_mean']          = (['sample_mean.f90'],env_normal)
variants['build_smc_nvt_lj']           = (['smc_nvt_lj.f90','smc_lj_module.f90','lrc_lj_module.f90']+utilities,env_normal)
variants['build_test_pot_at']          = (['test_pot_atom.f90','test_pot_at.f90','maths_module.f90'],env_normal)
variants['build_test_pot_bend']        = (['test_pot_atom.f90','test_pot_bend.f90','maths_module.f90'],env_normal)
variants['build_test_pot_twist']       = (['test_pot_atom.f90','test_pot_twist.f90','maths_module.f90'],env_normal)
variants['build_test_pot_dd']          = (['test_pot_linear.f90','test_pot_dd.f90','maths_module.f90'],env_normal)
variants['build_test_pot_dq']          = (['test_pot_linear.f90','test_pot_dq.f90','maths_module.f90'],env_normal)
variants['build_test_pot_qq']          = (['test_pot_linear.f90','test_pot_qq.f90','maths_module.f90'],env_normal)
variants['build_test_pot_gb']          = (['test_pot_linear.f90','test_pot_gb.f90','maths_module.f90'],env_normal)
variants['build_t_tensor']             = (['t_tensor.f90','maths_module.f90'],env_normal)
variants['build_wl_hist']              = (['wl_hist.f90'],env_normal)

# Build each variant in appropriate variant directory
for variant_dir,(sources,env) in variants.iteritems():
    SConscript('SConscript', variant_dir=variant_dir,exports={'env':env,'sources':sources},duplicate=1)
