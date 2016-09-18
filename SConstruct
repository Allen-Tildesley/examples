# SCons build file for the various programs
# Run this typing 'scons'
# Clean the program by 'scons -c' (i.e., remove object, module and executable)

import os, sys

env_normal=Environment(ENV=os.environ)
env_lapack=Environment(ENV=os.environ)
env_fftw=Environment(ENV=os.environ)
env_mpi=Environment(ENV=os.environ,F90='mpif90',LINK='mpif90')

# Assume that gfortran will be used. 
#Tool('gfortran')(env_normal)
#Tool('mpif90')(env_mpi)
#MY_F90FLAGS='-O2 -finline-functions -funswitch-loops -fwhole-file'
MY_F90FLAGS='-fdefault-real-8 -fall-intrinsics -std=f2008 -Wall'
FFTW_F90FLAGS='-fdefault-real-8 -fall-intrinsics -std=f2008 -Wall -I/opt/local/include'
MY_LINKFLAGS=''
LAPACK_LINKFLAGS='-L/opt/local/lib/lapack -llapack'
FFTW_LINKFLAGS='-L/opt/local/lib -lfftw3'

env_normal.Append(F90FLAGS=MY_F90FLAGS,LINKFLAGS=MY_LINKFLAGS,FORTRANMODDIRPREFIX='-J',FORTRANMODDIR = '${TARGET.dir}',F90PATH='${TARGET.dir}')
env_lapack.Append(F90FLAGS=MY_F90FLAGS,LINKFLAGS=LAPACK_LINKFLAGS,FORTRANMODDIRPREFIX='-J',FORTRANMODDIR = '${TARGET.dir}',F90PATH='${TARGET.dir}')
env_fftw.Append(F90FLAGS=FFTW_F90FLAGS,LINKFLAGS=FFTW_LINKFLAGS,FORTRANMODDIRPREFIX='-J',FORTRANMODDIR = '${TARGET.dir}',F90PATH='${TARGET.dir}')
env_mpi.Append(F90FLAGS=MY_F90FLAGS,LINKFLAGS=MY_LINKFLAGS,FORTRANMODDIRPREFIX='-J',FORTRANMODDIR = '${TARGET.dir}',F90PATH='${TARGET.dir}')

variants={}
variants['build_bd_nvt_lj']            = (['bd_nvt_lj.f90','md_lj_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_dpd']                  = (['dpd.f90','dpd_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_cluster']              = (['cluster.f90','config_io_module.f90'],env_normal)
variants['build_corfun']               = (['corfun.f90','maths_module.f90'],env_fftw)
variants['build_diffusion']            = (['diffusion.f90','config_io_module.f90'],env_normal)
variants['build_fft3dwrap']            = (['fft3dwrap.f90'],env_fftw)
variants['build_hit_and_miss']         = (['hit_and_miss.f90'],env_normal)
variants['build_initialize']           = (['initialize.f90','initialize_module.f90','config_io_module.f90','maths_module.f90'],env_normal)
variants['build_mc_chain_nvt_cbmc_lj'] = (['mc_chain_nvt_cbmc_lj.f90','mc_chain_lj_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_mc_chain_nvt_sw']      = (['mc_chain_nvt_sw.f90','mc_chain_sw_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_mc_chain_wl_sw']       = (['mc_chain_wl_sw.f90','mc_chain_sw_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_mc_npt_lj']            = (['mc_npt_lj.f90','mc_lj_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_mc_npt_lj_ll']         = (['mc_npt_lj.f90','mc_lj_ll_module.f90','link_list_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_mc_nvt_hs']            = (['mc_nvt_hs.f90','mc_hs_module.f90','config_io_module.f90','averages_module.f90'],env_normal)
variants['build_mc_nvt_lj']            = (['mc_nvt_lj.f90','mc_lj_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_mc_nvt_lj_ll']         = (['mc_nvt_lj.f90','mc_lj_ll_module.f90','link_list_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_mc_nvt_lj_re']         = (['mc_nvt_lj_re.f90','mc_lj_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_mpi)
variants['build_mc_nvt_poly_lj']       = (['mc_nvt_poly_lj.f90','mc_poly_lj_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_mc_nvt_sc']            = (['mc_nvt_sc.f90','mc_sc_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_mc_zvt_lj']            = (['mc_zvt_lj.f90','mc_lj_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_mc_zvt_lj_ll']         = (['mc_zvt_lj.f90','mc_lj_ll_module.f90','link_list_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_md_chain_mts_lj']      = (['md_chain_mts_lj.f90','md_chain_lj_module.f90','config_io_module.f90','averages_module.f90'],env_lapack)
variants['build_md_chain_nve_lj']      = (['md_chain_nve_lj.f90','md_chain_lj_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_lapack)
variants['build_md_lj_mts']            = (['md_lj_mts.f90','md_lj_mts_module.f90','config_io_module.f90','averages_module.f90'],env_normal)
variants['build_md_nve_hs']            = (['md_nve_hs.f90','md_nve_hs_module.f90','config_io_module.f90','averages_module.f90'],env_normal)
variants['build_md_nve_lj']            = (['md_nve_lj.f90','md_lj_module.f90','config_io_module.f90','averages_module.f90'],env_normal)
variants['build_md_nve_lj_vl']         = (['md_nve_lj.f90','md_lj_vl_module.f90','verlet_list_module.f90','config_io_module.f90','averages_module.f90'],env_normal)
variants['build_md_nve_lj_ll']         = (['md_nve_lj.f90','md_lj_ll_module.f90','link_list_module.f90','config_io_module.f90','averages_module.f90'],env_normal)
variants['build_md_nvt_lj_le']         = (['md_nvt_lj_le.f90','md_lj_le_module.f90','config_io_module.f90','averages_module.f90'],env_normal)
variants['build_md_nvt_lj_llle']       = (['md_nvt_lj_le.f90','md_lj_llle_module.f90','link_list_module.f90','config_io_module.f90','averages_module.f90'],env_normal)
variants['build_mesh']                 = (['mesh.f90'],env_normal)
variants['build_pair_distribution']    = (['pair_distribution.f90','config_io_module.f90'],env_normal)
variants['build_qmc_pi_sho']           = (['qmc_pi_sho.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_qmc_pi_lj']            = (['qmc_pi_lj.f90','qmc_pi_lj_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_qmc_walk_sho']         = (['qmc_walk_sho.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_smc_nvt_lj']           = (['smc_nvt_lj.f90','smc_lj_module.f90','config_io_module.f90','averages_module.f90','maths_module.f90'],env_normal)
variants['build_test_pot_at']          = (['test_pot_atom.f90','test_pot_at.f90','maths_module.f90'],env_normal)
variants['build_test_pot_bend']        = (['test_pot_atom.f90','test_pot_bend.f90','maths_module.f90'],env_normal)
variants['build_test_pot_twist']       = (['test_pot_atom.f90','test_pot_twist.f90','maths_module.f90'],env_normal)
variants['build_test_pot_dd']          = (['test_pot_linear.f90','test_pot_dd.f90','maths_module.f90'],env_normal)
variants['build_test_pot_dq']          = (['test_pot_linear.f90','test_pot_dq.f90','maths_module.f90'],env_normal)
variants['build_test_pot_qq']          = (['test_pot_linear.f90','test_pot_qq.f90','maths_module.f90'],env_normal)
variants['build_test_pot_gb']          = (['test_pot_linear.f90','test_pot_gb.f90','maths_module.f90'],env_normal)
variants['build_t_tensor']             = (['t_tensor.f90','maths_module.f90'],env_normal)

# Build each variant in appropriate variant directory
for variant_dir,(sources,env) in variants.iteritems():
    SConscript('SConscript', variant_dir=variant_dir,exports={'env':env,'sources':sources},duplicate=1)
