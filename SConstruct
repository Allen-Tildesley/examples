# SConstruct: SCons build file for the various programs. See https://scons.org

#------------------------------------------------------------------------------------------------#
# This software was written                                                                      #
# by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>                        #
# and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),                             #
# to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),     #
# published by Oxford University Press ("the publishers").                                       #
#                                                                                                #
# LICENCE                                                                                        #
# Creative Commons CC0 Public Domain Dedication.                                                 #
# To the extent possible under law, the authors have dedicated all copyright and related         #
# and neighboring rights to this software to the PUBLIC domain worldwide.                        #
# This software is distributed without any warranty.                                             #
# You should have received a copy of the CC0 Public Domain Dedication along with this software.  #
# If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.                               #
#                                                                                                #
# DISCLAIMER                                                                                     #
# The authors and publishers make no warranties about the software, and disclaim liability       #
# for all uses of the software, to the fullest extent permitted by applicable law.               #
# The authors and publishers do not recommend use of this software for any purpose.              #
# It is made freely available, solely to clarify points made in the text. When using or citing   #
# the software, you should not imply endorsement by the authors or publishers.                   #
#                                                                                                #
# If you see an error in this file, or a specific improvement that would improve robustness      #
# or portability, please feel free to raise an issue on the GitHub examples repository.          #
#------------------------------------------------------------------------------------------------#

# Build the programs, each in its own subdirectory of the build directory, by typing
# scons
# Clean the programs (i.e., remove object, module and executable files) by typing
# scons -c
# or just remove the build directory and all its subdirectories

import os, sys

# SCons default Fortran compiler (currently) is gfortran, and we assume this in setting flags
# Filetypes are .f90 but the language standard is more modern
# NB by default we do not invoke any optimization
# You could add, for example, -O2 to the F90FLAGS string below

env=Environment(ENV=os.environ.copy())
env.Append(F90FLAGS='-fdefault-real-8 -fall-intrinsics -std=f2018 -Wall')
env.Append(FORTRANMODDIR='${TARGET.dir}',F90PATH='${TARGET.dir}')

# The LAPACK library is required for some of the programs
env_lapack=env.Clone()
env_lapack.ParseConfig("pkg-config lapack --cflags --libs")
# If pkg-config cannot find lapack, look in places such as /usr/lib, /usr/local/lib, /opt/local/lib
# or subdirectories thereof, and do something like
# env_lapack.Append(LIBPATH='/opt/local/lib/lapack',LIBS='lapack')

# The FFTW3 library is required for some of the programs, with an include file
# Workaround applied, see https://github.com/SCons/scons/issues/4835
env_fftw=env.Clone()
env_fftw.ParseConfig("pkg-config fftw3 --cflags --libs")

# If FFTW3 was installed in /usr/include, we need additional flags
if env_fftw.get("CPPPATH") is None:
    env_fftw.ParseConfig("pkg-config fftw3 --cflags --libs --keep-system-cflags")
env_fftw.Append(F90PATH=env_fftw['CPPPATH'])
# If pkg-config cannot find fftw3, look in places such as /usr/lib, /usr/local/lib, /opt/local/lib
# and similar include directories, or subdirectories thereof, and do something like
# env_fftw.Append(F90PATH='/opt/local/include')
# env_fftw.Append(LIBPATH='/opt/local/lib',LIBS='fftw3')

# The OpenMPI library is required for one program, should just require the compiler wrapper
env_mpi=env.Clone(F90='mpifort',LINK='mpifort')

# The OpenMP library is required for one program, but we only need the compiler flags
env_omp=env.Clone()
env_omp.Append(F90FLAGS='-fopenmp',LINKFLAGS='-fopenmp')

# Abbreviations for commonly-used combinations of modules
conavgmat=['config_io_module.f90','averages_module.f90','maths_module.f90']
conavg=['config_io_module.f90','averages_module.f90']
conmat=['config_io_module.f90','maths_module.f90']

# Programs that are compiled in the default env environment
builds={}
builds['adjust']               = ['adjust.f90']+conmat
builds['bd_nvt_lj']            = ['bd_nvt_lj.f90','md_lj_module.f90','lrc_lj_module.f90']+conavgmat
builds['dpd']                  = ['dpd.f90','dpd_module.f90']+conavgmat
builds['cluster']              = ['cluster.f90','config_io_module.f90']
builds['diffusion']            = ['diffusion.f90','config_io_module.f90']
builds['diffusion_test']       = ['diffusion_test.f90']+conmat
builds['eos_hs']               = ['eos_hs.f90','maths_module.f90']
builds['eos_lj']               = ['eos_lj.f90','eos_lj_module.f90','lrc_lj_module.f90']
builds['error_calc']           = ['error_calc.f90','maths_module.f90']
builds['grint']                = ['grint.f90','grint_module.f90','config_io_module.f90']
builds['hit_and_miss']         = ['hit_and_miss.f90']
builds['initialize']           = ['initialize.f90','initialize_module.f90']+conmat
builds['mc_chain_nvt_cbmc_lj'] = ['mc_chain_nvt_cbmc_lj.f90','mc_chain_lj_module.f90']+conavgmat
builds['mc_chain_nvt_sw']      = ['mc_chain_nvt_sw.f90','mc_chain_sw_module.f90']+conavgmat
builds['mc_chain_wl_sw']       = ['mc_chain_wl_sw.f90','mc_chain_sw_module.f90']+conavgmat
builds['mc_gibbs_lj']          = ['mc_gibbs_lj.f90','mc_gibbs_lj_module.f90','lrc_lj_module.f90']+conavgmat
builds['mc_npt_hs']            = ['mc_npt_hs.f90','mc_hs_module.f90']+conavgmat
builds['mc_npt_lj']            = ['mc_npt_lj.f90','mc_lj_module.f90','lrc_lj_module.f90']+conavgmat
builds['mc_npt_lj_ll']         = ['mc_npt_lj.f90','mc_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+conavgmat
builds['mc_npt_sc']            = ['mc_npt_sc.f90','mc_sc_module.f90']+conavgmat
builds['mc_nvt_hs']            = ['mc_nvt_hs.f90','mc_hs_module.f90']+conavgmat
builds['mc_nvt_lj']            = ['mc_nvt_lj.f90','mc_lj_module.f90','lrc_lj_module.f90']+conavgmat
builds['mc_nvt_lj_ll']         = ['mc_nvt_lj.f90','mc_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+conavgmat
builds['mc_nvt_poly_lj']       = ['mc_nvt_poly_lj.f90','mc_poly_lj_module.f90']+conavgmat
builds['mc_nvt_sc']            = ['mc_nvt_sc.f90','mc_sc_module.f90']+conavgmat
builds['mc_zvt_lj']            = ['mc_zvt_lj.f90','mc_lj_module.f90','lrc_lj_module.f90']+conavgmat
builds['mc_zvt_lj_ll']         = ['mc_zvt_lj.f90','mc_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+conavgmat
builds['md_lj_mts']            = ['md_lj_mts.f90','md_lj_mts_module.f90','lrc_lj_module.f90']+conavg
builds['md_npt_lj']            = ['md_npt_lj.f90','md_lj_module.f90','lrc_lj_module.f90']+conavgmat
builds['md_npt_lj_ll']         = ['md_npt_lj.f90','md_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+conavgmat
builds['md_nve_hs']            = ['md_nve_hs.f90','md_nve_hs_module.f90']+conavg
builds['md_nve_lj']            = ['md_nve_lj.f90','md_lj_module.f90','lrc_lj_module.f90']+conavg
builds['md_nve_lj_vl']         = ['md_nve_lj.f90','md_lj_vl_module.f90','lrc_lj_module.f90','verlet_list_module.f90']+conavg
builds['md_nve_lj_ll']         = ['md_nve_lj.f90','md_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+conavg
builds['md_nvt_lj']            = ['md_nvt_lj.f90','md_lj_module.f90','lrc_lj_module.f90']+conavgmat
builds['md_nvt_lj_ll']         = ['md_nvt_lj.f90','md_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90']+conavgmat
builds['md_nvt_lj_le']         = ['md_nvt_lj_le.f90','md_lj_le_module.f90']+conavg
builds['md_nvt_lj_llle']       = ['md_nvt_lj_le.f90','md_lj_llle_module.f90','link_list_module.f90']+conavg
builds['md_nvt_poly_lj']       = ['md_nvt_poly_lj.f90','md_poly_lj_module.f90']+conavgmat
builds['mesh']                 = ['mesh.f90','mesh_module.f90']
builds['pair_distribution']    = ['pair_distribution.f90','config_io_module.f90']
builds['qmc_pi_sho']           = ['qmc_pi_sho.f90','averages_module.f90','maths_module.f90']
builds['qmc_pi_lj']            = ['qmc_pi_lj.f90','qmc_pi_lj_module.f90','lrc_lj_module.f90']+conavgmat
builds['qmc_walk_sho']         = ['qmc_walk_sho.f90','averages_module.f90','maths_module.f90']
builds['sample_mean']          = ['sample_mean.f90']
builds['smc_nvt_lj']           = ['smc_nvt_lj.f90','smc_lj_module.f90','lrc_lj_module.f90']+conavgmat
builds['test_pot_at']          = ['test_pot_atom.f90','test_pot_at.f90','maths_module.f90']
builds['test_pot_bend']        = ['test_pot_atom.f90','test_pot_bend.f90','maths_module.f90']
builds['test_pot_twist']       = ['test_pot_atom.f90','test_pot_twist.f90','maths_module.f90']
builds['test_pot_dd']          = ['test_pot_linear.f90','test_pot_dd.f90','maths_module.f90']
builds['test_pot_dq']          = ['test_pot_linear.f90','test_pot_dq.f90','maths_module.f90']
builds['test_pot_qq']          = ['test_pot_linear.f90','test_pot_qq.f90','maths_module.f90']
builds['test_pot_gb']          = ['test_pot_linear.f90','test_pot_gb.f90','maths_module.f90']
builds['t_tensor']             = ['t_tensor.f90','maths_module.f90']
builds['wl_hist']              = ['wl_hist.f90']

for name, sources in builds.items():
    SConscript('SConscript', variant_dir='build/'+name,exports={'env':env,'name':name,'sources':sources},duplicate=0)

# Programs that are compiled in the env_lapack environment
builds={}
builds['md_chain_mts_lj'] = ['md_chain_mts_lj.f90','md_chain_lj_module.f90']+conavg
builds['md_chain_nve_lj'] = ['md_chain_nve_lj.f90','md_chain_lj_module.f90']+conavgmat

for name, sources in builds.items():
    SConscript('SConscript', variant_dir='build/'+name,exports={'env':env_lapack,'name':name,'sources':sources},duplicate=0)

# Programs that are compiled in the env_fftw environment
builds={}
builds['corfun']    = ['corfun.f90','maths_module.f90']
builds['ewald']     = ['ewald.f90','ewald_module.f90','mesh_module.f90','config_io_module.f90']
builds['fft3dwrap'] = ['fft3dwrap.f90']

for name, sources in builds.items():
    SConscript('SConscript', variant_dir='build/'+name,exports={'env':env_fftw,'name':name,'sources':sources},duplicate=0)

# Program that is compiled in the env_mpi environment
name, sources = 'mc_nvt_lj_re', ['mc_nvt_lj_re.f90','mc_lj_module.f90','lrc_lj_module.f90']+conavgmat
SConscript('SConscript', variant_dir='build/'+name,exports={'env':env_mpi,'name':name,'sources':sources},duplicate=0)

# Program that is compiled in the env_omp environment
name, sources = 'md_nve_lj_omp', ['md_nve_lj.f90','md_lj_omp_module.f90','lrc_lj_module.f90']+conavg
SConscript('SConscript', variant_dir='build/'+name,exports={'env':env_omp,'name':name,'sources':sources},duplicate=0)
