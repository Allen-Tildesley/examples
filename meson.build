# meson.build: Meson build file for the various programs. See https://mesonbuild.com

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
# meson setup build
# followed by
# meson compile -C build
# Clean the programs (i.e., remove object, module and executable files) by typing
# meson compile --clean -C build
# or just remove the build directory and all its subdirectories

# We assume that the Fortran compiler is gfortran in setting some flags
# Filetypes are .f90 but the language standard is more modern
# NB by default we do not invoke any optimization
# You could add, for example, 'optimization=2' to the list of default_options below
# or alternatively use the command
# meson setup --optimization 2 build
# instead of meson setup build

project('build_all', 'fortran', meson_version: '>=1.10',default_options:['fortran_std=f2018','warning_level=1'])
add_project_arguments('-fdefault-real-8','-fall-intrinsics',language:'fortran')

# Abbreviations for commonly-used combinations of modules
conavgmat=['config_io_module.f90','averages_module.f90','maths_module.f90']
conavg=['config_io_module.f90','averages_module.f90']
conmat=['config_io_module.f90','maths_module.f90']

# Build executables that don't link with external libraries
builds = {
'adjust'               : ['adjust.f90',conmat],
'bd_nvt_lj'            : ['bd_nvt_lj.f90','md_lj_module.f90','lrc_lj_module.f90',conavgmat],
'dpd'                  : ['dpd.f90','dpd_module.f90',conavgmat],
'cluster'              : ['cluster.f90','config_io_module.f90'],
'diffusion'            : ['diffusion.f90','config_io_module.f90'],
'diffusion_test'       : ['diffusion_test.f90',conmat],
'eos_hs'               : ['eos_hs.f90','maths_module.f90'],
'eos_lj'               : ['eos_lj.f90','eos_lj_module.f90','lrc_lj_module.f90'],
'error_calc'           : ['error_calc.f90','maths_module.f90'],
'grint'                : ['grint.f90','grint_module.f90','config_io_module.f90'],
'hit_and_miss'         : ['hit_and_miss.f90'],
'initialize'           : ['initialize.f90','initialize_module.f90',conmat],
'mc_chain_nvt_cbmc_lj' : ['mc_chain_nvt_cbmc_lj.f90','mc_chain_lj_module.f90',conavgmat],
'mc_chain_nvt_sw'      : ['mc_chain_nvt_sw.f90','mc_chain_sw_module.f90',conavgmat],
'mc_chain_wl_sw'       : ['mc_chain_wl_sw.f90','mc_chain_sw_module.f90',conavgmat],
'mc_gibbs_lj'          : ['mc_gibbs_lj.f90','mc_gibbs_lj_module.f90','lrc_lj_module.f90',conavgmat],
'mc_npt_hs'            : ['mc_npt_hs.f90','mc_hs_module.f90',conavgmat],
'mc_npt_lj'            : ['mc_npt_lj.f90','mc_lj_module.f90','lrc_lj_module.f90',conavgmat],
'mc_npt_lj_ll'         : ['mc_npt_lj.f90','mc_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90',conavgmat],
'mc_npt_sc'            : ['mc_npt_sc.f90','mc_sc_module.f90',conavgmat],
'mc_nvt_hs'            : ['mc_nvt_hs.f90','mc_hs_module.f90',conavgmat],
'mc_nvt_lj'            : ['mc_nvt_lj.f90','mc_lj_module.f90','lrc_lj_module.f90',conavgmat],
'mc_nvt_lj_ll'         : ['mc_nvt_lj.f90','mc_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90',conavgmat],
'mc_nvt_poly_lj'       : ['mc_nvt_poly_lj.f90','mc_poly_lj_module.f90',conavgmat],
'mc_nvt_sc'            : ['mc_nvt_sc.f90','mc_sc_module.f90',conavgmat],
'mc_zvt_lj'            : ['mc_zvt_lj.f90','mc_lj_module.f90','lrc_lj_module.f90',conavgmat],
'mc_zvt_lj_ll'         : ['mc_zvt_lj.f90','mc_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90',conavgmat],
'md_lj_mts'            : ['md_lj_mts.f90','md_lj_mts_module.f90','lrc_lj_module.f90',conavg],
'md_npt_lj'            : ['md_npt_lj.f90','md_lj_module.f90','lrc_lj_module.f90',conavgmat],
'md_npt_lj_ll'         : ['md_npt_lj.f90','md_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90',conavgmat],
'md_nve_hs'            : ['md_nve_hs.f90','md_nve_hs_module.f90',conavg],
'md_nve_lj'            : ['md_nve_lj.f90','md_lj_module.f90','lrc_lj_module.f90',conavg],
'md_nve_lj_vl'         : ['md_nve_lj.f90','md_lj_vl_module.f90','lrc_lj_module.f90','verlet_list_module.f90',conavg],
'md_nve_lj_ll'         : ['md_nve_lj.f90','md_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90',conavg],
'md_nvt_lj'            : ['md_nvt_lj.f90','md_lj_module.f90','lrc_lj_module.f90',conavgmat],
'md_nvt_lj_ll'         : ['md_nvt_lj.f90','md_lj_ll_module.f90','lrc_lj_module.f90','link_list_module.f90',conavgmat],
'md_nvt_lj_le'         : ['md_nvt_lj_le.f90','md_lj_le_module.f90',conavg],
'md_nvt_lj_llle'       : ['md_nvt_lj_le.f90','md_lj_llle_module.f90','link_list_module.f90',conavg],
'md_nvt_poly_lj'       : ['md_nvt_poly_lj.f90','md_poly_lj_module.f90',conavgmat],
'mesh'                 : ['mesh.f90','mesh_module.f90'],
'pair_distribution'    : ['pair_distribution.f90','config_io_module.f90'],
'qmc_pi_sho'           : ['qmc_pi_sho.f90','averages_module.f90','maths_module.f90'],
'qmc_pi_lj'            : ['qmc_pi_lj.f90','qmc_pi_lj_module.f90','lrc_lj_module.f90',conavgmat],
'qmc_walk_sho'         : ['qmc_walk_sho.f90','averages_module.f90','maths_module.f90'],
'sample_mean'          : ['sample_mean.f90'],
'smc_nvt_lj'           : ['smc_nvt_lj.f90','smc_lj_module.f90','lrc_lj_module.f90',conavgmat],
'test_pot_at'          : ['test_pot_atom.f90','test_pot_at.f90','maths_module.f90'],
'test_pot_bend'        : ['test_pot_atom.f90','test_pot_bend.f90','maths_module.f90'],
'test_pot_twist'       : ['test_pot_atom.f90','test_pot_twist.f90','maths_module.f90'],
'test_pot_dd'          : ['test_pot_linear.f90','test_pot_dd.f90','maths_module.f90'],
'test_pot_dq'          : ['test_pot_linear.f90','test_pot_dq.f90','maths_module.f90'],
'test_pot_qq'          : ['test_pot_linear.f90','test_pot_qq.f90','maths_module.f90'],
'test_pot_gb'          : ['test_pot_linear.f90','test_pot_gb.f90','maths_module.f90'],
't_tensor'             : ['t_tensor.f90','maths_module.f90'],
'wl_hist'              : ['wl_hist.f90']
}
foreach name, source : builds
	executable(name,source,build_subdir:name)
endforeach

# Build executables that link with the lapack library
lapackdep = dependency('lapack')
# If meson cannot find lapack, look in places such as /usr/lib, /usr/local/lib, /opt/local/lib
# or subdirectories thereof, and do something like
# lapackdep = declare_dependency(link_args:['-L/opt/local/lib/lapack','-llapack'])

builds = {
'md_chain_mts_lj' : ['md_chain_mts_lj.f90','md_chain_lj_module.f90',conavg],
'md_chain_nve_lj' : ['md_chain_nve_lj.f90','md_chain_lj_module.f90',conavgmat]
}
foreach name, source : builds
	executable(name,source,build_subdir:name,dependencies:lapackdep)
endforeach

# Build executables that link with the fftw3 library
fftw3dep = dependency('fftw3')
# If meson cannot find fftw3, look in places such as /usr/lib, /usr/local/lib, /opt/local/lib
# and similar include directories, or subdirectories thereof, and do something like
# fftw3dep = declare_dependency(link_args:['-L/opt/local/lib','-lfftw3'],include_directories:['/opt/local/include'])
builds = {
'corfun'    : ['corfun.f90','maths_module.f90'],
'ewald'     : ['ewald.f90','ewald_module.f90','mesh_module.f90','config_io_module.f90'],
'fft3dwrap' : ['fft3dwrap.f90']
}
foreach name, source : builds
	executable(name,source,build_subdir:name,dependencies:fftw3dep)
endforeach

# Build MPI executable (should find compiler wrapper, e.g. mpifort)
mpidep = dependency('mpi', language: 'fortran')
name   = 'mc_nvt_lj_re'
source = ['mc_nvt_lj_re.f90','mc_lj_module.f90','lrc_lj_module.f90',conavgmat]
executable(name,source,build_subdir:name,dependencies:mpidep)

# Build OpenMP executable (should just need compiler and linker options -fopenmp)
name   = 'md_nve_lj_omp'
source = ['md_nve_lj.f90','md_lj_omp_module.f90','lrc_lj_module.f90',conavg]
executable(name,source,build_subdir:name,fortran_args:'-fopenmp',link_args:'-fopenmp')
