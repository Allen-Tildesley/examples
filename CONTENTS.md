# Contents
This is a list of the Fortran code examples appearing in the text, in order,
indicating the corresponding online program files.
The online [GUIDE](GUIDE.md) provides more information about running these programs,
but there they are grouped together by system studied,
for easier comparison of the typical test results.
Although links are provided here to the individual files,
it is expected that you will have downloaded them all together.

Not all the examples have online files: in some cases the code was provided in the text.

The utility module files
[averages_module.f90](averages_module.f90),
[config_io_module.f90](config_io_module.f90), and
[maths_module.f90](maths_module.f90),
and a small module containing long-range-correction formulae
[lrc_lj_module.f90](lrc_lj_module.f90)
are widely used by the simulation programs.
They are not explicitly listed in most cases below:
the [SConstruct](SConstruct) file gives a full list of the
build dependencies in each case.

Several other small programs are not described in the text,
but are mentioned in the online [GUIDE](GUIDE.md):
* [adjust.f90](adjust.f90)
* [eos_lj.f90](eos_lj.f90) with [eos_lj_module.f90](eos_lj_module.f90)
* [eos_hs.f90](eos_hs.f90)

An additional molecular dynamics program using quaternions
[md_nvt_poly_lj.f90](md_nvt_poly_lj.f90) and [md_poly_lj_module.f90](md_poly_lj_module.f90)
is described in the online [GUIDE](GUIDE.md),
but does not appear in the text.

### 1.1 Calculation of T tensors
[t_tensor.f90](t_tensor.f90).

### 1.2 Double loop for Lennard-Jones potential
This code snippet appears in the text, not online.
However working examples of a similar kind may be found in
[mc_lj_module.f90](mc_lj_module.f90) and
[md_lj_module.f90](md_lj_module.f90).

### 1.3 Site–site potential energy calculation
This code snippet appears in the text, not online.
However a working example of a similar kind may be found in
[mc_poly_lj_module.f90](mc_poly_lj_module.f90).

### 1.4 Periodic boundaries for truncated octahedron
This code snippet appears in the text, not online.

### 1.5 Periodic boundaries for rhombic dodecahedron
This code snippet appears in the text, not online.

### 1.6 Periodic boundaries for rhombus
This code snippet appears in the text, not online.

### 3.1 Velocity Verlet algorithm
This code snippet appears in the text, not online.
However a working example of a similar kind may be found in
[md_nve_lj.f90](md_nve_lj.f90),
and in several other MD program files.

### 3.2 Velocity Verlet algorithm (original)
This code snippet appears in the text, not online.

### 3.3 Classic Verlet algorithm
This code snippet appears in the text, not online.

### 3.4 Molecular dynamics, NVE-ensemble, Lennard-Jones
[md_nve_lj.f90](md_nve_lj.f90) and [md_lj_module.f90](md_lj_module.f90).

### 3.5 Constraint algorithms for a chain molecule
[md_chain_nve_lj.f90](md_chain_nve_lj.f90) and [md_chain_lj_module.f90](md_chain_lj_module.f90).

### 3.6 Multiple-timestep algorithm
[md_chain_mts_lj.f90](md_chain_mts_lj.f90) and [md_chain_lj_module.f90](md_chain_lj_module.f90).

### 3.7 Molecular dynamics of hard spheres
[md_nve_hs.f90](md_nve_hs.f90) and [md_nve_hs_module.f90](md_nve_hs_module.f90).

### 3.8 Measure-preserving constant-NVT MD algorithm
[md_nvt_lj.f90](md_nvt_lj.f90) and [md_lj_module.f90](md_lj_module.f90).

### 3.9 Measure-preserving constant-NPT MD algorithm
[md_npt_lj.f90](md_npt_lj.f90) and [md_lj_module.f90](md_lj_module.f90).

### 4.1 Hit-and-miss integration
[hit_and_miss.f90](hit_and_miss.f90) (also appears in full in the text).

### 4.2 Sample mean integration
[sample_mean.f90](sample_mean.f90) (also appears in full in the text).

### 4.3 Monte Carlo NVT-ensemble for Lennard-Jones atoms
[mc_nvt_lj.f90](mc_nvt_lj.f90) and [mc_lj_module.f90](mc_lj_module.f90).

### 4.4 Monte Carlo of hard spheres
[mc_nvt_hs.f90](mc_nvt_hs.f90) and [mc_hs_module.f90](mc_hs_module.f90).
We also supply a constant-pressure version [mc_npt_hs.f90](mc_npt_hs.f90).

### 4.5 Monte Carlo NPT-ensemble for Lennard-Jones atoms
[mc_npt_lj.f90](mc_npt_lj.f90) and [mc_lj_module.f90](mc_lj_module.f90).

### 4.6 Monte Carlo &mu;VT-ensemble for Lennard-Jones atoms
[mc_zvt_lj.f90](mc_zvt_lj.f90) and  [mc_lj_module.f90](mc_lj_module.f90).

### 4.7 Monte Carlo program using quaternions
[mc_nvt_poly_lj.f90](mc_nvt_poly_lj.f90) and [mc_poly_lj_module.f90](mc_poly_lj_module.f90).
We also supply a molecular dynamics program
[md_nvt_poly_lj.f90](md_nvt_poly_lj.f90) and [md_poly_lj_module.f90](md_poly_lj_module.f90)
to simulate the same model.

### 4.8 Monte Carlo of hard spherocylinders
[mc_nvt_sc.f90](mc_nvt_sc.f90) and [mc_sc_module.f90](mc_sc_module.f90).
We also supply a constant-pressure version [mc_npt_sc.f90](mc_npt_sc.f90).

### 5.1 Force routine using Verlet neighbour list
[md_lj_vl_module.f90](md_lj_vl_module.f90) and [verlet_list_module.f90](verlet_list_module.f90).
These may be used instead of [md_lj_module.f90](md_lj_module.f90)
in combination with [md_nve_lj.f90](md_nve_lj.f90).

### 5.2 Building a cell structure with linked lists
This code snippet appears in the text, not online.
However a working version of a similar kind may be found in
[link_list_module.f90](link_list_module.f90) of Code 5.3.

### 5.3 Force routine using linked lists
[md_lj_ll_module.f90](md_lj_ll_module.f90) and [link_list_module.f90](link_list_module.f90).
These may be used instead of [md_lj_module.f90](md_lj_module.f90)
in combination with [md_nve_lj.f90](md_nve_lj.f90) or [md_nvt_lj.f90](md_nvt_lj.f90).

### 5.4 Monte Carlo routines using linked lists
[mc_lj_ll_module.f90](mc_lj_ll_module.f90) and [link_list_module.f90](link_list_module.f90).
These may be used instead of [mc_lj_module.f90](mc_lj_module.f90)
in combination with [mc_nvt_lj.f90](mc_nvt_lj.f90),
[mc_npt_lj.f90](mc_npt_lj.f90), or
[mc_zvt_lj.f90](mc_zvt_lj.f90).
However, for simplicity the routines have not been designed to be especially robust against
large changes in system dimensions or number of particles: this could be fixed
by including some memory reallocation statements in the Fortran.

### 5.5 Multiple-timestep algorithm, Lennard-Jones atoms
[md_lj_mts.f90](md_lj_mts.f90) and [md_lj_mts_module.f90](md_lj_mts_module.f90).

### 5.6 Initialization of a crystal lattice
[initialize.f90](initialize.f90) and [initialize_module.f90](initialize_module.f90).
This program also has options to initialize a random configuration, and a chain of atoms.

### 6.1 Force routine using the Ewald sum
[ewald.f90](ewald.f90) with [ewald_module.f90](ewald_module.f90) and [mesh_module.f90](mesh_module.f90).
As well as the real-space and reciprocal-space routines mentioned in the text,
we also illustrate the particle-mesh-Ewald method.

### 6.2 Assignment of charge to a uniform mesh
[mesh.f90](mesh.f90) and [mesh_module.f90](mesh_module.f90).

### 7.1 Parallelized double-loop, shared memory
This code snippet appears in the text, not online.
However a working example of a similar kind may be found in
[md_lj_omp_module.f90](md_lj_omp_module.f90) of Code 7.2.

### 7.2 Parallelized force routine, shared memory
[md_lj_omp_module.f90](md_lj_omp_module.f90).
This may be used instead of [md_lj_module.f90](md_lj_module.f90)
in combination with [md_nve_lj.f90](md_nve_lj.f90).

### 7.3 Replica exchange, by message-passing
[mc_nvt_lj_re.f90](mc_nvt_lj_re.f90) and [mc_lj_module.f90](mc_lj_module.f90).

### 8.1 Calculating the pair distribution function
[pair_distribution.f90](pair_distribution.f90).
Part of the code also appears in the text.

### 8.2 Calculating a time correlation function
This code snippet appears in the text, not online.
However working examples of a similar kind may be found in
[diffusion.f90](diffusion.f90) of Code 8.3 and
[corfun.f90](corfun.f90) of Code 8.4.

### 8.3 Program to compute diffusion coefficients
[diffusion.f90](diffusion.f90).
Also supplied is a program [diffusion_test.f90](diffusion_test.f90) to generate
test data.

### 8.4 Calculating time correlation functions
[corfun.f90](corfun.f90).

### 8.5 Calculating statistical inefficiency and errors
[error_calc.f90](error_calc.f90).

### 9.1 Wang–Landau simulation of chain molecule
[mc_chain_wl_sw.f90](mc_chain_wl_sw.f90) and [mc_chain_sw_module.f90](mc_chain_sw_module.f90).
Also provided is a small program
[wl_hist.f90](wl_hist.f90)
to do post-processing of the histograms produced by the Wang-Landau simulation
and, for comparison,
a constant-temperature MC program
[mc_chain_nvt_sw.f90](mc_chain_nvt_sw.f90).

### 9.2 Configuration-biased simulation of chain molecule
[mc_chain_nvt_cbmc_lj.f90](mc_chain_nvt_cbmc_lj.f90) and [mc_chain_lj_module.f90](mc_chain_lj_module.f90).

### 9.3 Gibbs ensemble simulation
[mc_gibbs_lj.f90](mc_gibbs_lj.f90) and [mc_gibbs_lj_module.f90](mc_gibbs_lj_module.f90).

### 11.1 Molecular dynamics using Lees–Edwards boundaries
[md_nvt_lj_le.f90](md_nvt_lj_le.f90) and [md_lj_le_module.f90](md_lj_le_module.f90).

### 11.2 Cell structure and linked lists in sheared boundaries
[md_lj_llle_module.f90](md_lj_llle_module.f90) and [link_list_module.f90](link_list_module.f90).
These may be used instead of [md_lj_le_module.f90](md_lj_le_module.f90)
in combination with [md_nvt_lj_le.f90](md_nvt_lj_le.f90).

### 12.1 Brownian dynamics program
[bd_nvt_lj.f90](bd_nvt_lj.f90) and [md_lj_module.f90](md_lj_module.f90).

### 12.2 Smart Monte Carlo simulation
[smc_nvt_lj.f90](smc_nvt_lj.f90) and [smc_lj_module.f90](smc_lj_module.f90).

### 12.3 Dissipative particle dynamics
[dpd.f90](dpd.f90) and [dpd_module.f90](dpd_module.f90).

### 13.1 Path-integral Monte Carlo
[qmc_pi_sho.f90](qmc_pi_sho.f90) for the simple harmonic oscillator.
[qmc_pi_lj.f90](qmc_pi_lj.f90) and [qmc_pi_lj_module.f90](qmc_pi_lj_module.f90) for Lennard-Jones.

### 13.2 Quantum Monte Carlo using a diffusion equation
[qmc_walk_sho.f90](qmc_walk_sho.f90).

### 14.1 Calculating the temperature and density profiles
This code snippet appears in the text, not online.
However, a working example of the density profile calculation may be found in
[grint.f90](grint.f90) of Code 14.2.

### 14.2 Radial distribution function in a planar interface
[grint.f90](grint.f90) and [grint_module.f90](grint_module.f90).
We also provide a (comparatively large) file `grint_data.zip` containing test data
in the [Data repository](https://github.com/Allen-Tildesley/data).

### 14.3 Cluster analysis
[cluster.f90](cluster.f90).
We also provide a test configuration file [cluster.inp](cluster.inp).

### A.1 Utility modules
[averages_module.f90](averages_module.f90),
[config_io_module.f90](config_io_module.f90), and
[maths_module.f90](maths_module.f90),
used by many of the programs listed above.

### C.1 Programs to test forces and torques
[test_pot_atom.f90](test_pot_atom.f90) together with any one of
[test_pot_at.f90](test_pot_at.f90),
[test_pot_bend.f90](test_pot_bend.f90), or
[test_pot_twist.f90](test_pot_twist.f90).

[test_pot_linear.f90](test_pot_linear.f90) together with any one of
[test_pot_dd.f90](test_pot_dd.f90),
[test_pot_dq.f90](test_pot_dq.f90),
[test_pot_qq.f90](test_pot_qq.f90), or
[test_pot_gb.f90](test_pot_gb.f90).

### D.1 3D Fourier transform example
[fft3dwrap.f90](fft3dwrap.f90).

### F.1 Calculating the configurational temperature
This code snippet appears in the text, not online.
However, a working example may be found in
[mc_nvt_lj.f90](mc_nvt_lj.f90) and its module [mc_lj_module.f90](mc_lj_module.f90).

### F.2 Correction to configurational temperature
This code snippet appears in the text, not online.
However, a working example may be found in
[md_nve_lj.f90](md_nve_lj.f90) and its module [md_lj_module.f90](md_lj_module.f90).
