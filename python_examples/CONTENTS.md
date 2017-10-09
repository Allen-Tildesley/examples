# Contents
This is a list of Python versions of the code examples appearing in the text, in order,
indicating the corresponding online program files.
The online [GUIDE](GUIDE.md) provides more information about running these programs,
but there they are grouped together by system studied,
for easier comparison of the typical test results.
Although links are provided here to the individual files,
it is expected that you will have downloaded them all together.

Not all the examples have online files: in some cases the code was provided in the text.
We have only provided a *selection* of examples in Python,
whereas the full set is provided in Fortran.

The utility module files
[averages_module.py](averages_module.py),
[config_io_module.py](config_io_module.py), and
[maths_module.py](maths_module.py),
and a small module containing long-range-correction formulae
[lrc_lj_module.py](lrc_lj_module.py)
are widely used by the simulation programs.
They are not explicitly listed in most cases below.

Several other small programs are not described in the text,
but are mentioned in the online [GUIDE](GUIDE.md):
* [adjust.py](adjust.py)
* [eos_lj.py](eos_lj.py) with [eos_lj_module.py](eos_lj_module.py)
* [eos_hs.py](eos_hs.py)

An additional molecular dynamics program using quaternions
[md_nvt_poly_lj.py](md_nvt_poly_lj.py) and [md_poly_lj_module.py](md_poly_lj_module.py)
is described in the online [GUIDE](GUIDE.md),
but does not appear in the text.

### 1.1 Calculation of T tensors
[t_tensor.py](t_tensor.py).

### 1.2 Double loop for Lennard-Jones potential
A Fortran version of this code snippet appears in the text, not online.
However working Python examples of a similar kind may be found in
[mc_lj_module.py](mc_lj_module.py) and
[md_lj_module.py](md_lj_module.py).

### 1.3 Site–site potential energy calculation
A Fortran version of this code snippet appears in the text, not online.

### 1.4 Periodic boundaries for truncated octahedron
A Fortran version of this code snippet appears in the text, not online.

### 1.5 Periodic boundaries for rhombic dodecahedron
A Fortran version of this code snippet appears in the text, not online.

### 1.6 Periodic boundaries for rhombus
A Fortran version of this code snippet appears in the text, not online.

### 3.1 Velocity Verlet algorithm
A Fortran version of this code snippet appears in the text, not online.
However a working Python example of a similar kind may be found in
[md_nve_lj.py](md_nve_lj.py).

### 3.2 Velocity Verlet algorithm (original)
A Fortran version of this code snippet appears in the text, not online.

### 3.3 Classic Verlet algorithm
A Fortran version of this code snippet appears in the text, not online.

### 3.4 Molecular dynamics, NVE-ensemble, Lennard-Jones
[md_nve_lj.py](md_nve_lj.py) and [md_lj_module.py](md_lj_module.py).

### 3.5 Constraint algorithms for a chain molecule
[md_chain_nve_lj.py](md_chain_nve_lj.py) and [md_chain_lj_module.py](md_chain_lj_module.py).

### 3.6 Multiple-timestep algorithm
[md_chain_mts_lj.py](md_chain_mts_lj.py) and [md_chain_lj_module.py](md_chain_lj_module.py).

### 3.7 Molecular dynamics of hard spheres
[md_nve_hs.py](md_nve_hs.py) and [md_nve_hs_module.py](md_nve_hs_module.py).

### 3.8 Measure-preserving constant-NVT MD algorithm
[md_nvt_lj.py](md_nvt_lj.py) and [md_lj_module.py](md_lj_module.py).

### 3.9 Measure-preserving constant-NPT MD algorithm
[md_npt_lj.py](md_npt_lj.py) and [md_lj_module.py](md_lj_module.py).

### 4.1 Hit-and-miss integration
[hit_and_miss.py](hit_and_miss.py)
(a Fortran version also appears in full in the text).

### 4.2 Sample mean integration
[sample_mean.py](sample_mean.py)
(a Fortran version also appears in full in the text).

### 4.3 Monte Carlo NVT-ensemble for Lennard-Jones atoms
[mc_nvt_lj.py](mc_nvt_lj.py) and [mc_lj_module.py](mc_lj_module.py).

### 4.4 Monte Carlo of hard spheres
[mc_nvt_hs.py](mc_nvt_hs.py) and [mc_hs_module.py](mc_hs_module.py).
We also supply a constant-pressure version [mc_npt_hs.py](mc_npt_hs.py).

### 4.5 Monte Carlo NPT-ensemble for Lennard-Jones atoms
[mc_npt_lj.py](mc_npt_lj.py) and [mc_lj_module.py](mc_lj_module.py).

### 4.6 Monte Carlo &mu;VT-ensemble for Lennard-Jones atoms
[mc_zvt_lj.py](mc_zvt_lj.py) and  [mc_lj_module.py](mc_lj_module.py).

### 4.7 Monte Carlo program using quaternions
[mc_nvt_poly_lj.py](mc_nvt_poly_lj.py) and [mc_poly_lj_module.py](mc_poly_lj_module.py).
We also supply a molecular dynamics program
[md_nvt_poly_lj.py](md_nvt_poly_lj.py) and [md_poly_lj_module.py](md_poly_lj_module.py)
to simulate the same model.

### 4.8 Monte Carlo of hard spherocylinders
[mc_nvt_sc.py](mc_nvt_sc.py) and [mc_sc_module.py](mc_sc_module.py).
We also supply a constant-pressure version [mc_npt_sc.py](mc_npt_sc.py).

### 5.1 Force routine using Verlet neighbour list
We do not provide a Python version of this example.
A Fortran version is available.

### 5.2 Building a cell structure with linked lists
A Fortran version of this code snippet appears in the text, not online.

### 5.3 Force routine using linked lists
[md_lj_ll_module.py](md_lj_ll_module.py).
This may be used
in combination with [md_nve_lj.py](md_nve_lj.py), [md_nvt_lj.py](md_nvt_lj.py), or [md_npt_lj.py](md_npt_lj.py).
In each case, the program needs to be edited so as to import routines from this module
instead of [md_lj_module.py](md_lj_module.py).
Despite the name, linked lists are _not_ used, as there is no advantage in Python.
However, the module illustrates the cell structure, and the calculation of forces between neighbouring cells,
analogous to the Fortran example.

### 5.4 Monte Carlo routines using linked lists
We do not provide a Python version of this example.
A Fortran version is available.

### 5.5 Multiple-timestep algorithm, Lennard-Jones atoms
We do not provide a Python version of this example.
A Fortran version is available.

### 5.6 Initialization of a crystal lattice
[initialize.py](initialize.py).
This program also has options to initialize a random configuration, and a chain of atoms.

### 6.1 Force routine using the Ewald sum
[ewald.py](ewald.py) with [ewald_module.py](ewald_module.py) and [mesh_module.py](mesh_module.py).
As well as the real-space and reciprocal-space routines mentioned in the text,
we also illustrate the particle-mesh-Ewald method.

### 6.2 Assignment of charge to a uniform mesh
[mesh.py](mesh.py) and [mesh_module.py](mesh_module.py).

### 7.1 Parallelized double-loop, shared memory
A Fortran version of this code snippet appears in the text, not online.

### 7.2 Parallelized force routine, shared memory
We do not provide a Python version of this example.
A Fortran version is available.

### 7.3 Replica exchange, by message-passing
[mc_nvt_lj_re.py](mc_nvt_lj_re.py) and [mc_lj_module.py](mc_lj_module.py).

### 8.1 Calculating the pair distribution function
[pair_distribution.py](pair_distribution.py).
A Fortran version of part of the code also appears in the text.

### 8.2 Calculating a time correlation function
A Fortran version of this code snippet appears in the text, not online.
However working Python examples of a similar kind may be found in
[diffusion.py](diffusion.py) of Code 8.3 and
[corfun.py](corfun.py) of Code 8.4.

### 8.3 Program to compute diffusion coefficients
[diffusion.py](diffusion.py).
Also supplied is a program [diffusion_test.py](diffusion_test.py) to generate
test data.

### 8.4 Calculating time correlation functions
[corfun.py](corfun.py).

### 8.5 Calculating statistical inefficiency and errors
[error_calc.py](error_calc.py).

### 9.1 Wang–Landau simulation of chain molecule
[mc_chain_wl_sw.py](mc_chain_wl_sw.py) and [mc_chain_sw_module.py](mc_chain_sw_module.py).
Also provided is a small program
[wl_hist.py](wl_hist.py)
to do post-processing of the histograms produced by the Wang-Landau simulation
and, for comparison,
a constant-temperature MC program
[mc_chain_nvt_sw.py](mc_chain_nvt_sw.py).

### 9.2 Configuration-biased simulation of chain molecule
[mc_chain_nvt_cbmc_lj.py](mc_chain_nvt_cbmc_lj.py) and [mc_chain_lj_module.py](mc_chain_lj_module.py).

### 9.3 Gibbs ensemble simulation
[mc_gibbs_lj.py](mc_gibbs_lj.py).

### 11.1 Molecular dynamics using Lees–Edwards boundaries
[md_nvt_lj_le.py](md_nvt_lj_le.py) and [md_lj_le_module.py](md_lj_le_module.py).

### 11.2 Cell structure and linked lists in sheared boundaries
[md_lj_llle_module.py](md_lj_llle_module.py).
This may be used
in combination with [md_nvt_lj_le.py](md_nvt_lj_le.py),
edited so as to import routines from this module
instead of [md_lj_le_module.py](md_lj_le_module.py).
Despite the name, linked lists are _not_ used, as there is no advantage in Python.
However, the module illustrates the cell structure, and the calculation of forces between neighbouring cells,
analogous to the Fortran example.

### 12.1 Brownian dynamics program
[bd_nvt_lj.py](bd_nvt_lj.py) and [md_lj_module.py](md_lj_module.py).

### 12.2 Smart Monte Carlo simulation
[smc_nvt_lj.py](smc_nvt_lj.py) and [smc_lj_module.py](smc_lj_module.py).

### 12.3 Dissipative particle dynamics
[dpd.py](dpd.py) and [dpd_module.py](dpd_module.py).

### 13.1 Path-integral Monte Carlo
[qmc_pi_sho.py](qmc_pi_sho.py) for the simple harmonic oscillator.
[qmc_pi_lj.py](qmc_pi_lj.py) and [qmc_pi_lj_module.py](qmc_pi_lj_module.py) for Lennard-Jones.

### 13.2 Quantum Monte Carlo using a diffusion equation
[qmc_walk_sho.py](qmc_walk_sho.py).

### 14.1 Calculating the temperature and density profiles
A Fortran version of this code snippet appears in the text, not online.

### 14.2 Radial distribution function in a planar interface
[grint.py](grint.py).
We also provide a (comparatively large) file `grint_data.zip` containing test data
in the [Data repository](https://github.com/Allen-Tildesley/data).

### 14.3 Cluster analysis
[cluster.py](cluster.py).
We also provide a test configuration file [cluster.inp](cluster.inp).

### A.1 Utility modules
[averages_module.py](averages_module.py),
[config_io_module.py](config_io_module.py), and
[maths_module.py](maths_module.py),
used by many of the programs listed above.

### C.1 Programs to test forces and torques
[test_pot_atom.py](test_pot_atom.py) together with any one of
[test_pot_at.py](test_pot_at.py),
[test_pot_bend.py](test_pot_bend.py), or
[test_pot_twist.py](test_pot_twist.py).

[test_pot_linear.py](test_pot_linear.py) together with any one of
[test_pot_dd.py](test_pot_dd.py),
[test_pot_dq.py](test_pot_dq.py),
[test_pot_qq.py](test_pot_qq.py), or
[test_pot_gb.py](test_pot_gb.py).

### D.1 3D Fourier transform example
[fft3dwrap.py](fft3dwrap.py).

### F.1 Calculating the configurational temperature
A Fortran version of this code snippet appears in the text, not online.
However, a working Python example may be found in
[mc_nvt_lj.py](mc_nvt_lj.py) and its module [mc_lj_module.py](mc_lj_module.py).

### F.2 Correction to configurational temperature
A Fortran version of this code snippet appears in the text, not online.
However, a working Python example may be found in
[md_nve_lj.py](md_nve_lj.py) and its module [md_lj_module.py](md_lj_module.py).
