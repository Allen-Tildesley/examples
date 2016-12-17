# Brief Guide
Here are some notes to assist in running the programs.

## Data Input
Most of the Fortran codes use a `namelist` to input a few parameters from standard input.
This gives an easy way to specify default values in the program itself, and to use a
keyword-based syntax to specify values different from the default ones at run-time.
The input file, or string, should usually begin with `&nml` and end with `/`.
As a minimum, the program will expect to read `&nml /` (an empty list), but to
change the parameters, typical input might be
```
&nml nblock=20, nstep=10000, dt=0.002 /
```
and the `key=value` pairs may be set out on different lines if you wish.

The program `md_nve_lj_vl` expects the input file to contain a second namelist,
`&nml_list`, which should either be empty or contain the desired value of `r_list_factor`.
See the file `verlet_list_module.f90` for details.

## Initial Configuration
Simulation runs require a starting configuration which can usually be prepared using
the `initialize` program (built, by the default SConstruct file, in `build_initialize/`).
The default parameters produce an FCC configuration of 256 atoms at reduced density &rho;=0.75,
writing out just the positions (for an MC program) to a file `cnf.inp`.
If the parameter `velocities=.true.` is supplied, then positions and velocities are
written to the file, corresponding to a reduced temperature _T_ = 1.0.
These values of &rho; and _T_ (see below) lie in the liquid region of the Lennard-Jones phase diagram.
Non-default values may, of course, be supplied for this or other models.
The `cnf.inp` file may then be copied to the directory in which the run is carried out.
Typically, runs produce a final configuration `cnf.out`
(which may be renamed to `cnf.inp` as a starting point for further runs)
and intermediate configurations `cnf.001`, `cnf.002` etc during the run.

A utility program,
`adjust` takes in an MC or MD configuration and
scales the velocities to change the kinetic energy per atom by a specified amount,
and/or the positions (and the box length) to change the density by a specified amount.

## Lennard-Jones simulation programs
A large number of the examples simulate the Lennard-Jones liquid.
Before discussing these in detail,
we consider the different variants of the model that appear in the programs,
and the state points used for testing.

### State points for different Lennard-Jones models
For most of the examples, we use a cutoff of _R_<sub>c</sub>=2.5&sigma;.
For MC programs the cut (but not shifted) potential (c) is used in the simulation,
while for MD programs the cut-and-shifted potential (cs) is used.
In all cases, an attempt is made to produce results without, and with, long-range corrections,
the latter giving an estimate for the full potential (f).
Differences between the different models are discussed in various places,
see e.g.

* A Trokhymchuk, J Alejandre, _J Chem Phys,_ __111,__ 8510 (1999).

Using their table V as a guide, we take the critical point to be roughly located at:

LJ model                 | _T_<sub>crit</sub> | &rho;<sub>crit</sub> | _P_<sub>crit</sub>
-----                    | ----               | ----                 | ----
f: full, with LRC        | 1.31               | 0.31                 | 0.13
c: cut (but not shifted) | 1.19               | 0.32                 | 0.11
cs: cut-and-shifted      | 1.08               | 0.32                 | 0.09

At any temperature below _T_<sub>crit</sub>, the liquid state is bounded below by the
liquid-gas coexistence density, and using Tables II-IV of the same reference as a guide,
we take the values for three example temperatures as

LJ model                 |  _T_ | &rho;<sub>L</sub> | _P_
-----                    | ---- |  ----             | ----
f: full, with LRC        |  0.8 | 0.793             | 0.005
c: cut (but not shifted) |  0.8 | 0.765             | 0.008
cs: cut-and-shifted      |  0.8 | 0.730             | 0.013

LJ model                 |  _T_ | &rho;<sub>L</sub> | _P_
-----                    | ---- |  ----             | ----
f: full, with LRC        |  0.9 | 0.746             | 0.012
c: cut (but not shifted) |  0.9 | 0.714             | 0.020
cs: cut-and-shifted      |  0.9 | 0.665             | 0.030

LJ model                 |  _T_ | &rho;<sub>L</sub> | _P_
-----                    | ---- |  ----             | ----
f: full, with LRC        |  1.0 | 0.695             | 0.026
c: cut (but not shifted) |  1.0 | 0.652             | 0.036
cs: cut-and-shifted      |  1.0 | 0.578             | 0.062

For Lennard-Jones, the state point (&rho;,_T_)=(0.75,1.0)
lies in the liquid region of the phase diagram for all these variants of the model.

Approximate equations of state for the full LJ potential have been presented by

* JK Johnson, JA Zollweg, KE Gubbins, _Mol Phys_, __78,__ 591 (1993)
* J Kolafa and I Nezbeda, _Fluid Phase Equil,_ __100,__ 1 (1994)
* M Mecke, A Muller, J Winkelmann, J Vrabec, J Fischer, R Span, W Wagner,
_Int J Thermophys,_ __17,__ 391 (1996) and __19,__ 1493(E) (1998)

For testing our programs we used the more recent fitted equations of state presented by

* M Thol, G Rutkai, R Span, J Vrabec and R Lustig, _Int J Thermophys,_ __36,__ 25 (2015)
* M Thol, G Rutkai, A Koester, R Lustig, R Span, J Vrabec, _J Phys Chem Ref Data,_ __45,__ 023101 (2016)

for both the cut-and-shifted potential (denoted cs below), at _R_<sub>c</sub>=2.5&sigma;,
and the full potential (denoted f).
The formulae are implemented in the supplied program `eos_lj`.
For completeness, note that Thol et al also supply C++ programs, and tables of data,
in the Supplementary Information associated with their papers.
They are not responsible for our (Fortran) program!

### Lennard-Jones MD programs
First we look at the MD (and related) programs, which use the cut-and-shifted potential.
Here we compare with typical test runs from our programs using default parameters, _N_=256, except where stated.
Note that _E_ is the total internal energy per atom,
that _C_ is short for _C<sub>v</sub>_ (per atom) and _P_ is the pressure,
all including the ideal gas contributions.
The Smart Monte Carlo code `smc_nvt_lj` is included here since it uses the
cut-and-shifted potential which corresponds to the force calculation
(although it is not essential to do so).

Source                 | &rho;    | _T_       | _E_ (cs)   | _P_ (cs) | _C_ (cs)  | _E_ (f)    | _P_ (f)  | _C_ (f)  
------                 | -----    | -----     | --------   | -------- | --------- | -------    | -------  | --------
Thol et al (2015) (cs) | 0.75     | 1.00      | -2.9286    | 0.9897   |  2.2787   |            |          |          
Thol et al (2016) (f)  | 0.75     | 1.00      |            |          |           | -3.7212    | 0.3996   | 2.2630  
`bd_nvt_lj`            | 0.75     | 1.00      | -2.934(4)  | 0.974(7) |  2.26(8)  | -3.733(4)  | 0.373(7) | 2.27(8)
`md_nvt_lj`            | 0.75     | 1.00      | -2.940(4)  | 0.965(6) |  2.27(12) | -3.740(4)  | 0.363(6) | 2.27(12)
`md_npt_lj`            | 0.749(1) | 1.00      | -2.920(8)  | 0.99     |           | -3.719(9)  | 0.395(1) |
`md_nve_lj`            | 0.75     | 1.0022(3) | -2.9289    | 0.987(2) |  2.24(1)  | -3.7284    | 0.386(2) |          
`md_nve_lj_omp`        | 0.75     | 1.0027(2) | -2.9278    | 0.986(2) |  2.28(1)  | -3.7273    | 0.385(2) |          
`md_nve_lj_vl`         | 0.75     | 1.0023(3) | -2.9278    | 0.992(2) |  2.24(1)  | -3.7274    | 0.391(2) |          
`md_nve_lj_ll`&Dagger; | 0.75     | 1.0010(1) | -2.9272    | 0.992(1) |  2.28(1)  | -3.7268    | 0.391(1) |          
`smc_nvt_lj`           | 0.75     | 1.00      | -2.9300(5) | 0.971(2) |  2.263(5) | -3.7296(5) | 0.369(2) | 2.270(5)
`smc_nvt_lj` (100%)    | 0.75     | 1.00      | -2.928(2)  | 0.99(1)  |  2.26(2)  | -3.728(2)  | 0.39(1)  | 2.27(2)
`smc_nvt_lj` (30%)     | 0.75     | 1.00      | -2.930(3)  | 0.98(2)  |  2.26(3)  | -3.729(3)  | 0.38(2)  | 2.27(3)

* &Dagger; indicates a larger system size, _N_=864, needed to make the link-list method viable.
* The `md_nvt_lj` program seems to give a low energy and pressure, maybe needs looking at???
* The `md_npt_lj` program does not conserve well,
MSD of order 10<sup>-3</sup> rather than 10<sup>-8</sup> which is what we see for `md_nvt_lj`.
Are the barostat parameters poorly chosen??? Calculated _C<sub>p</sub>_ (cs)=5.0(4) while EOS gives 4.84.
* The `smc_nvt_lj` program was tested in default, single-particle-move, mode, with &delta;t=0.1,
in multi-particle mode, moving 100% of particles, with &delta;t=0.02,
and in multi-particle mode, moving 30% of particles, with &delta;t=0.03.
These values give acceptance rates in the 45% &ndash; 55% range.

Results for `md_lj_mts` are not directly comparable,
because a larger cutoff (by default _R<sub>c</sub>_=4.0&sigma;) is used to illustrate the method.
The program was tested with _N_=400 (box length 8.1).
The usual state point is simulated: &rho;=0.75 throughout.
No fitted EOS for the cs potential for this cutoff is available; obviously the estimates for the full potential
should be comparable with the values given above from Thol et al (2016).
Smallest timestep &delta;t (called `dt1` in the program)
and multiple-timestep-multipliers (the `n_mts` array) are given in column 2.
In all cases the run length was equal to 10 blocks of 20000 steps of length 0.005.
So, for example, 0.005 (142) means 10 blocks of 2500 steps of length 0.04,
each subdivided into 2 steps of length 0.02,
each subdivided again into 4 steps of length 0.005.
We also tried 0.002 (142) meaning 10 blocks of 6250 steps of length 0.016,
each subdivided into 2 steps of length 0.008,
each subdivided again into 4 steps of length 0.002.
The first line in the table below is from a run of `md_nve_lj` with the same system, for reference.

Source          | &delta;t    | _T_       | _E_ (cs)   | _P_ (cs) | _C_ (cs)  | _E_ (f)    | _P_ (f)  | CPU (mins) | _E_ (MSD)
-------         | --------    | -------   | ---------  | -------- | --------- | -------    | -------  | ------     | ------
`md_nve_lj`     | 0.005       | 1.0038(1) | -3.5199    | 0.557(2) | 2.26(1)   | -3.7157    | 0.410(2) | 30         | 1.7x10<sup>-8</sup>
`md_lj_mts` (0) | 0.005 (111) | 1.002(3)  | -3.5199    | 0.58(1)  | 2.4(1)    | -3.7157    | 0.43(2)  | 75         | 1.6x10<sup>-8</sup>
`md_lj_mts` (1) | 0.002 (142) | 1.0040(2) | -3.5196(2) | 0.559(1) | 2.26(1)   | -3.7153(2) | 0.412(1) | 50         | 1.1x10<sup>-7</sup>
`md_lj_mts` (2) | 0.005 (142) | 1.017(2)  | -3.491(4)  | 0.610(7) | 2.26(1)   | -3.686(4)  | 0.463(7) | 20         | 6.8x10<sup>-6</sup>
`md_lj_mts` (3) | 0.005 (142) | 1.0094(8) | -3.508(2)  | 0.576(3) | 2.26(1)   | -3.703(2)  | 0.429(3) | 30         | 7.8x10<sup>-7</sup>

* Run (0) was run with all the timesteps the same length, as a check of the program book-keeping.
The energy conservation is excellent, but the statistical errors on measured properties seem to be an order of magnitude
larger than for `md_nve_lj` and the pressure seems a bit off. Why is this so bad??
* Run (1), program default parameters, has quite good energy conservation, and statistical errors,
but this is partly a result of the smaller timestep.
* Run (2) uses the 0.005 timestep,
has worse energy MSD, and the average energies and pressure are deviating significantly from the expected values.
That's rather disappointing.
* Run (3) is identical to run (2) except that the switching length lambda is increased from 0.1 to 0.15.
This improves the energy conservation, but the pressure still looks a bit off.

Perhaps this program needs looking at???
Be aware that there have been some cosmetic changes to the introductory output since these first few test runs.

### Lennard-Jones MC programs
Our MC programs use the cut (but not shifted) potential
(which means that there is a delta correction, in the program, for the pressure).
In this case, the value of _C<sub>v</sub>_ (reported as _C_ below)
should be equal to the value for the full potential,
since the energy LRC is independent of temperature.
The Thol et al (2016) EOS for the full potential
is used to predict results for the cut (but not shifted) potential (denoted c),
at _R_<sub>c</sub>=2.5&sigma;, using the same LRC and delta corrections as in the MC codes.
Once again, all values in the table include the ideal gas contribution.
Except where indicated, tests are performed for _N_=256.

Source                 | &rho;     | _T_   | _E_ (c)    | _P_ (c)  | _E_ (f)    | _P_ (f)  | _C_ (f)
------                 | -----     | ----- | -------    | -------  | -------    | -------  | --------
Thol et al (2016) (f)  | 0.75      | 1.00  | -3.3197    | 0.7008   | -3.7212    | 0.3996   |  2.2630  
`mc_nvt_lj`            | 0.75      | 1.00  | -3.332(1)  | 0.651(3) | -3.734(1)  | 0.350(3) |  2.28(1)
`mc_nvt_lj_re`         | 0.75      | 1.00  | -3.332(1)  | 0.648(2) | -3.734(1)  | 0.347(2) |  2.258(4)
`mc_nvt_lj_ll`&Dagger; | 0.75      | 1.00  | -3.3230(3) | 0.669(1) | -3.7246(3) | 0.367(1) |  2.27(1)
`mc_npt_lj`            | 0.7501(2) | 1.00  | -3.331(1)  | 0.69     | -3.733(1)  | 0.364(2) |          
`mc_npt_lj_ll`&Dagger; | 0.7506(4) | 1.00  | -3.332(3)  | 0.69     | -3.734(3)  | 0.358(3) |          
`mc_zvt_lj`            | 0.7504(4) | 1.00  | -3.333(3)  | 0.668(4) | -3.735(3)  | 0.366(4) |          
`mc_zvt_lj_ll`&Dagger; | 0.7501(3) | 1.00  | -3.328(2)  | 0.669(2) | -3.729(2)  | 0.368(2) |          

* &Dagger; indicates a larger system size, _N_=864 (or approximately so for `mc_zvt_lj_ll`)
* The `mc_nvt_lj` program seems to give a low pressure, needs investigating? or longer run?
* The `mc_nvt_lj_re` program was run for four temperatures, see below for details.
* The `mc_npt_lj` _measured_ pressure (c) is 0.666(2) which is a little low. Needs checking?
Measured _C<sub>p</sub>_ (full) is 5.28(7) compared with Thol et al (2016) EOS giving 5.22
* The `mc_npt_lj_ll` program was run with `db_max`=0.015 to give a volume acceptance ratio around 9%.
Measured pressure (c) is 0.660(3) which is again a little low. Is the delta correction wrong somehow???
Or the constant-pressure algorithm???
Measured _C<sub>p</sub>_ (full) is 5.04(16) compared with Thol et al (2016) EOS value of 5.22.
The program probably needs making more resilient against changes in box size (list array reallocation).
* The `mc_zvt_lj` program was run at activity _z_=0.0795, the default value in the program, in a box of length 7&sigma;.
The Thol et al (2016) LRC-corrected value to give &rho;=0.75 would be _z_=0.080627.
Acceptance rate of creation/destruction moves is quite small, at about 0.3%.
For other state points see below.
We could look at including a reallocate_arrays routine to cope better with varying _N_.
* The `mc_zvt_lj_ll` program has the same acceptance ratio of moves.
It seems to run very slowly, which needs looking into???
Again, it would be more satisfying to use list array reallocation to make the program resilient to _N_ increasing.
* In principle, there should be a delta correction for the configurational temperature.
Long-range corrections are discussed by A Baranyai _J Chem Phys,_ __112,__ 3964 (2000) and by
A Lervik, O Wilhelmsen, TT Trinh, HR Nagel, _J Chem Phys,_ __143,__ 114106 (2015),
but they do not seem to discuss the truncation discontinuity.
This needs looking into ???

Tests for the grand canonical MC program were initially conducted at a slightly lower density,
very close to the liquid-vapour coexistence line (see Gibbs simulations below).
A box length of 7&sigma; was used, and creation/destruction acceptance ratios were around 1.5%.
Comparison was made with the Thol et al (2016) equation of state, with corrections for the cutoff.
The corresponding density is lower than the liquid coexistence density for the full potential,
so there is no guarantee that the EOS will be accurate
(it is only fitted in the single-phase regions of the full potential).

Source                |  z     | &rho;     | _T_  | _E_ (c)   | _P_ (c)    | _E_ (f)   | _P_ (f)  
-------               | ----   | -----     | ---- | --------- | -------    | -------   | -------  
Thol et al (2016) (c) | 0.032  | 0.65325   | 1.0  | -2.7212   | 0.0457     | -3.0710   | -0.1828  
`mc_zvt_lj`           | 0.032  | 0.6532(5) | 1.0  | -2.728(3) | 0.0325(25) | -3.078(4) | -0.196(2)

### Brownian dynamics program
The program `bd_nvt_lj` carries out a Brownian dynamics simulation for a set of atoms
interacting through the cut-and-shifted Lennard-Jones potential.
An initial configuration may be prepared, at a typical Lennard-Jones state point,
using the `initialize` program in the usual way.
As well as the usual run parameters, similar to a molecular dynamics code,
the user specifies a friction coefficient.
The calculated average thermodynamic quantities should be as expected for an
equilibrium simulation of this model at the chosen state point (see e.g. the table above).

### Gibbs Monte Carlo program
The program `mc_gibbs_lj` carries out Gibbs ensemble Monte Carlo,
and to test it we selected a temperature _T_=1.0,
which is below the critical point for the cut (but not shifted) LJ potential
(see tables above).
It was found convenient to start from a lower temperature,
with configurations at gas and liquid densities, with roughly equal numbers of particles,
and slowly work upwards in temperature, to equilibrate.
Note that the program expects two starting configurations: `cnf1.inp` and `cnf2.inp`.
Exchanges of box identity are expected as the critical temperature is approached,
and so one should not place blind trust in the separate box averages reported by the program,
but refer to histograms of density, energy etc.
At _T_=1.0, however, these exchanges of box identity are quite infrequent,
and the averages corresponded well to literature values for the coexistence parameters.
The production run corresponded to default parameters in the program.

Source  | &rho;<sub>L</sub> | &rho;<sub>G</sub> | _P_<sub>L</sub> | _P_<sub>G</sub> | _E_<sub>L</sub> (c) | _E_<sub>G</sub> (c)
-------              | -------- | -------- | -------  | -------- | --------------  | --------------
Trokhymchuk et al MC | 0.6542   | 0.0439   | 0.0336   | 0.0336   |                 |
Trokhymchuk et al MD | 0.6507   | 0.0500   | 0.0380   | 0.0380   | -2.713 &Dagger; | 1.047 &Dagger;
`mc_gibbs_lj`        | 0.652(1) | 0.050(1) | 0.028(1) | 0.038(1) | -2.730(5)       | 1.054(8)

* There is a small discrepancy between pressures in the two boxes. Is this expected?
* &Dagger; indicates values for given &rho; and _T_ from the Thol et al (2016) EOS (f) with cutoff correction.

### Replica exchange program
The `mc_nvt_lj_re` program conducts runs at several temperatures: four were used in testing.
The default program values include _T_=1.0, which is reported above, and here is the complete set,
with expected values from the Thol et al (2016) equation of state (f) corrected for cutoff.
As usual the program employed the cut (but not shifted) potential.
All runs are for density &rho;=0.75.
At the lowest temperature, the full-potential system would lie in the coexistence region.

Source                 | _T_    | _E_ (c)   | _P_ (c)  | _E_ (f)   | _P_ (f)   | _C<sub>v</sub>_ (f)
------                 | -----  | -------   | -------  | -------   | -------   | --------
Thol et al (2016) (f)  | 0.8772 | -3.6001   | 0.1942   | -4.0017   | -0.1070   |  2.3081  
`mc_nvt_lj_re`         | 0.8772 | -3.613(1) | 0.140(2) | -4.014(1) | -0.161(2) |  2.31(1)
Thol et al (2016) (f)  | 1.0000 | -3.3197   | 0.7008   | -3.7212   | 0.3996    |  2.2630  
`mc_nvt_lj_re`         | 1.0000 | -3.332(1) | 0.648(2) | -3.734(1) |  0.347(2) |  2.258(4)
Thol et al (2016) (f)  | 1.1400 | -3.0055   | 1.2571   | -3.4070   |  0.9559   |  2.2278  
`mc_nvt_lj_re`         | 1.1400 | -3.016(1) | 1.212(2) | -3.417(1) |  0.911(2) |  2.233(4)
Thol et al (2016) (f)  | 1.2996 | -2.6523   | 1.8667   | -3.0539   |  1.5655   |  2.1989  
`mc_nvt_lj_re`         | 1.2996 | -2.662(1) | 1.820(3) | -3.063(1) |  1.519(3) |  2.214(5)

## Lees-Edwards programs
The programs `md_nvt_lj_le` and `md_nvt_lj_llle` are intended to illustrate:
the moving boundaries used in nonequilibrium shear flow simulations;
an algorithm for integrating the SLLOD equations of motion with constrained kinetic energy;
and an adapted link-cell method required to handle the modified boundaries.
These programs both use the short-ranged WCA Lennard-Jones potential,
in order to compare results with the following papers:

* G Pan, JF Ely, C McCabe, DJ Isbister, _J Chem Phys,_ __122,__ 094114 (2005)
* KP Travis, DJ Searles, DJ Evans, _Mol Phys,_ __95,__ 195 (1998)

Testing was performed at the state point used in those papers: &rho;=0.8442, _T_=0.722.
A system size _N_=256 was used.
The given program defaults, including a time step of 0.005, were used throughout,
except for the strain rate which was varied.

Strain rate | _E_       | _P_       | &eta;
-----       | -----     | -----     | -----
0.04        | 1.8043(1) | 6.3907(7) | 2.30(4)
0.04        | 1.8039(3) | 6.389(2)  | 2.31(4)
0.16        | 1.8103(2) | 6.431(1)  | 2.254(8)
0.16        | 1.8098(2) | 6.427(1)  | 2.23(1)
0.64        | 1.8649(2) | 6.778(1)  | 1.939(2)
0.64        | 1.8646(2) | 6.776(1)  | 1.935(2)

In the table above, for each strain rate,
the first line comes from `md_nvt_lj_le`
and the second from `md_nvt_lj_llle`
(essentially identical, but roughly twice as fast for _N_=256).
In all cases the kinetic energy was conserved very accurately by the algorithm.
The results, particularly the increase in _E_ and _P_,
and the decrease in shear viscosity &eta;,
as the strain rate increases,
are in good agreement with the above papers.

* Although the link-list program seems to be working fine, it might be worth double-checking
the logic in the force routine, as it is a bit fiddly.

## Hard-particle programs
The programs `mc_nvt_hs` and `md_nve_hs` illustrate, respectively,
the simplest MC and MD methods for the basic hard-sphere model.
The temperature is not important in the first case: a factor _kT_ is used to normalize the energies.
The energy, in the second case, is identical with the (exactly conserved) kinetic energy,
and hence closely related to the temperature.
Equations of state for this model have been reported many times.
Here we refer to some fairly recent, useful, sources of data and/or fitted equations

* H Hansen-Goos, _J Chem Phys,_ __144,__ 164506 (2016)
* MN Bannerman, L Lue, LV Woodcock, _J Chem Phys,_ __132,__ 084507 (2010)
* J Kolafa, S Labik, A Malijevsky, _Phys Chem Chem Phys,_ __6,__ 2335 (2004)

The paper of Kolafa et al (2004) is particularly careful to discuss corrections
due to different ensembles and system size. Here we just present the raw results
for a small system, _N_=256; programs are run with default parameters.
Starting fcc lattice configurations may be prepared using `initialize` in
the usual way.
The EOS is taken from the Hansen-Goos (2016) paper, and a program to evaluate it
may be found in `eos_hs.f90`.

&rho; | _P_ (EOS) | _P_ `mc_nvt_hs`| _P_ `md_nve_hs` | &rho; `mc_npt_hs`
----- | -----     | -----    | ----- | -----
0.75  | 4.9910    | 4.960(7) | 4.985(4) | 0.749(2)
0.70  | 4.0087    | 3.996(8) | 4.005(3) | 0.700(2)
0.65  | 3.2171    | 3.210(8) | 3.215(2) | 0.651(3)
0.60  | 2.5769    | 2.573(4) | 2.573(1) | 0.600(3)
0.55  | 2.0574    | 2.051(3) | 2.055(1) | 0.553(3)
0.50  | 1.6347    | 1.634(2) | 1.632(1) | 0.502(2)

We must remember that _P_ is calculated by a box-scaling method in the _NVT_ simulation,
which may introduce a small systematic error. This can be reduced by reducing the
scaling factor, at the expense of worsening the statistics.
We also provide a program `mc_npt_hs` to illustrate the constant-_NPT_ method.
For the averages of &rho; reported above, the input pressure was that given by
the corresponding EOS entry.
With default parameters, volume move acceptance ratio was nearly 5% at the highest pressure,
and around 11% at the lowest pressure studied here.

We also provide two programs to simulate the hard spherocylinder model,
of cylinder length _L_ and diameter _D_:
`mc_npt_sc` and `mc_nvt_sc`.
In preparing configurations for these programs,
one must not allow overlap between the hard particles.
A simple approach is to
run `initialize` with `molecules="linear", random_positions=.t.`,
and to request a very low density.
For the default `length=5` spherocylinders
(_L_=5, _D_=1) a value of `density=0.05` is suitable.
Then,
the constant-pressure program may be used to compress the system to higher density.
This is a fairly slow process,
requiring the density &rho; and nematic order parameter _S_ to be carefully monitored.
Once suitable high-density state points have been prepared,
a configuration at a precisely specified density, for use in the constant-volume program,
may be obtained by a small expansion (using the `adjust` program).

For testing we use `N=256`;
such a small system is not recommended for serious work,
but allows us to explore up to &rho;=0.148
(box length 12 _D_, twice the interaction range of 6 _D_)
which is more than sufficient for our purposes.
Equations of state from MC simulations are presented in two papers

* SC McGrother, DC Williamson, G Jackson, _J Chem Phys,_ __104,__ 6755 (1996)
* PG Bolhuis, D Frenkel, _J Chem Phys,_ __106,__ 666 (1997)

In making comparisons, care must be taken with the units.
McGrother et al (1996) quote
densities in the form of the packing fraction &eta;=&rho; _v_<sub>mol</sub>
and pressures as _P_ _v_<sub>mol</sub>,
where _v_<sub>mol</sub> is the molecular volume.
We translate selected values from their Table V
(denoted (M) below)
into our reduced units based on _D_=1 below;
for _L/D_=5, _v_<sub>mol</sub>=4.4506.
(Bolhuis and Frenkel (1997) define reduced densities relative to
the density of closest packing of spherocylinders,
while reporting pressures the same way as McGrother et al.
We do not use the Bolhuis-Frenkel results below.)

_P_ _v_<sub>mol</sub> | &rho; _v_<sub>mol</sub> | _P_ | &rho; | _S_ | &rho; | _S_ | _P_ | _S_
----- | ----- | ----- | ----- | ----- | -----       | -----       | -----       | -----
(M)   | (M)   | (M)   | (M)   | (M)   | `mc_npt_sc` | `mc_npt_sc` | `mc_nvt_sc` | `mc_nvt_sc`
2.53  | 0.310 | 0.568 | 0.070 | 0.041 | 0.0698(2)   | 0.081(7)    | 0.579(2)    | 0.073(6)
3.63  | 0.352 | 0.816 | 0.079 | 0.053 | 0.0799(2)   | 0.098(6)    | 0.810(3)    | 0.098(6)
4.89  | 0.397 | 1.099 | 0.089 | 0.136 | 0.0895(2)   | 0.25(1)     | 1.126(5)    | 0.23(1)
5.05  | 0.400 | 1.135 | 0.090 | 0.170 | 0.0903(2)   | 0.32(3)     | 1.136(3)    | 0.149(7)&Dagger;
5.05  | 0.400 | 1.135 | 0.090 | 0.170 | 0.0923(3)   | 0.502(7)    | 1.061(3)    | 0.592(8)&Dagger;
5.40  | 0.419 | 1.213 | 0.094 | 0.574 | 0.0966(3)   | 0.764(6)    | 1.193(6)    | 0.596(3)
5.80  | 0.436 | 1.303 | 0.098 | 0.714 | 0.0984(2)   | 0.783(8)    | 1.325(6)    | 0.751(6)
6.20  | 0.448 | 1.393 | 0.101 | 0.754 | 0.1021(3)   | 0.839(3)    | 1.406(6)    | 0.795(4)

The `mc_npt_sc` runs use pressures from column 3 above;
the `mc_nvt_sc` runs are at densities taken from column 4.
At the highest pressure, using default parameters,
move acceptance ratio was around 30%,
and volume acceptance ratio around 10%.
These values rose to 50% and 15% respectively at the lowest pressure.
Several successive runs (10 blocks of 10000 steps each)
were undertaken for state points
near the isotropic-nematic transition,
where very slow evolution of the nematic order parameter could be observed.
Considerable sluggishness is expected around this transition,
giving irreproducible results, an example indicated by &Dagger; in the table.
Much longer runs are needed to achieve consistency!
Also the system size is about 25% that used by McGrother,
which has a direct effect on the measured nematic order parameter.
With these caveats in mind,
which apply mainly to the middle three state points in the table,
agreement between the two programs, and with the results of McGrother,
is reasonable.

* Perhaps we should increase the default number of steps in the programs? and/or increase the system size?

## Chain simulation programs
The program `mc_chain_nvt_cbmc_lj` simulates a single Lennard-Jones chain,
where the atoms are linked by harmonic springs.
There are, therefore,
both bonded and non-bonded interactions,
the former being used to select atom positions,
and the latter appearing in Rosenbluth weights,
which govern the acceptance/rejection of moves.
For comparison with the paper of Calvo, Doye and Wales, _J Chem Phys,_ __116,__ 2642 (2002),
test runs were carried out using _N_=13 atoms, a bond length of 1.122462&sigma;
(prepared using `initialize` with random non-overlapping atom positions)
and a rather low spring potential _k_<sub>spring</sub>=20.
We only use CBMC moves in this code: for a practical application it would be advisable
to include other kinds of move, for example crankshaft, pivot, and bridging moves.
Replica exchange (as used by Calvo et al) would also improve the sampling at low temperature.
Below we report the excess, potential, energy per atom _PE_,
and the excess heat capacity per atom _C<sub>v</sub>_(ex),
as well as the radius of gyration _R_<sub>g</sub>.
The program default run length is 10 blocks of 100000 steps.
For lower temperatures (below 0.40), longer runs (10 blocks of 1000000
steps) were used.
For temperatures higher than 1.0,
the bond length fluctuations become unphysically large for this value of _k_<sub>spring</sub>.

_T_   | _PE_      | _R_<sub>g</sub> | _C<sub>v</sub>_(ex)
----- | ------    | ------          | ------
0.26  | -2.044(7) | 1.074(1)        | 3.3(2)  
0.28  | -1.967(5) | 1.086(1)        | 4.16(9)
0.29  | -1.905(4) | 1.096(1)        | 4.6(1)  
0.30  | -1.893(7) | 1.098(1)        | 4.71(9)
0.31  | -1.819(3) | 1.1105(8)       | 4.34(6)
0.32  | -1.789(2) | 1.1162(5)       | 4.2(1)  
0.33  | -1.742(3) | 1.1254(5)       | 4.05(1)
0.34  | -1.705(4) | 1.133(1)        | 3.7(1)  
0.35  | -1.672(3) | 1.140(1)        | 3.49(8)
0.40  | -1.50(1)  | 1.173(3)        | 2.51(17)
0.45  | -1.40(1)  | 1.201(3)        | 2.26(8)  
0.50  | -1.297(8) | 1.224(1)        | 1.99(2)  
1.00  | -0.438(2) | 1.538(2)        | 1.371(3)

At the lowest temperatures, the acceptance rate of CBMC moves (with the default parameters) was around 2%,
while at _T_=0.35 it was around 11%, increasing further at higher temperatures.
The results are broadly in agreement with Calvo et al (2002) showing a similar sized peak in _C<sub>v</sub>_,
although at a somewhat lower temperature (0.30 as opposed to 0.35).

Here we give analogous results for the program default spring constant of _k_<sub>spring</sub>=400.

_T_   | _PE_      | _R_<sub>g</sub> | _C<sub>v</sub>_(ex)
----- | -------   | --------        | --------
0.26  | -2.076(3) | 1.069(1)        | 2.13(5)
0.28  | -2.027(5) | 1.075(1)        | 2.91(15)
0.29  | -1.993(4) | 1.081(1)        | 3.41(8)
0.30  | -1.958(4) | 1.087(1)        | 3.79(8)
0.31  | -1.915(4) | 1.094(1)        | 4.40(6)
0.32  | -1.859(6) | 1.105(1)        | 4.74(8)
0.33  | -1.816(4) | 1.114(1)        | 4.94(6)
0.34  | -1.765(5) | 1.125(1)        | 4.90(9)
0.35  | -1.719(4) | 1.135(1)        | 4.75(7)
0.40  | -1.505(5) | 1.183(1)        | 3.06(14)
0.45  | -1.379(3) | 1.212(1)        | 2.32(6)
0.50  | -1.266(3) | 1.240(1)        | 2.03(3)
1.00  | -0.459(1) | 1.512(1)        | 1.233(6)
2.00  |  0.387(2) | 1.850(1)        | 0.591(3)
5.00  |  1.986(3) | 2.035(2)        | 0.465(2)

Similar models were employed in `md_chain_nve_lj` and `md_chain_mts_lj`:
_N_=13 atoms and equilibrium bond length of 1.122462&sigma;.
Here we report results for constrained bond lengths, using the first program,
and for _k_<sub>spring</sub>=400 and 10000 (the program default value), using the second program.
In all cases, the primary indicator of a correctly-functioning program is energy conservation,
and this was checked in all cases.


constrained

_E_     | _T_       | _R_<sub>g</sub> | _C<sub>v</sub>_(ex)
-----   | -----     | -----           | -----
-2.0246 | 0.2485(2) | 1.06374(4)      | 2.176(6)
-1.9145 | 0.296(2)  | 1.073(1)        | 2.38(8)  
-1.6145 | 0.345(4)  | 1.125(2)        | 3.14(8)  
-1.3495 | 0.404(1)  | 1.182(2)        | 2.39(1)  
-1.2195 | 0.451(1)  | 1.207(1)        | 2.36(2)  
-1.0968 | 0.499(2)  | 1.234(1)        | 2.28(1)  
-0.1244 | 1.009(5)  | 1.471(5)        | 2.04(2)  
 1.0456 | 2.008(5)  | 1.754(9)        | 1.653(3)
 3.6459 | 4.996(4)  | 1.889(7)        | 1.534(1)

_k_<sub>spring</sub>=10000

_E_      | _T_       | _R_<sub>g</sub> | _C<sub>v</sub>_(ex)
-----    | -----     | -----           | -----
-1.7734  | 0.2496(2) | 1.0695(1)       | 2.96(4)
-1.3444  | 0.301(2)  | 1.144(2)        | 5.5(3)  
-1.1394  | 0.350(2)  | 1.173(2)        | 3.70(5)
-0.9494  | 0.399(1)  | 1.199(1)        | 3.19(5)
-0.7694  | 0.448(1)  | 1.230(2)        | 3.09(3)
-0.5943  | 0.497(2)  | 1.262(3)        | 3.04(5)
 0.7857  | 1.000(4)  | 1.467(12)       | 2.50(3)
 2.8858  | 1.98(2)   | 1.752(14)       | 2.08(6)
 8.3859  | 5.04(4)   | 1.904(2)        | 1.94(2)

_k_<sub>spring</sub>=400

_E_     | _T_       | _R_<sub>g</sub> | _C<sub>v</sub>_(ex)
-----   | -----     | -----           | -----
-1.7934 | 0.2487(2) | 1.0646(1)       | 2.89(5)
-1.6250 | 0.299(1)  | 1.0753(4)       | 3.24(7)
-1.2900 | 0.346(5)  | 1.127(3)        | 5.7(3)
-0.9942 | 0.401(2)  | 1.179(2)        | 3.58(8)
-0.8042 | 0.451(2)  | 1.208(2)        | 3.21(4)
-0.6558 | 0.497(2)  | 1.224(2)        | 3.03(3)
 0.7565 | 0.995(3)  | 1.447(6)        | 2.55(3)
 2.9036 | 2.006(6)  | 1.757(9)        | 2.08(1)
 8.3488 | 5.00(1)   | 1.92(1)         | 1.97(1)

When comparing results with the MC program, several points should be remembered.

1. Constraining the bond lengths affects average potential energy, kinetic energy, and heat capacity.
2. While we use _k_<sub>spring</sub>=10000 to highlight the multiple timestep method,
it is quite likely that energy flow between bond vibrations and other degrees of freedom will be inefficient,
due to the timescale separation.
3. The constant-_NVE_ and constant-_NVT_ ensembles are expected to yield different behaviour around the collapse transition.
4. Molecular dynamics is not expected to thoroughly explore the energy landscape at low temperatures,
giving instead (typically) quasi-harmonic vibrations in a single basin.
The default run lengths are fairly modest here: 10 blocks,
each consisting of 100000 steps of length &delta;t=0.002.

For the hard-sphere square-well chain, the aim was to show the operation of the Wang-Landau method.
Here we used pivot and crankshaft moves as well as CBMC regrowth.
In a practical application it would be advisable to include some bridging moves as well.
Reasonably long chains have been studied by Taylor, Paul and Binder, _J Chem Phys,_ __131,__ 114907 (2009),
who provide references to earlier simulation work, as well as exact results for very short chains.
Here we choose _N_=13, bond length equal to &sigma;, and a nonbonded interaction range of 1.5&sigma;.
The starting chain configuration can be prepared using `initialize` in the usual way.

As a reference for comparison, we ran a set of canonical ensemble calculations with `mc_chain_nvt_sw`.
The program default is to run for 10 blocks, each of 100000 steps;
this was increased to 10 blocks of 500000 steps for temperatures below 0.25.
The results are shown on the left of the following table.

_T_   |  _PE_ | _R_<sub>g</sub> | _C<sub>v</sub>_(ex) | _PE_ | _R_<sub>g</sub> | _C<sub>v</sub>_(ex)
----- | ------    | ------    | ------     | ------  | ------ | ------
method | _NVT_    | _NVT_     | _NVT_      |   WL    |   WL   |  WL
0.15  | -2.81(1)  |  1.070(2) |  1.1(2)    | -2.814  |  1.068 |  2.053
0.18  | -2.759(8) |  1.072(2) |  2.2(2)    | -2.744  |  1.073 |  2.498
0.20  | -2.699(8) |  1.077(2) |  2.4(1)    | -2.694  |  1.078 |  2.366
0.22  | -2.645(4) |  1.082(1) |  2.16(7)   | -2.649  |  1.082 |  2.226
0.25  | -2.586(8) |  1.090(2) |  2.15(9)   | -2.584  |  1.089 |  2.121
0.30  | -2.482(6) |  1.104(2) |  1.97(6)   | -2.481  |  1.103 |  2.010
0.50  | -2.127(2) |  1.161(1) |  1.50(2)   | -2.128  |  1.161 |  1.491
1.00  | -1.547(2) |  1.318(1) |  1.024(6)  | -1.543  |  1.319 |  1.020
2.00  | -0.885(1) |  1.658(1) |  0.330(2)  | -0.883  |  1.660 |  0.332
5.00  | -0.566(1) |  1.889(1) |  0.0318(1) | -0.565  |  1.890 |  0.032

Results obtained from a run of the Wang-Landau program `mc_chain_wl_sw`,
using the same model, are given on the right of the table above.
The program was run with default parameters:
the flatness criterion was set at 80% and there were 20 stages during which the entropy
modification constant `ds` was halved at each stage.
The results are from the histograms produced in the 20th stage.
This analysis can also be performed (for any desired temperature) by the program `wl_hist`, after the run.
The results are generally in good agreement with the canonical ensemble test runs.
The most significant discrepancies are in the heat capacities at the lowest two temperatures,
which reflects the poor sampling of the canonical ensemble program (with these basic MC moves),
and (probably to a lesser extent) the sampling problems about to be discussed.

This particular test run illustrated one drawback of the simplest Wang-Landau implementation:
two low-lying energies (corresponding to _q_=38 and 39 square-well interactions) were discovered
during the very last stage, in which `ds` is very small.
Accordingly, the system remained stuck in these low-energy states for a very long time,
until their tabulated entropy _s(q)_ reached a high enough value to allow _q_ to change;
even then, the final weight of the lowest state in the final "flat" histogram
could not be considered completely reliable.
The overall run length was of order 90000 blocks of 10000 steps each,
most of it spent in stage 20.

If the run were repeated with the same parameters,
one cannot guarantee that the same result will be obtained:
indeed, it is more likely that the run would conclude much earlier, within a few thousand blocks,
without ever discovering these low-lying energies.
This would affect the results, particularly at the lower temperatures,
and of course there would be no indication of anything wrong.
It is always a danger with any Monte Carlo method, including Wang-Landau,
that inaccessible configurations will not be sampled.
Various improvements of the method may be found in the literature.
(The reader is invited to work out the structure of the 13-atom chain configuration
with hard sphere diameter and bond length both equal to &sigma;
having 39 nonbonded interactions within range 1.5&sigma;.
_Hint:_ it is not one of the highest-symmetry clusters such as the icosahedron;
however, it is based on a fairly symmetric cluster,
with small distortions due to the fixed bond length and the need to maximise interactions.)

* We might consider adding a fixed-length production run at the end of `mc_chain_wl_sw`.

## Polyatomic Lennard-Jones program
The program `mc_nvt_poly_lj` conducts Monte Carlo simulations of a system of rigid molecules
composed of Lennard-Jones interaction sites.
For simplicity the sites are taken to be identical, although the program is easily generalized.
Molecular orientations are represented by quaternions,
which are used to calculate the rotation matrix
and hence the LJ site positions.

We test this with the three-site model of orthoterphenyl, a fragile glassformer,
described in the following publications amongst others.

* LJ Lewis, G Wahnstrom, _Sol State Commun,_ __86,__ 295 (1993)
* LJ Lewis, G Wahnstrom, _Phys Rev E,_ __50,__ 3865 (1994)
* S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, _Phys Rev E,_ __65,__ 041205 (2002)
* E La Nave, S Mossa, F Sciortino, P Tartaglia, _J Chem Phys,_ __120,__ 6128 (2004)

The sites are arranged at the vertices of an isosceles triangle with bond angle 75 degrees,
LJ parameters &epsilon;/_k_<sub>B</sub> = 600K
or &epsilon; = 5 kJ mol<sup>-1</sup> to a good approximation,
&sigma;=0.483nm,
and bond lengths equal to &sigma;.

Tests were performed at &rho;=0.32655 which is equivalent to &rho<sub>4</sub>=1.108g cm<sup>-3</sup>
in Mossa et al (2002).
Comparisons of potential energy (_PE_=_E_-3_T_ converted to kJ/mol with a factor 5)
were made with the fit given by eqn (6) of that paper.

_T_   | _E_       | _P_      | _T_ (K) | _PE_ (kJ/mol) | eqn (6)
----- | -----     | -----    | -----   | -----         | -----
0.5   | -14.30(1) | 1.63(3)  |  300    | -79.00(5)     | -77.52
1.0   | -11.21(1) | 5.52(3)  |  600    | -71.05(5)     | -68.55
1.5   | -8.262(7) | 8.95(2)  |  900    | -63.81(4)     | -61.19
2.0   | -5.52(1)  | 11.86(3) | 1200    | -57.60(5)     | -54.69

Exact agreement is not expected because the potential of Mossa et al (2002) has a different
cutoff correction, but the agreement is reasonable.
A second set of tests were performed at _T_=0.63333=380K
at the specified densities &rho;<sub>1</sub>, &hellip; &rho;<sub>5</sub>.
Here the excess pressure (_P_(ex)=_P_-&rho;_T_ converted to MPa
with a factor 73.54 based on the values of &epsilon; and &sigma;)
is compared with the fit given by eqn (28) and the coefficients in Table III of Mossa et al (2002).
NB the volumes to insert into the equation are those of their Table I,
which are specific to their system size.

Id    | &rho;   | _E_        | _P_     | _P_(ex) (MPa) | eqn (28)
----- | -----   | -----      | -----   | -----         | -----
1     | 0.30533 | -12.737(7) | 0.35(2) | 12(2)         | 19.077
2     | 0.31240 | -13.025(7) | 0.99(2) | 58(2)         | 60.143
3     | 0.31918 | -13.293(8) | 1.66(2) | 107(2)        | 112.798
4     | 0.32655 | -13.50(1)  | 2.60(2) | 176(1)        | 177.222
5     | 0.33451 | -13.65(1)  | 3.95(2) | 274(1)        | 253.510

Although not perfect at the ends of the range, the agreement is not bad;
once more, the difference in cutoff correction should be borne in mind.

* The comparisons are perhaps not the best, since Mossa et al (2002) used a different value
&epsilon;=5.276 kJ/mol which, with their cutoff term, gives a well depth 4.985 kJ/mol.
Maybe better to compare with another paper,
or modify the potential to compare more precisely with this one??
* Interestingly, in the book, it seems we never offer an example of a polyatomic molecular dynamics program,
for either linear or nonlinear molecules. It is tempting to write examples of both,
but I cannot give this a high priority.

## DPD program
For the `dpd` example, we recommend generating an initial configuration
using the `initialize` program, with namelist input similar to the following
```
&nml n = 100, density = 3.0, random_positions = .true.,
velocities = .true., soft=.true. /
```
The above value of the density is typical when using this method to model water.

For testing we compare with an approximate DPD equation of state for _P_.

* RD Groot, PB Warren, _J Chem Phys,_ __107,__ 4423 (1997)
* TP Liyana-Arachchi, SN Jamadagni, D Eike, PH Koenig, JI Siepmann,
_J Chem Phys,_ __142,__ 044902 (2015)

The paper of Liyana-Arachchi et al (2015) is an improvement of the original
EOS of Groot and Warren (1997), which is more accurate and
applicable over a wider range of state points.
The function is included in the `dpd` program,
and the expected value of _P_ (labelled EOS below)
is printed for comparison at the end.
We give results obtained by both
the Lowe thermostat (L) and the Shardlow algorithm (S).
We take the default values of _a_ &rho;/T=75, and of other parameters not mentioned below.

 _T_   | &rho; | _P_ (EOS) | _P_ (L)   | _P_ (S)
 ----- | ----- | -----     | -----     | -----
 0.5   | 3.0   | 11.864    | 11.814(2) | 11.819(2)
 1.0   | 3.0   | 23.587    | 23.637(2) | 23.635(2)
 1.5   | 3.0   | 35.276    | 35.449(3) | 35.455(4)
 2.0   | 3.0   | 46.951    | 47.257(4) | 47.265(5)
 1.0   | 2.0   | 14.187    | 14.320(2) | 14.316(2)
 1.0   | 4.0   | 32.811    | 32.622(3) | 32.628(3)
 1.0   | 5.0   | 41.887    | 41.539(4) | 41.533(3)

## Test programs for potentials, forces and torques
Two program files are provided: `test_pot_atom.f90` and `test_pot_linear.f90`,
for pair potentials between, respectively, atoms and linear molecules.
These are combined, as appropriate, with modules which contain a subroutine to calculate
the necessary potential, forces and torques.
The aim is to demonstrate the numerical testing of the analytical derivatives
which go into the forces and torques:
small displacements and rotations are applied in order to do this.
The test is performed for a randomly selected configuration.
Some parameters are used to prevent serious overlap,
which might produce numerical overflow,
while keeping the particles close enough together to give non-zero results.
The values of these parameters may be adjusted via the namelist in individual cases;
to run the programs without any tweaking,
simply give an empty namelist `&nml /` to standard input in the usual way.
The supplied examples are:

* `test_pot_at` the Axilrod-Teller three-body potential
* `test_pot_bend` the angle-bending part of a polymer chain potential
* `test_pot_dd` the dipole-dipole potential
* `test_pot_dq` the dipole-quadrupole and quadrupole-dipole potential
* `test_pot_gb` the Gay-Berne potential
* `test_pot_qq` the quadrupole-quadrupole potential
* `test_pot_twist` the angle-torsion part of a polymer chain potential

In all cases, the SConstruct file builds these programs in a directory
whose name is taken from the module name above,
but the executable file is named `test_pot_atom` or `test_pot_linear`
as appropriate.

## T-tensor program
The program `t_tensor` compares the calculation of multipole energies by two methods:
using explicit formulae based on trigonometric functions of the Euler angles,
and via the Cartesian T-tensors.
Two linear molecules are placed in random positions and orientations,
within a specified range of separations,
and some of the contributions to the electrostatic energies and forces are calculated.
The program may be run using an empty namelist `&nml /`,
so as to take the program defaults,
or various parameters may be specified in the namelist.

* How easy would it be to add quadrupole-quadrupole energy, quadrupole-dipole forces,
quadrupole-quadrupole forces, and all the torques, calculated both ways??

## Cluster program
The `cluster` program is self contained. It reads in a configuration of atomic positions
and produces a circular linked list for each cluster identified within it.
The best value of the critical separation `r_cl` depends on the particular physical system
being considered. The supplied default will, for a liquid-state Lennard-Jones configuration,
most likely identify almost all atoms as belonging to a single cluster, with at most
a few atoms being isolated on their own. However, this is unlikely to be the type of system
for which this analysis is useful.

* Should we provide a more illuminating test case?

## Correlation function program
The aim of the program `corfun` is to illustrate the direct method, and the FFT method,
for calculating time correlation functions.
The program is self contained: it generates the time dependent data itself,
using a generalized Langevin equation,
for which the time correlation function is known.
The default parameters produce a damped, oscillatory, correlation function,
but these can be adjusted to give monotonic decay,
or to make the oscillations more prominent.
If the `origin_interval` parameter is left at its default value of 1,
then the direct and FFT methods should agree with each other to within numerical precision.
The efficiency of the direct method may be improved,
by selecting origins less frequently,
and in this case the results obtained by the two methods may differ a little.

Sample results using default program parameters are shown here.
The direct method is indicated in black, plotting only every fifth point for clarity.
The FFT result is shown as a red line: it matches the direct method as expected.
The exactly known function is a blue line.
There are very small discrepancies with the results of the simulation,
due to the finite length of the latter.

![alt text](corfun.png "corfun test results")

## Diffusion program
The program `diffusion` reads in a sequence of configurations and calculates
the velocity auto correlation function (vacf),
the mean square displacement (msd), and
the cross-correlation between velocity and displacement (rvcf).
Any of these may be used to estimate the diffusion coefficient,
as described in the text.
The output appears in `diffusion.out`
It is instructive to plot all three curves vs time.

The input trajectory is handled in a crude way,
by reading in successive snapshots with filenames `cnf.000`, `cnf.001`, etc.
These might be produced by a molecular dynamics program,
at the end of each block,
choosing to make the blocks fairly small (perhaps 10 steps).
As written, the program will only handle up to `cnf.999`.
Obviously, in a practical application,
a proper trajectory file would be used instead of these separate files.

It is up to the user to provide the time interval between successive configurations.
This will typically be a small multiple of the timestep used in the original simulation.
This value `delta` is only used to calculate the time, in the first column of
the output file.
A default value of 0.05 is provided as a place-holder, but
the user really should specify a physically meaningful value;
forgetting to do so could cause confusion when one attempts
to quantify the results.

To make it easier to test this program,
we have also supplied a self-contained program `diffusion_test`,
which generates an appropriate trajectory by numerically solving
the simple Langevin equation for _N_ non-interacting atoms (_N_=250 by default).
For this model, one specifies the temperature and friction coefficient,
which dictates the rate of exponential decay of the vacf,
and hence the diffusion coefficient.
The exact results for the vacf, rvcf and msd are written out to `diffusion_exact.out`
for easy comparison with `diffusion.out`.
Here are some typical results using default program parameters throughout.
The vacf is in red, rvcf in blue, and msd in green;
every fifth point is shown for the results of `diffusion`,
while the exact results are indicated as lines.
For the default program parameters, the diffusion coefficient is _D_=1.

![alt text](diffusion.png "diffusion test results")

## Interface pair correlation function
The program `grint.f90` reads in a set of configurations and calculates
the pair correlation function for a system that is inhomogeneous in the
z direction. It is assumed that the configurations consist of a liquid slab,
surrounded by gas, with the interfaces lying in the _xy_-plane.
The two interfaces are located by fitting the instantaneous density profile
to a difference of two tanh functions. Then the single-particle density function,
relative to the interface locations, is calculated.
To make the process as robust as possible, an initial guess at the mid-point
of the liquid slab should be provided, and this is updated automatically as
successive configurations are read in, so as to shift the liquid slab into
the middle of the periodic box, before fitting.
Also, the results of one fit are passed on as the starting point of the next one.
The program handles cubic boxes only;
the modifications necessary to handle non-cubic boxes are fairly easy to make,
but we will not do that here.

* This program is incomplete. The essential parts, of reading in and
fitting the profiles, and calculating the single-particle density in coordinates
relative to the interface position, are done, and should not require any changes.
A first attempt has been made to calculate the two-particle density: however
this calculation, including the scaling and normalization factors, need to be
checked. It was also promised, in the book, to convert this to a two-particle
correlation function; however, this requires calculating the single-body density
at a wider range of z than that of the two-body function (we need rho at z+r*c for
all tabulated z, r and c). Over to you Dominic!
* The program, as it stands, has been tested on a set of 101 configurations
of 3200 atoms in a 20x20x20 box. It was produced by `mc_nvt_lj_ll`, at a temperature
_T_=0.70: 100 blocks of 5000 steps each. These configurations are available if
needed for further testing.

## Error calculation
The program `error_calc` is a self-contained illustration of the effects of
correlations on the estimation of errors for a time series.
We produce the series using a generalized Langevin equation,
in the same manner as for the correlation function program (see above).
Since the correlation time of the GLE is exactly known,
we can predict the effects, and compare with the empirical estimates
obtained by different methods.
The program contains extensive comments to explain what is being calculated at each stage.

## FFT program
The aim of `fft3dwrap` is to illustrate the way a standard Fast Fourier Transform
library routine is wrapped in a user program.
We numerically transform a 3D Gaussian function,
and compare with the analytically, exactly, known result,
User input defines the number of grid points and the box size;
sensible defaults are provided.
The library that we use for this example is [FFTW](http://www.fftw.org/).

## Hit-and-miss and sample-mean
The two programs `hit_and_miss` and `sample_mean` illustrate two very simple
Monte Carlo methods to estimate the volume of a 3D object.
They are both described in detail at the start of Chapter 4.
No user input is required.
For the built-in values defining the geometry, the exact result is 5/3.
