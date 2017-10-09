# Brief Guide
Here are some notes to assist in running the programs.
Some details of tests carried out on the programs are also given.
Do not expect to duplicate these results,
they are simply a guide as to the kind of behaviour to expect.
If you find a definite programming error,
please report it via the [Issues](https://github.com/Allen-Tildesley/examples/issues) tab above
(you will need to be signed into GitHub to do this).

## Data Input
Most of the Fortran codes use a `namelist` to input a few parameters from standard input.
This gives an easy way to specify default values in the program itself, and to use a
keyword-based syntax to specify values different from the default ones at run-time.
The input file, or string, should usually begin with `&nml` and end with `/`.
As a minimum, the program will expect to read `&nml /` (an empty list),
and on a Linux/Unix/bash system a typical command line would be
```
echo '&nml /' | ./mc_nvt_lj
```
or similar.
To change the parameters, the namelist input might include the new values like this
```
&nml nblock=20, nstep=10000, temperature=2.0 /
```
Alternatively the namelist may be supplied in a file,
for instance `run.inp` containing
```
&nml
nblock=20,
nstep=10000,
temperature=2.0
/
```
and the program run like this
```
./mc_nvt_lj < run.inp
```
As indicated, the `key=value` pairs may be set out on different lines if you wish.

## Initial Configuration
Simulation runs for bulk liquids require a starting configuration which can usually be prepared using
the `initialize` program (built, by the default SConstruct file, in `build_initialize/`).
The default parameters produce an FCC configuration of 256 atoms at reduced density &rho;=0.75,
writing out just the positions (for an MC program) to a file `cnf.inp`.
You would run the program as described above,
for instance on a Linux/Unix/bash system like this
```
echo '&nml /' | ./initialize
```
If the parameter `velocities=.true.` is supplied within the namelist,
then positions and velocities are written to the file,
corresponding to a reduced temperature _T_ = 1.0.
These values of &rho; and _T_ (see below) lie in the liquid region of the Lennard-Jones phase diagram.
Non-default values may, of course, be supplied for this or other models.
The `cnf.inp` file may then be copied to the directory in which the run is carried out.
Typically, runs produce a final configuration `cnf.out`
(which may be renamed to `cnf.inp` as a starting point for further runs)
and intermediate configurations `cnf.001`, `cnf.002` etc during the run.

Some of the programs simulate a single chain of atoms, without periodic boundary conditions.
Initial configurations for these may also be prepared using the `initialize` program
selecting `molecules="chain"`, an appropriate number of atoms, for example `n=13`,
and `velocities=.true.` if required. There is an option `constraints=.true.` if the velocities
should be chosen with constraints applied relative to the bonds between neighbouring atoms in the chain.

A utility program, `adjust`, takes in an MC or MD configuration and
scales the velocities to change the kinetic energy per atom by a specified amount,
and/or the positions (and the box length) to change the density by a specified amount.
You may prefer to write your own program or script to perform these types of operation.

## Visualizing configurations
Our simulation configuration files have a very simple format:
the first line is the number of molecules,
the second line is the periodic box length (we always use cubic boxes),
or occasionally the bond length (for chain molecules which are simulated without a box),
and the third and subsequent lines each contain the coordinates of one atom or molecule.
The first three numbers on each of these lines are always the (x,y,z) position,
in simulation units (e.g. Lennard-Jones &sigma;=1).
The rest of the numbers on each line contain velocities, orientation vectors etc.,
as appropriate.

This format is not compatible with most molecular visualization programs,
such as [JMOL](http://jmol.sourceforge.net/)
or [VMD](http://www.ks.uiuc.edu/Research/vmd/).
However, conversion into the basic [XYZ](https://en.wikipedia.org/wiki/XYZ_file_format) format
is easily achieved using a simple program, script, or editor. For example,
on most Linux/Unix systems, the `awk` language will be available:
```
awk '(NR==1) {print} (NR==2) {print "Comment line"} (NR>2) {printf "%5s%15.6f%15.6f%15.6f\n", ("Ar"),($1*3.4),($2*3.4),($3*3.4)}' cnf.inp > cnf.xyz
```
This produces a file which should be recognized as a set of Argon atoms with positions
(in Angstroms) appropriate to their van der Waals diameter.
This can be read into a molecular visualizer,
which will typically have various options for representing the atoms.
No attempt is made here to represent the periodicity of the system in the `cnf.xyz` file.

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

### Lennard-Jones MD, BD and SMC programs
First we look at the MD (and related) programs, which use the cut-and-shifted potential.

The first test of the MD codes is that energy, or the appropriate energy-like variable, is conserved.
The following table uses runs of 10 blocks,
each consisting of a number of steps to give 16 units of time per block
with the indicated timestep (e.g. 1000&times;0.016, 2000&times;0.008 etc).
We report the MSD values of the conserved variable for each program.

&delta;t | `md_nve_lj` | `md_nvt_lj` | `md_npt_lj`
-------- | --------    | --------    | --------
0.016    | 4.3284&times;10<sup>-6</sup>  | 4.0409&times;10<sup>-6</sup>  | 4.5875&times;10<sup>-6</sup>
0.008    | 1.8430&times;10<sup>-7</sup>  | 2.0036&times;10<sup>-7</sup>  | 2.4402&times;10<sup>-7</sup>
0.004    | 1.3121&times;10<sup>-8</sup>  | 1.4494&times;10<sup>-8</sup>  | 1.4956&times;10<sup>-8</sup>
0.002    | 1.2621&times;10<sup>-9</sup>  | 1.3526&times;10<sup>-9</sup>  | 1.5914&times;10<sup>-9</sup>
0.001    | 1.8530&times;10<sup>-10</sup>  | 2.1371&times;10<sup>-10</sup>  | 2.1005&times;10<sup>-10</sup>

Log-log plots show the expected dependence MSD &prop; &delta;t<sup>4</sup>,
hence RMSD &prop; &delta;t<sup>2</sup>,
except for some small deviations at the smallest timestep.

Now we compare EOS data
with typical test runs from our programs using default parameters, _N_=256, except where stated.
Note that _E_ is the total internal energy per atom,
that _C_ is short for _C<sub>v</sub>_ (per atom) and _P_ is the pressure,
all including the ideal gas contributions.
The Smart Monte Carlo code `smc_nvt_lj` is included here since it uses the
cut-and-shifted potential which corresponds to the force calculation
(although it is not essential to do so).
Similarly, we include here the Brownian dynamics program `bd_nvt_lj`.

Numbers in parentheses (here and in the following tables)
indicate errors in the last quoted digit, estimated from block averages.
Results without error estimates are fixed (such as the temperature or density) or conserved.

Source                 | &rho;     | _T_       | _E_ (cs)   | _P_ (cs) | _C_ (cs)  | _E_ (f)    | _P_ (f)  | _C_ (f)  
------                 | -----     | -----     | --------   | -------- | --------- | -------    | -------  | --------
Thol et al (2015) (cs) | 0.75      | 1.00      | -2.9286    | 0.9897   |  2.2787   |            |          |          
Thol et al (2016) (f)  | 0.75      | 1.00      |            |          |           | -3.7212    | 0.3996   | 2.2630  
`md_nvt_lj`            | 0.75      | 1.00      | -2.940(4)  | 0.965(6) |  2.27(12) | -3.740(4)  | 0.363(6) | 2.27(12)
`md_npt_lj`&sect;      | 0.7514(6) | 1.00      | -2.947(6)  | 0.995(1) |           | -3.748(7)  | 0.391(1) |
`md_nve_lj`            | 0.75      | 1.0022(3) | -2.9289    | 0.987(2) |  2.24(1)  | -3.7284    | 0.386(2) |          
`md_nve_lj_omp`        | 0.75      | 1.0027(2) | -2.9278    | 0.986(2) |  2.28(1)  | -3.7273    | 0.385(2) |          
`md_nve_lj_vl`         | 0.75      | 1.0023(3) | -2.9278    | 0.992(2) |  2.24(1)  | -3.7274    | 0.391(2) |          
`md_nve_lj_ll`&Dagger; | 0.75      | 1.0010(1) | -2.9272    | 0.992(1) |  2.28(1)  | -3.7268    | 0.391(1) |          
`md_nvt_lj_ll`&Dagger; | 0.75      | 1.00      | -2.927(2)  | 0.994(3) |  2.3(1)   | -3.727(2)  | 0.392(3) | 2.3(1)         
`smc_nvt_lj`&sharp;(a) | 0.75      | 1.00      | -2.9300(5) | 0.971(2) |  2.263(5) | -3.7296(5) | 0.369(2) | 2.270(5)
`smc_nvt_lj`&sharp;(b) | 0.75      | 1.00      | -2.928(2)  | 0.99(1)  |  2.26(2)  | -3.728(2)  | 0.39(1)  | 2.27(2)
`smc_nvt_lj`&sharp;(c) | 0.75      | 1.00      | -2.930(3)  | 0.98(2)  |  2.26(3)  | -3.729(3)  | 0.38(2)  | 2.27(3)
`bd_nvt_lj`            | 0.75      | 1.00      | -2.934(4)  | 0.974(7) |  2.26(8)  | -3.733(4)  | 0.373(7) | 2.27(8)

&Dagger; Indicates a larger system size, _N_=864, needed to make the link-list method viable. Note that
the speedup is not enormous for this system size, corresponding to 4&times;4&times;4 cells.

&sect; The constant-pressure simulation was run at _P_=0.99, the program default.

&sharp; The `smc_nvt_lj` program was tested (a) in default, single-particle-move, mode, with &delta;t=0.1;
(b) in multi-particle mode, moving 100% of particles, with &delta;t=0.02;
and (c) in multi-particle mode, moving 30% of particles, with &delta;t=0.03.
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
This example is just to illustrate the idea:
most of the test runs are actually slower, not faster, than `md_nve_lj`.

Source          | &delta;t    | _T_       | _E_ (cs)   | _P_ (cs) | _C_ (cs)  | _E_ (f)    | _P_ (f)  |  _E_ (MSD)
-------         | --------    | -------   | ---------  | -------- | --------- | -------    | -------  |  ------
`md_nve_lj`     | 0.005       | 1.0038(1) | -3.5199    | 0.557(2) | 2.26(1)   | -3.7157    | 0.410(2) | 1.7&times;10<sup>-8</sup>
`md_lj_mts`&dagger; | 0.005 (111) | 1.002(3)  | -3.5199    | 0.58(1)  | 2.4(1)  | -3.7157    | 0.43(2)  | 1.6&times;10<sup>-8</sup>
`md_lj_mts`&Dagger; | 0.002 (142) | 1.0040(2) | -3.5196(2) | 0.559(1) | 2.26(1) | -3.7153(2) | 0.412(1) | 1.1&times;10<sup>-7</sup>
`md_lj_mts`&sect; | 0.005 (142) | 1.017(2)  | -3.491(4)  | 0.610(7) | 2.26(1)   | -3.686(4)  | 0.463(7) | 6.8&times;10<sup>-6</sup>
`md_lj_mts`&para; | 0.005 (142) | 1.0094(8) | -3.508(2)  | 0.576(3) | 2.26(1)   | -3.703(2)  | 0.429(3) | 7.8&times;10<sup>-7</sup>

&dagger; All the timesteps the same length, as a check of the program book-keeping.

&Dagger; Program default parameters; note the smaller timestep.

&sect; Program defaults, except for the timestep.

&para; Identical to &sect; except that the switching length lambda is increased from 0.1 to 0.15.

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

Source                       | &rho;     | _T_   | _E_ (c)    | _P_ (c)  | _E_ (f)    | _P_ (f)  | _C_ (f)
------                       | -----     | ----- | -------    | -------  | -------    | -------  | --------
Thol et al (2016) (f)        | 0.75      | 1.00  | -3.3197    | 0.7008   | -3.7212    | 0.3996   |  2.2630  
`mc_nvt_lj`                  | 0.75      | 1.00  | -3.332(1)  | 0.651(3) | -3.734(1)  | 0.350(3) |  2.28(1)
`mc_nvt_lj_re`&sharp;        | 0.75      | 1.00  | -3.332(1)  | 0.648(2) | -3.734(1)  | 0.347(2) |  2.258(4)
`mc_nvt_lj_ll`&Dagger;       | 0.75      | 1.00  | -3.3230(3) | 0.669(1) | -3.7246(3) | 0.367(1) |  2.27(1)
`mc_npt_lj`&sect;            | 0.7501(2) | 1.00  | -3.331(1)  | 0.666(2) | -3.733(1)  | 0.364(2) |          
`mc_npt_lj_ll`&Dagger;&sect; | 0.7506(4) | 1.00  | -3.332(3)  | 0.660(3) | -3.734(3)  | 0.358(3) |          
`mc_zvt_lj`&para;            | 0.7504(4) | 1.00  | -3.333(3)  | 0.668(4) | -3.735(3)  | 0.366(4) |          
`mc_zvt_lj_ll`&Dagger;&para; | 0.7501(3) | 1.00  | -3.328(2)  | 0.669(2) | -3.729(2)  | 0.368(2) |          

&Dagger; Indicates a larger system size, _N_=864 (or approximately so for `mc_zvt_lj_ll`).
Note that the linked lists do not give an enormous speedup for this system size,
which corresponds to 4&times;4&times;4 cells.

&sect; The constant pressure simulations were run at _P_=0.69, the program default.
The measured _C<sub>p</sub>_ (full) values were 5.28(7) for `mc_npt_lj` and 5.04(16) for `mc_npt_lj_ll`,
compared with Thol et al (2016) EOS giving 5.22.
The `mc_npt_lj_ll` program was run with non-default value `db_max`=0.015 to give a volume acceptance ratio around 9%.

&para; The grand canonical programs were run at activity _z_=0.0795, the program default value.
The Thol et al (2016) LRC-corrected value to give &rho;=0.75 would be _z_=0.080627.
For `mc_zvt_lj` the box length was _L_=7&sigma;; for `mc_zvt_lj_ll` _L_=10.5&sigma;.
Acceptance rate of creation/destruction moves is quite small, at about 0.3%.
For other state points see below.

&sharp; The `mc_nvt_lj_re` program was run for four temperatures, see below for details.

Several of these programs could be improved to use array reallocation (available in Fortran)
to make them more resilient against changes in box size or number of particles.
For simplicity we have not included these features.

The measured pressures _P_ (c) are systematically a little low;
this is particularly noticeable for the constant-pressure programs,
where they might be expected to agree with the user-defined value of _P_.
This reflects the approximate nature of
the delta correction applied to the virial pressure,
to account for the discontinuous potential at _R<sub>c</sub>_.
At the density &rho;=0.75, with _R<sub>c</sub>_=2.5,
the pressure correction is &Delta; _P_&asymp;-0.3,
which is substantial.
However, this estimate is based on the assumption
that the pair distribution function _g(R<sub>c</sub>)_=1.
In fact, the choice _R<sub>c</sub>_=2.5 is a poor one in this regard,
lying near a local minimum where _g(R<sub>c</sub>)_&asymp; 0.91
(an illustration of _g(r)_ appears below in the __Pair distribution function__ section).
Consequently the applied correction is slightly too large,
and the resulting estimated pressure is systematically too low by &asymp; 0.03.
This serves as a reminder to always make clear what the cutoff is,
and what corrections (for discontinuities or long-range interactions)
have been applied.

In principle, there should be a delta correction for the configurational temperature.
Long-range corrections to _T_<sub>c</sub> are discussed by
A Baranyai _J Chem Phys,_ __112,__ 3964 (2000) and by
A Lervik, O Wilhelmsen, TT Trinh, HR Nagel, _J Chem Phys,_ __143,__ 114106 (2015).

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

### Gibbs Monte Carlo program
The program `mc_gibbs_lj` carries out Gibbs ensemble Monte Carlo,
and to test it we selected a temperature _T_=1.0,
which is below the critical point for the cut (but not shifted) LJ potential
(see tables above).
It was found convenient to start from a lower temperature,
with configurations at gas and liquid densities, with roughly equal numbers of particles,
and slowly work upwards in temperature, to equilibrate.
Note that the program expects two starting configurations: `cnf1.inp` and `cnf2.inp`.
The total number of atoms was fixed at _N_<sub>L</sub>+_N_<sub>G</sub>=512
and total volume _V_<sub>L</sub>+_V_<sub>G</sub>&asymp;5514.
Exchanges of box identity are expected as the critical temperature is approached,
and so one should not place blind trust in the separate box averages reported by the program,
but refer to histograms of density, energy etc.,
illustrative examples of which appear here.

![mc_gibbs_lj histograms](mc_gibbs_lj_his.png "mc_gibbs_lj histograms")

At _T_=1.0, however, these exchanges of box identity
are expected to be infrequent, were not observed in the test runs,
and the averages corresponded well to literature values for the coexistence parameters.
The production run corresponded to default parameters in the program.

Source  | &rho;<sub>L</sub> | &rho;<sub>G</sub> | _P_<sub>L</sub> | _P_<sub>G</sub> | _E_<sub>L</sub> (c) | _E_<sub>G</sub> (c)
-------              | -------- | -------- | -------  | -------- | --------------  | --------------
Trokhymchuk et al MC | 0.6542   | 0.0439   | 0.0336   | 0.0336   |                 |
Trokhymchuk et al MD | 0.6507   | 0.0500   | 0.0380   | 0.0380   | -2.713 &Dagger; | 1.047 &Dagger;
`mc_gibbs_lj`        | 0.653(1) | 0.050(1) | 0.031(2) | 0.038(1) | -2.731(5)       | 1.049(9)

&Dagger; Indicates values for given &rho; and _T_ from the Thol et al (2016) EOS (f) with cutoff correction.

The small discrepancy between measured pressures in the two phases reflects the approximate nature
of the delta correction for potential discontinuity, particularly in the liquid phase (see above).
For a density &rho;&asymp; 0.65 and _R<sub>c</sub>_=2.5
the pressure correction is &Delta; _P_&asymp;-0.23.
However, this assumes _g(R<sub>c</sub>)_=1,
whereas actually _g(R<sub>c</sub>)_&asymp; 0.95 at this density.
Hence the correction is too large by approximately 0.01.

### Replica exchange program
The `mc_nvt_lj_re` program uses MPI to handle communications between processes.
Here are some notes on the way the code is written.

We have only attempted to handle the most obvious errors at the start of the program,
such as missing configuration files and incorrect user data,
by closing down all the processes.
A production code would take more care to handle exceptions during the run.
Unhandled exceptions could possibly lead to processes hanging or becoming deadlocked,
so you should be aware of the dangers in running this example.

In the program, all processes write to their standard output `output_unit`, but the default in MPI is
for all this output to be collated (in an undefined order) and written to a single channel. Testing
was carried out using Open MPI, which allows the program to be run with a command line which includes
an option for each process to write to separate files, similar to the following:
```
mpirun -np 4 -output-filename out ./mc_nvt_lj_re < mc.inp
```
whereby the standard output files are named `out##`, the `##` part being determined by the process rank.
If your implementation does not have this option, you should edit the code to explicitly open a file for
standard output, with a process-rank-dependent name, and associate the `output_unit` with it.

The `mc_nvt_lj_re` program conducts runs at several temperatures: four were used in testing.
The default program values include _T_=1.0, which is reported above, and here is the complete set,
with expected values from the Thol et al (2016) equation of state (f) corrected for cutoff.
As usual the program employed the cut (but not shifted) potential.
All runs are for density &rho;=0.75, _N_=256, as usual.
At the lowest temperature, the full-potential system would lie in the coexistence region,
and the estimated pressure is negative.

Source                 | _T_    | _E_ (c)   | _P_ (c)  | _E_ (f)   | _P_ (f)   | _C<sub>v</sub>_ (f)
------                 | -----  | -------   | -------  | -------   | -------   | --------
Thol et al (2016) (f)  | 0.8772 | -3.6001   | 0.1942   | -4.0017   | -0.1070   |  2.3081  
`mc_nvt_lj_re`         | 0.8772 | -3.613(1) | 0.140(2) | -4.014(1) | -0.161(2) |  2.31(1)
Thol et al (2016) (f)  | 1.0000 | -3.3197   | 0.7008   | -3.7212   |  0.3996   |  2.2630  
`mc_nvt_lj_re`         | 1.0000 | -3.332(1) | 0.648(2) | -3.734(1) |  0.347(2) |  2.258(4)
Thol et al (2016) (f)  | 1.1400 | -3.0055   | 1.2571   | -3.4070   |  0.9559   |  2.2278  
`mc_nvt_lj_re`         | 1.1400 | -3.016(1) | 1.212(2) | -3.417(1) |  0.911(2) |  2.233(4)
Thol et al (2016) (f)  | 1.2996 | -2.6523   | 1.8667   | -3.0539   |  1.5655   |  2.1989  
`mc_nvt_lj_re`         | 1.2996 | -2.662(1) | 1.820(3) | -3.063(1) |  1.519(3) |  2.214(5)

The above (default) temperatures are chosen to give swap acceptance ratios all fairly close to 20% here
(of course, the set of temperatures, and all other run parameters, may be chosen by the user in a
namelist contained in the input file).
It should be noted that process `m` reports the swap acceptance ratio for exchanges with process `m+1`,
and the output file for the process with highest rank will report a zero swap ratio.

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
A system size _N_=256 was used
(for such a short-range potential, this system size was also suitable for the
link-cell program).
The given program defaults, including a time step of 0.005, were used throughout,
except for the strain rate which was varied.

Strain rate | _E_       | _P_      | &eta;    | _E_       | _P_       | &eta;
-----       | -----     | -----    | -----    | -----     | -----     | -----
0.04        | 1.8042(1) | 6.390(1) | 2.38(4)  | 1.8039(3) | 6.389(2)  | 2.31(4)
0.16        | 1.8095(2) | 6.426(1) | 2.228(8) | 1.8098(2) | 6.427(1)  | 2.23(1)
0.64        | 1.8648(2) | 6.777(2) | 1.938(2) | 1.8646(2) | 6.776(1)  | 1.935(2)

In the table above, for each strain rate,
the results in columns 2-4 come from `md_nvt_lj_le`
and those in columns 5-7 from `md_nvt_lj_llle`
(essentially identical, but roughly twice as fast for _N_=256).
In all cases the kinetic energy was conserved very accurately by the algorithm.
The results, particularly the increase in _E_ and _P_,
and the decrease in shear viscosity &eta;,
as the strain rate increases,
are in good agreement with the above papers.
Incidentally, at the highest strain rate 0.64,
the configurational temperature is systematically about 1% lower
than the (constrained) kinetic temperature.

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
----- | -----     | -----          | -----           | -----
0.50  | 1.6347    | 1.634(2)       | 1.632(1)        | 0.502(2)
0.55  | 2.0574    | 2.051(3)       | 2.055(1)        | 0.553(3)
0.60  | 2.5769    | 2.573(4)       | 2.573(1)        | 0.600(3)
0.65  | 3.2171    | 3.210(8)       | 3.215(2)        | 0.651(3)
0.70  | 4.0087    | 3.996(8)       | 4.005(3)        | 0.700(2)
0.75  | 4.9910    | 4.960(7)       | 4.985(4)        | 0.749(2)

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
run `initialize` with `molecules="linear", lattice=.false.`,
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
which is sufficient for our purposes.
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
(prepared using `initialize` with `molecules="chain"` to give random non-overlapping atom positions)
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
The default run lengths are fairly modest here: 10 blocks,
each consisting of 100000 steps of length &delta;t=0.002.
The primary indicator of a correctly-functioning program is energy conservation,
and this was checked in all cases.
Energies were chosen to give average temperatures close to the values used in
the MC simulations above.

Results for constrained system (columns 2:4 RATTLE, columns 5:7 MILC-SHAKE):

_E_     | _T_       | _R_<sub>g</sub> | _C<sub>v</sub>_ | _T_       | _R_<sub>g</sub> | _C<sub>v</sub>_
-----   | -----     | -----           | -----           | -----     | -----           | -----
-2.0246 | 0.2485(2) | 1.06374(4)      | 2.176(6)        | 0.2475(1) | 1.06450(3)      | 2.172(7)
-1.9145 | 0.296(2)  | 1.073(1)        | 2.38(8)         | 0.2989(3) | 1.0724(1)       | 2.27(2)
-1.6145 | 0.345(4)  | 1.125(2)        | 3.14(8)         | 0.347(2)  | 1.125(1)        | 3.16(6)
-1.3495 | 0.404(1)  | 1.182(2)        | 2.39(1)         | 0.404(2)  | 1.183(2)        | 2.50(2)
-1.2195 | 0.451(1)  | 1.207(1)        | 2.36(2)         | 0.449(2)  | 1.210(1)        | 2.34(2)
-1.0968 | 0.499(2)  | 1.234(1)        | 2.28(1)         | 0.503(2)  | 1.231(2)        | 2.28(2)
-0.1244 | 1.009(5)  | 1.471(5)        | 2.04(2)         | 1.024(6)  | 1.455(7)        | 2.01(2)
 1.0456 | 2.008(5)  | 1.754(9)        | 1.653(3)        | 2.009(7)  | 1.753(8)        | 1.652(3)
 3.6459 | 4.996(4)  | 1.889(7)        | 1.534(1)        | 4.988(3)  | 1.901(6)        | 1.534(1)

Results for _k_<sub>spring</sub>=10000 system using MTS:

_E_      | _T_       | _R_<sub>g</sub> | _C<sub>v</sub>_
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

Results for _k_<sub>spring</sub>=400 system using MTS:

_E_     | _T_       | _R_<sub>g</sub> | _C<sub>v</sub>_
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

For the hard-sphere square-well chain, the aim was to show the operation of the Wang-Landau method.
In `mc_chain_wl_sw` we use pivot and crankshaft moves as well as CBMC regrowth.
In a practical application it would be advisable to include some bridging moves as well.
Reasonably long chains, _N_=128, have been studied by this method,
and exact results are available for very short chains;
see, for example,

* MP Taylor,  _J Chem Phys,_ __118,__ 883 (2003),
* JE Magee, L Lue, RA Curtis, _Phys Rev E,_ __78,__ 031803 (2008),
* MP Taylor, W Paul, K Binder, _J Chem Phys,_ __131,__ 114907 (2009),

who provide references to other simulation work.

For testing purposes our aims are quite modest:
we choose _N_=6, bond length equal to &sigma;, and a nonbonded interaction range of 1.5&sigma;.
The starting chain configuration can be prepared using `initialize` in the usual way
(note the non-default value of the bond length).
Default parameters are used in `mc_chain_wl_sw`,
including a flatness criterion of 0.9.
The entropy modification constant `ds` is halved at each stage,
and there are 20 stages.
For this system, the energy range (in units of the well depth) is
_E_ = 0 &hellip; -10.
The principal result is the histogram of entropies _S(E)_ produced at the final stage.
For convenience we (arbitrarily) define _S_(0)=0.
We conduct a set of nine independent WL runs,
and report the results from the two runs with the highest and lowest values of _S_(-10),
which bracket all the other results in the set,
as a rough indication of the errors.
We compare with the exact values calculated from the density of states
of Taylor (2003), normalized in the same way to make _S_(0)=0.

_E_    | _S(E)_ (exact) | _S(E)_ (WL) | _S(E)_ (WL)
------ | ------         | ------      | ------
  0.0  |   0.0000       |   0.0000    |   0.0000
 -1.0  |   0.7521       |   0.7629    |   0.7518
 -2.0  |   0.6661       |   0.7014    |   0.6683
 -3.0  |   0.2108       |   0.2308    |   0.2108
 -4.0  |  -0.4433       |  -0.4152    |  -0.4449
 -5.0  |  -1.3484       |  -1.3316    |  -1.3444
 -6.0  |  -2.4438       |  -2.4256    |  -2.4322
 -7.0  |  -3.6832       |  -3.6634    |  -3.6733
 -8.0  |  -5.8548       |  -5.8440    |  -5.8620
 -9.0  |  -8.4766       |  -8.4050    |  -8.4733
 -10.0 | -14.9981       | -14.6824    | -15.0295

As a further check, we ran a set of canonical ensemble calculations for the same system
with `mc_chain_nvt_sw` at selected temperatures.
The program default is to run for 10 blocks, each of 100000 steps;
this was increased to 10 blocks of 500000 steps for temperatures
below 0.25.
The results may be compared with values reconstructed using the
`wl_hist` program from the simulation histograms.
Below we show the heat capacity per atom from the above two WL runs (red),
from the exact density of states of Taylor (black),
and from the canonical ensemble calculations (blue error bars).

![Wang-Landau test results](wl.png "Wang-Landau test results")

It is also straightforward to compare average energies and radii of gyration,
but we do not do that here.

As the chain length increases, the energy landscape becomes more challenging.
For _N_=13, with the same bond length of &sigma;,
and nonbonded interaction range of 1.5&sigma;,
sensible results may still be achieved
with the simple example program `mc_chain_wl_sw`.
Once more, as a reference for comparison,
we ran a set of canonical ensemble calculations with `mc_chain_nvt_sw`,
using the same parameters described above.
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
The program was run with default parameters,
except that the flatness criterion was set at 80%.
The results are from the histograms produced in the 20th stage.
This analysis can also be performed (for any desired temperature) by the program `wl_hist`, after the run.
The results are generally in good agreement with the canonical ensemble test runs.
The most significant discrepancies are in the heat capacities at the lowest two temperatures,
which reflects the poor sampling of the canonical ensemble program (with these basic MC moves),
and the sampling problems about to be discussed.

This particular test run illustrated one drawback of the simplest Wang-Landau implementation:
two low-lying energies (corresponding to _q_=38 and 39 square-well interactions) were discovered
during the very last stage, in which `ds` is very small.
Accordingly, the system remained stuck in these low-energy states for a very long time,
until their tabulated entropy _S(q)_ reached a high enough value to allow _q_ to change;
even then, the final weight of the lowest state in the final "flat" histogram
could not be considered completely reliable.
The overall run length was of order 90000 blocks of 10000 steps each,
most of it spent in stage 20.

Repeating the run with the same parameters typically produces different results,
depending on whether, and when, these low-lying states are discovered.
This affects the canonical ensemble results reconstructed from the histograms
through `wl_hist`, particularly at the lower temperatures,
while the higher temperatures are largely unaffected.
From a single run, or a few runs, there might be no indication of anything wrong.
It is always a danger with any Monte Carlo method, including Wang-Landau,
that inaccessible configurations will not be sampled.
Various improvements of the method may be found in the literature.
(The reader is invited to work out the structure of the 13-atom chain configuration
with hard sphere diameter and bond length both equal to &sigma;
having 39 nonbonded interactions within range 1.5&sigma;.
_Hint:_ it is not one of the highest-symmetry clusters such as the icosahedron;
however, it is based on a fairly symmetric cluster,
with small distortions due to the fixed bond length and the need to maximise interactions.)

To obtain more reliable results it would be advisable to add a
production run at the end of the `mc_chain_wl_sw` program,
during which the weights are no longer adjusted,
allowing averages to be generated using a proper Markov chain.

## Polyatomic Lennard-Jones program
The program `mc_nvt_poly_lj` conducts Monte Carlo simulations of a system of rigid molecules
composed of Lennard-Jones interaction sites.
For simplicity the sites are taken to be identical, although the program is easily generalized.
Molecular orientations are represented by quaternions,
which are used to calculate the rotation matrix
and hence the interaction site positions.

We test this with the three-site model of orthoterphenyl, a fragile glassformer,
described in the following publications amongst others.

* LJ Lewis, G Wahnstrom, _Sol State Commun,_ __86,__ 295 (1993)
* LJ Lewis, G Wahnstrom, _Phys Rev E,_ __50,__ 3865 (1994)
* S Mossa, E La Nave, HE Stanley, C Donati, F Sciortino, P Tartaglia, _Phys Rev E,_ __65,__ 041205 (2002)
* E La Nave, S Mossa, F Sciortino, P Tartaglia, _J Chem Phys,_ __120,__ 6128 (2004)

We compare with the results of Mossa et al (2002).
The sites are arranged at the vertices of an isosceles triangle with bond angle 75 degrees,
LJ parameters &epsilon; = 5.276 kJ mol<sup>-1</sup>,
&sigma;=0.483nm,
and two equal bonds of length &sigma;.
The program employs the usual reduced units based on &epsilon; and &sigma;
and in these units the potential cutoff of Mossa et al (2002) is _R_<sub>c</sub>=2.612;
the pair potential is Lennard-Jones with a shifted-force correction term, linear in _r_,
to make the potential and its derivative vanish at _r_=_R_<sub>c</sub>.
Apart from the temperatures,
default program parameters were used throughout the tests.

Tests were performed at &rho;=0.32655 which is equivalent to &rho;<sub>4</sub>=1.108g cm<sup>-3</sup>
in Mossa et al (2002).
Comparisons of potential energy (_PE_=_E_-3 _T_ converted to kJ/mol with a factor 5.276)
were made with the fit given by eqn (23) of that paper.
Note that &epsilon;/k<sub>B</sub>&asymp;635 K.

_T_   | _E_        | _P_       | _T_ (K) | _PE_ (kJ/mol) | _PE_ (kJ/mol) eqn (23)
----- | -----      | -----     | -----   | -----         | -----
0.5   | -12.993(3) |  1.773(8) |  317    | -76.47(2)     | -76.945
1.0   | -9.817(9)  |  5.86(2)  |  635    | -67.62(5)     | -67.634
1.5   | -6.84(1)   |  9.38(4)  |  952    | -59.83(5)     | -60.011
2.0   | -4.07(1)   | 12.37(4)  | 1270    | -53.13(5)     | -53.265

A second set of tests was performed at _T_=0.6&asymp;380K
at the specified densities &rho;<sub>1</sub>, &hellip; &rho;<sub>5</sub> of Mossa et al (2002).
A set of starting configurations is provided in the [Data repository](https://github.com/Allen-Tildesley/data).
Here the excess pressure (_P_(ex)=_P_-&rho;_T_ converted to MPa
with a factor 77.75 based on the values of &epsilon; and &sigma;)
is compared with the fit given by eqn (28) and the coefficients in Table III of Mossa et al (2002).
NB the volumes to insert into the equation are those of their Table I,
which are specific to their system size.
In addition their eqn (29) with coefficients in Table V is a fit to their potential energy,
which we calculate from the simulation as described above.

Id    | &rho;   | _E_        | _P_     | _P_(ex) (MPa) | _P_(ex) (MPa) eqn (28) | _PE_ (kJ/mol) | _PE_ (kJ/mol) eqn (29)
----- | -----   | -----      | -----   | -----         | -----                  | ------        | ------
1     | 0.30533 | -11.625(4) | 0.45(1) | 20.7(7)       | 19.077                 | -70.83(2)     | -70.818
2     | 0.31240 | -11.914(6) | 1.01(2) | 64(2)         | 60.143                 | -72.36(3)     | -72.289
3     | 0.31918 | -12.143(5) | 1.74(1) | 120(2)        | 112.798                | -73.56(3)     | -73.601
4     | 0.32655 | -12.400(4) | 2.53(1) | 181(1)        | 177.222                | -74.92(2)     | -74.886
5     | 0.33451 | -12.487(4) | 3.84(1) | 283(1)        | 253.510                | -75.38(2)     | -75.825

In making these comparisons,
our default run length (10 blocks of 20000 sweeps each) should be borne in mind,
since this system can show sluggish behaviour.
The MD simulations of Mossa et al (2002) are reported to extend to several hundred nanoseconds
(of order 10<sup>7</sup> MD timesteps) at the lowest temperatures.

For comparison we provide a molecular dynamics code `md_nvt_poly_lj` for the same model.
The program takes the molecular mass _M_ to be unity.
Mossa et al (2002) ascribe a notional mass of 78u to each of the three LJ sites,
so _M_&asymp;3.9&times;10<sup>-25</sup>kg.
Combined with the above values of &epsilon; and &sigma;,
this gives a time scale (_M_/&epsilon;)<sup>1/2</sup>&sigma; &asymp; 3.22 ps.
The timestep of &delta;t=0.01 ps used by Mossa et al (2002)
corresponds to the default value in the program `dt=0.003` in these units.
By default, the program simulates the constant-_NVE_ ensemble,
but there is an option to simulate at constant _NVT_ by velocity randomization (Andersen thermostat).
If the latter option is selected,
the program will read configurations in the same format as `mc_nvt_poly_lj` (positions and quaternions only),
selecting random initial velocities and angular momenta,
which can be convenient.

By default the program calculates the inertia tensor from the LJ site bond vectors,
assuming equal masses.
For simplicity it is assumed that the bond vectors are defined such that
the principal axes of the inertia tensor coincide with
the xyz axes of the molecular coordinate system,
with the centre of mass at the origin;
it is always possible to arrange this.
In general,
the three principal moments of inertia will all be different,
so the molecule is an asymmetric top.
The MD algorithm for rotation is a symplectic one
in which a `kick` propagator advances the space-fixed angular momenta,
using the torque on each molecule,
and a succession of `drift` steps implement free rotation about each of the principal axes.
This is described in the text, section 3.3; see

* A Dullweber, B Leimkuhler, R McLachlan, _J Chem Phys,_ __107,__ 5840 (1997),
* TF Miller, M Eleftheriou, P Pattnaik, A Ndirango, D Newns, GJ Martyna, _J Chem Phys,_ __116,__ 8649 (2002).

The results below are for test runs in both constant-_NVE_  and constant-_NVT_ ensembles,
at (approximately) the same state points as those given above.
All runs were 10&times;20000 steps in length and used program defaults,
except for `t_interval=1` and the specified temperature in the _NVT_ case.
For constant-_NVE_ runs we report RMS energy fluctuations,
and _T_ is the average translational temperature.

 &rho;   | _T_       | _E_        | _P_       | _E_(RMS)
 -----   | -----     | -----      | -----     | -----
 0.32655 | 0.5       | -12.984(5) |  1.81(1)  |
 0.32655 | 0.5082(4) | -12.9838   |  1.755(6) | 1.23&times;10<sup>-8</sup>
 0.32655 | 1.0       | -9.80(2)   |  5.87(4)  |
 0.32655 | 1.004(1)  | -9.8006    |  5.896(5) | 9.88&times;10<sup>-8</sup>
 0.32655 | 1.5       | -6.83(1)   |  9.35(3)  |
 0.32655 | 1.506(1)  | -6.8326    |  9.378(6) | 3.67&times;10<sup>-7</sup>
 0.32655 | 2.0       | -4.05(1)   | 12.42(4)  |
 0.32655 | 2.007(1)  | -4.0507    | 12.405(4) | 9.39&times;10<sup>-7</sup>

 &rho;   |  _T_       | _E_        | _P_      | _E_ (RMS)
 -----   | -----      | -----      | -----    | -----
 0.30533 |  0.6       | -11.613(5) | 0.45(2)  |
 0.30533 |  0.6013(8) | -11.6131   | 0.485(4) | 1.31&times;10<sup>-8</sup>
 0.31240 |  0.6       | -11.892(6) | 1.04(2)  |
 0.31240 |  0.6039(8) | -11.8923   | 1.058(5) | 1.52&times;10<sup>-8</sup>
 0.31918 |  0.6       | -12.147(4) | 1.75(1)  |
 0.31918 |  0.603(1)  | -12.1465   | 1.710(9) | 1.73&times;10<sup>-8</sup>
 0.32655 |  0.6       | -12.362(3) | 2.61(1)  |
 0.32655 |  0.604(2)  | -12.3616   | 2.58(1)  | 2.05&times;10<sup>-8</sup>
 0.33451 |  0.6       | -12.453(7) | 3.96(2)  |
 0.33451 |  0.612(1)  | -12.4532   | 3.87(1)  | 2.53&times;10<sup>-8</sup>

## DPD program
For the `dpd` example, we recommend generating an initial configuration
using the `initialize` program, with namelist input similar to the following
```
&nml n = 100, density = 3.0, lattice = .false.,
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

The force between the molecules is calculated from the analytical derivative of the
T-tensor with respect to the separation vector.
This naturally leads to formulae where the original T-tensor of rank n
is replaced by one of rank n+1.

The torque on each molecule is calculated by formulae similar to those used for
torques on multipoles in an external field, field gradient, etc., but in which the
field terms are replaced by tensors based on T and the multipoles on the other molecule.
This naturally leads to formulae involving the Levi-Civita (antisymmetric) symbol.

In practical applications, the formulae would usually be incorporated in a scheme
for handling long-range forces in periodic boundaries (e.g. Ewald sum).

## Ewald program
The k-space and r-space contributions to the Ewald sum are illustrated in `ewald_module`
and we provide a program `ewald` to test these.
The program reads in a configuration file `cnf.inp` in the usual format:
any of the Lennard-Jones or hard-sphere configurations would be suitable.
Charges are assigned to the atoms in an arbitrary way.
The program itself adds the surface term (the self term is included in the k-space routine).
Then, a comparison is made with the brute force summation over all pairs
in shells of periodic boxes surrounding the central box.
For default parameters, and test configurations with _N_=256,
reasonable convergence is obtained within 8-10 shells.
One can adjust the screening parameter kappa within reason
(and the number of k-vectors may need changing as well):
the contributions of r-space and k-space terms will change, but their sum should
remain approximately constant.

There is also a comparison with a simplified particle-mesh Ewald method.
As discussed in the text, the charge distribution is assigned to a cubic mesh,
Fourier transformed by FFT, and used to calculate the total potential energy,
using the solution of Poisson's equation in Fourier space.
In doing so, accuracy is improved by optimizing the so-called influence function G.
In this example, we use a simple sharpening function discussed by

* V Ballenegger, JJ Cerda, C Holm, _J Chem Theo Comp,_ __8,__ 936 (2012)

but more sophisticated optimized functions are possible. It is easy to comment out
this sharpening function, to see the extent of the correction; it is reasonably
significant for the default parameter values.

See below for more discussion of the mesh function, provided in `mesh_module`,
and of the FFT routine which is illustrated in `fft3dwrap`.

## Mesh program
The program `mesh` generates a random configuration of a small number of charges
and illustrates the way this may be assigned to a regular cubic mesh using the
triangular-shaped cloud distribution described in

* RW Hockney, JW Eastwood, _Computer simulation using particles_ (Adam Hilger, Bristol, 1988)

The function generating the charge density is provided in `mesh_module`. The mesh dimension
is, by default, kept small enough to print out the whole array for inspection afterwards.
The number of charges and mesh dimension may be adjusted by the user, via namelist parameters.

## Cluster program
The `cluster` program is self contained. It reads in a configuration of atomic positions
and produces a circular linked list for each cluster identified within the configuration.
The best value of the critical separation `r_cl` depends on the particular physical system
being considered. To illustrate, we have provided a file `cluster.inp` consisting of
_N_=256 atoms at low overall density, generated by a short quench from a disordered system
at high temperature to a low temperature inhomogeneous state (not yet at equilibrium).
This file should be copied into the working directory before running the program.
With the default value `r_cl`=1.5, six well-separated clusters should be identified.
The results in this case are moderately insensitive to the value of `r_cl`, but increasing
it above 3 includes all atoms in a single cluster, while reducing it below 1.15 will start to
separate isolated atoms into clusters of their own.

Clustering algorithms are part of the standard toolkit of data analysis, and in practical
applications it may be more efficient and convenient to use a packaged implementation of
an algorithm such as `dbscan`  

* M Ester, H-P Kriegel, J Sander, X Xu. (1996).
[Proc. Second Int. Conf. on Knowledge Discovery and Data Mining (KDD-96) p 226](https://www.aaai.org/Papers/KDD/1996/KDD96-037.pdf)
(Eds: E Simoudis, J Han, UM Fayyad; AAAI Press, 1996).

Fortran implementations of `dbscan` are available from various sources.
For systems in periodic boundaries, rather than supplying the atomic positions, the user should
compute a distance matrix using the minimum image convention, and supply that to the routine,
as suggested by [Turci](https://francescoturci.wordpress.com/2016/03/16/clustering-and-periodic-boundaries/)
in the context of a Python version.

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

![corfun test results](corfun.png "corfun test results")

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

![diffusion test results](diffusion.png "diffusion test results")

## Pair distribution function
The program `pair_distribution` reads in a set of configurations and calculates
the pair correlation function _g(r)_.
We limit the number of configurations to a maximum of 1000 (numbered from 000 to 999)
simply so as to use a fixed naming scheme for the input configurations;
in a practical application, a trajectory file would be used instead.
We have tested it on a set of 500 configurations
of _N_=256 Lennard-Jones atoms,
cut (but not shifted) at _R_<sub>c</sub>=2.5&sigma;,
at the usual state point &rho;=0.75, _T_=1.0.
The interval between configurations was 100 MC sweeps.
This data set is provided in the
file `pair_distribution_data.zip` in the [Data repository](https://github.com/Allen-Tildesley/data).
Using the default resolution of 0.02&sigma;,
the results shown below were obtained for _g(r)_.

![g(r) test results](gr.png "g(r) test results")

## Interface pair correlation function
The program `grint.f90` reads in a set of configurations and calculates
the pair correlation function for a system that is inhomogeneous in the
z direction. It is assumed that the configurations consist of a liquid slab,
surrounded by gas, with the interfaces lying in the _xy_-plane.
The two interfaces are located by fitting the instantaneous density profile
to a difference of two tanh functions. Then the single-particle density function,
relative to each of the interface locations, is calculated.
To make the process as robust as possible, an initial guess at the mid-point
of the liquid slab should be provided, and this is updated automatically as
successive configurations are read in, so as to shift the liquid slab into
the middle of the periodic box, before fitting.
Also, the results of one fit are passed on as the starting point of the next one.
The program handles cubic boxes only;
the modifications necessary to handle non-cubic boxes are fairly easy to make,
but we will not do that here.
No attempt is made to correct for capillary-wave fluctuations affecting the
width of the interface.

Having located, and combined, the two interfaces, the histograms for calculating
the one-body and two-body densities are accumulated, in a coordinate system
which has its origin at the interface. The one-body density is then fitted by
a single tanh curve: this could easily be changed by the user if desired.
For simplicity, in the normalization of the two-body density to give the
pair correlation function, we use the _fitted_ tanh form of the single particle
density. This is obviously an approximation, and could be replaced by a better fit,
or an interpolation scheme, to evaluate the function at arbitrary z.
Finally, the program writes out the average, overall, density profile, the
single-particle density (with fit) and slices at selected values of z and c
through the pair correlation function.

In testing this program, it is important to use a large enough system so that
all values of z<sub>1</sub> and z<sub>2</sub> of interest (measured relative to
the interface position) lie far from the _other_ interface position.

The program was tested on a system of _N_=10000 atoms, interacting through
the Lennard-Jones potential cut (but not shifted) at _R_<sub>c</sub>=2.5&sigma;,
in a cubic box of side 30&sigma;, at a temperature _T_=0.90.
For this system, &rho;<sub>G</sub> &asymp; 0.024, &rho;<sub>L</sub> &asymp; 0.713
(see Trokhymchuk, op. cit.). A set of 100 configurations from this run, together
with the output of `grint.f90` with default parameters, are provided in the
file `grint_data.zip` in the [Data repository](https://github.com/Allen-Tildesley/data).

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
and compare with the analytically, exactly, known result.
User input defines the number of grid points and the box size;
sensible defaults are provided.
The library that we use for this example is [FFTW](http://www.fftw.org/).

## Hit-and-miss and sample-mean
The two programs `hit_and_miss` and `sample_mean` illustrate two very simple
Monte Carlo methods to estimate the volume of a 3D object.
They are both described in detail at the start of Chapter 4.
No user input is required.
For the built-in values defining the geometry, the exact result is 5/3.

## Quantum simulation programs
The program `qmc_walk_sho` solves the diffusion equation in imaginary time
corresponding to the Schrodinger equation,
for a single simple harmonic oscillator.
Atomic units are chosen so that the effective diffusion coefficient is _D_=1/2.
A few hundred independent systems, or walkers, are simulated using a simple random walk
and a crude creation/destruction scheme based on the difference between the potential energy
and the trial energy.
The scheme is described in the text.
The value of the trial energy `et` is updated regularly,
and the hope is that, after convergence,
it will be equal to the correct ground-state energy for the system which, in this case, is 1/2.
The updating scheme, and several of the default parameters,
are taken from the following paper

* I Kostin, B Faber, K Schulten, _Amer J Phys,_ __64,__ 633 (1996).

Reasonable results for the energy and the ground-state wavefunction,
which is accumulated as a histogram of walker positions,
should be obtained using the default input values,
with an empty namelist `&nml /`;
these defaults include setting `et` initially to the exact ground state energy.
Other values such as `&nml et=0.6 /` may be supplied through the namelist in the usual way.
This type of simulation is sensitive to the initial value,
and quite noisy:
possible improvements are discussed in general terms in the text.

The program `qmc_pi_sho` carries out a path integral Monte Carlo simulation
for a single simple harmonic oscillator,
at a specified temperature and ring-polymer size _P_.
Larger values of _P_ give energies closer to the exact quantum mechanical canonical ensemble average.
For this simple model,
exact results can also be calculated for the finite values of _P_ used in the simulation

* KS Schweizer, RM Stratt, D Chandler, PG Wolynes, _J Chem Phys,_ __75,__ 1347 (1981),
* M Takahashi, M Imada, _J Phys Soc Japan,_ __53,__ 3765 (1984),

and a routine to evaluate these is included in the example.
No special techniques are used to accelerate the simulation;
standard Metropolis moves are employed.
Default parameters correspond to _P_=8, _T_=0.2.
The table below is for test runs at various values of _P_,
keeping the same temperature,
which generates a range of average energies between
the classical limit _E_=0.2
and the quantum limit _E_=0.506784;
in each case we compare with the exactly-known value for the same _P_.

_P_  | _E_ (MC)  | _E_ (exact)
---- | --------  | -----------
 2   | 0.3218(3) | 0.321951
 3   | 0.3933(4) | 0.392308
 4   | 0.4312(3) | 0.431618
 5   | 0.4543(4) | 0.454545
 6   | 0.4694(6) | 0.468708
 7   | 0.4778(6) | 0.477941
 8   | 0.4846(9) | 0.484244

The program `qmc_pi_lj` applies the path-integral method to the Lennard-Jones fluid.
The simplest, primitive, algorithm is used,
together with the crudest estimators for energy and pressure.
The program uses
single-particle Monte Carlo moves for the individual beads in the ring polymer,
along with translations of the centre-of-mass of each polymer.
As mentioned in the text, there are many improvements of all these aspects
of the algorithm, which are recommended for production work.

The program takes in configuration files `cnf##.inp` where the `##` reflects
the polymer bead number, in the range 1 to _P_.
These files have the same format as the classical Lennard-Jones configurations.
They may be prepared in the same way as usual,
from the `initialize` program, from a run of `mc_nvt_lj`
at the appropriate density and temperature,
or from a run of `qmc_pi_lj` for a different value of _P_.
It does no harm if these starting configurations are simply duplicates of each other,
provided that a preliminary run is conducted to allow the polymers to equilibrate,
after which all the output files `cnf##.out` may be renamed to `cnf##.inp`.

For testing, we compare with a set of simulations of neon,

* M Neumann, M Zoppi, _Phys Rev E,_ __65,__ 031203 (2002),

which are mainly based on an empirical pair potential,
but include selected results for Lennard-Jones for the case _P_=32.
The LJ parameters for neon are &epsilon;=36.8K, &sigma;=0.2789nm, atomic mass _m_=20.18u,
and hence a reduced de Boer parameter &lambda;=0.092&sigma;,
which is the default value in the program.
We choose their lowest-temperature state point,
(_T_,&rho;)=(25.8K,36.28nm<sup>-3</sup>)=(0.701087,0.787069) in reduced LJ units.
We use _N_=108 atoms, compared with Neumann and Zoppi's _N_=256,
and our runs are five times shorter (10 blocks of 10000 sweeps);
these differences should not be critical.
The maximum centre-of-mass displacement is 0.1&sigma;.
Because the intra-polymer springs increase in strength with _P_,
the maximum displacement parameter for individual bead moves
is reduced from 0.06&sigma; for _P_=2
down to 0.02&sigma; for _P_=32.
These choices give reasonable acceptance ratios for both kinds of move.
In the table below we report:
the rms radius _R_ of the ring polymers,
the average spring potential _E_(spring) per atom,
the kinetic energy _KE_ per atom,
total energy _E_ per atom, and pressure _p_
(the last two including the standard long-range correction
for the applied cutoff _R_<sub>c</sub>=2.5&sigma;).

_P_         | _R_        | _E_(spring) |  _KE_      |  _E_      |  _p_
------      | ------     | ------      | ------     | ------    | ------
1 &dagger;  | 0.0        |  0.0        | 1.052      | -4.692(1) | -0.756(5)
2           | 0.04043(1) |  0.8963(7)  | 1.2070(7)  | -4.454(1) | -0.408(6)
4           | 0.04942(1) |  2.906(1)   | 1.301(1)   | -4.320(3) | -0.231(9)
8           | 0.05181(2) |  7.073(3)   | 1.340(3)   | -4.261(3) | -0.150(8)
16          | 0.05247(2) | 15.479(4)   | 1.347(4)   | -4.252(4) | -0.148(5)
32          | 0.05261(4) | 32.287(8)   | 1.365(8)   | -4.233(8) | -0.140(6)
32 &Dagger; | 0.053      | 32.3008     | 1.352      | -4.227    | -0.039

&dagger; For completeness the _P_=1 runs were performed using `mc_nvt_lj`.

&Dagger; These results are taken from Table I of Neumann and Zoppi (2002),
converted into LJ reduced units.

A drawback of this state point is that the pressure is negative,
suggesting instability with respect to phase separation.
Nonetheless we have seen no evidence of crystalline structure
in the simulated configurations,
and assume that the liquid phase is metastable.
