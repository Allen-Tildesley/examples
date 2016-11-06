#Brief Guide
Here are some notes to assist in running the programs.

##Data Input
Most of the Fortran codes use a `namelist` to input a few parameters from standard input.
This gives an easy way to specify default values in the program itself, and to use a 
keyword-based syntax to specify values different from the default ones at run-time.
The input file, or string, should usually begin with `&nml` and end with `/`.
As a minimum, the program will expect to read `&nml /` (an empty list), but to
change the parameters, typical input might be 
```
&nml nblock=20, nstep=1000, dt=0.001 /
```
and the `key=value` pairs may be set out on different lines if you wish.

##Initial Configuration
Simulation runs require a starting configuration which can usually be prepared using
the `initialize` program (in `build_initialize/`).
The default parameters produce an FCC configuration of 256 atoms at reduced density &rho; = 0.75,
writing out just the positions (for an MC program) to a file `cnf.inp`.
If the parameter `velocities=.true.` is supplied, then positions and velocities are
written to the file, corresponding to a reduced temperature _T_ = 1.0.
The file may then be copied to the directory in which the run is carried out.
Typically, runs produce a final configuration `cnf.out`
(which may be renamed to `cnf.inp` as a starting point for further runs)
and intermediate configurations `cnf.001`, `cnf.002` etc during the run.

##State points for different Lennard-Jones models
For most of the examples, we use a cutoff of _Rc_ = 2.5 &sigma;.
For MC programs the cut (but not shifted) potential (c) is used in the simulation,
while for MD programs the cut-and-shifted potential (cs) is used.
In all cases, an attempt is made to produce results without, and with, long-range corrections,
the latter giving an estimate for the full potential (f).
Differences between the different models are discussed in various places,
see e.g. 

* A Trokhymchuk, J Alejandre, _J Chem Phys,_ __111,__ 8510 (1999).

Using their table V as a guide, we take the critical point to be roughly located at:

LJ model                 | _T_ | &rho; | _P_
-----                    | ---- | ---- | ----
f: full, with LRC        | 1.31 | 0.31 | 0.13
c: cut (but not shifted) | 1.19 | 0.32 | 0.11
cs: cut-and-shifted      | 1.08 | 0.32 | 0.09

At any temperature below _Tc_, the liquid state is bounded below by the
liquid-gas coexistence density, and using Tables II-IV of the same reference as a guide,
we take the values for three example temperatures as

LJ model                 |  _T_ | &rho; | _P_
-----                    | ---- |  ---- | ----
f: full, with LRC        |  0.8 | 0.793 | 0.005
c: cut (but not shifted) |  0.8 | 0.765 | 0.008
cs: cut-and-shifted      |  0.8 | 0.730 | 0.013

LJ model                 |  _T_ | &rho; | _P_
-----                    | ---- |  ---- | ----
f: full, with LRC        |  0.9 | 0.746 | 0.012
c: cut (but not shifted) |  0.9 | 0.714 | 0.020
cs: cut-and-shifted      |  0.9 | 0.665 | 0.030

LJ model                 |  _T_ | &rho; | _P_
-----                    | ---- |  ---- | ----
f: full, with LRC        |  1.0 | 0.695 | 0.026
c: cut (but not shifted) |  1.0 | 0.652 | 0.036
cs: cut-and-shifted      |  1.0 | 0.578 | 0.062

For Lennard-Jones, the default state point (&rho; = 0.75 and _T_ =1.0) 
lies in the liquid region of the phase diagram for all these variants of the model. 

For comparison with MD runs, a wide range of simulation data for LJ with 
_Rc_ = 2.5  &sigma; (cut-and-shifted, denoted cs below) is summarized by

* M Thol, G Rutkai, R Span, J Vrabec and R Lustig, _Int J Thermophys,_ __36,__ 25 (2015)

who present data (in Supplementary Information) and a fitted approximate 
equation of state (which they also provide as a C++ program). 

Approximate equations of state for the full LJ potential are available from

* J Kolafa and I Nezbeda, _Fluid Phase Equil,_ __100,__ 1 (1994) 
* M Mecke, A Muller, J Winkelmann, J Vrabec, J Fischer, R Span, W Wagner,
_Int J Thermophys,_ __17,__ 391 (1996) and __19,__ 1493(E) (1998)

Below we use the Kolafa & Nezbeda formulae (denoted f) and estimate the results
for the cut-but-not-shifted potential (denoted c) using the same corrections as in the MC codes.

Here we compare with typical test runs from our programs using default parameters except where stated.
Note that E/N is the total internal energy per atom, including the ideal gas part.

Source           | &rho;   | _T_         | _E/N_ (cs) | _P_ (cs) | _E/N_ (c) | _P_ (c) | _E/N_ (f) | _P_ (f) |
------           | ------- | ----------- | --------   | ------   | -------   | ------  | --------  | ------- |
Thol et al       |   0.75  |   1.00      | -2.9280    | 0.9909   |           |         |           |         |
Kolafa & Nezbeda |   0.75  |   1.00      |            |          | -3.3172   | 0.6951  | -3.7188   | 0.3939  |
`bd_nvt_lj`      |   0.75  |   1.00      | -2.93(1)   | 0.98(2)  |           |         | -3.73(1)  | 0.38(2) |
