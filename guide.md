<h1>Brief Guide</h1>
Here are some notes to assist in running the programs.
Most of the Fortran codes use a **namelist** to input a few parameters from standard input.
This gives an easy way to specify default values in the program itself, and to use a 
keyword-based syntax to specify values different from the default ones at run-time.
The input file, or string, should usually begin with **&nml** and end with **/**.
As a minimum, the program will expect to read **&nml /** (an empty list), but to
change the parameters, typical input might be **&nml nblock=20, nstep=1000, dt=0.001 /**,
and the key/value pairs may be set out on different lines if you wish.

Simulation runs require a starting configuration which can usually be prepared using
the **initialize** program (in **build_initialize/**).
The default parameters produce an FCC configuration of 256 atoms at reduced density 0.75,
writing out just the positions (for an MC program) to a file **cnf.inp**.
If the parameter **velocities=.true.** is supplied, then positions and velocities are
written to the file, corresponding to a reduced temperature 1.0.

For most of the examples, we use a cutoff of 2.5 (in reduced units);
for MC programs the cut (but not shifted) potential is used in the simulation,
while for MD programs the cut-and-shifted potential is used.
In all cases, an attempt is made to produce results with and without long-range corrections.
Differences between the different models are discussed in various places,
see e.g. A Trokhymchuk, J Alejandre, J Chem Phys, 111, 8510 (1999).
Using their table V as a guide, we take the critical point to be roughly located at:

model                 | temperature | density | pressure
-----                   -------   -----------   --------
full, with LRC        | 1.31 | 0.31 | 0.13
cut (but not shifted) | 1.19 | 0.32 | 0.11
cut-and-shifted       ! 1.08 | 0.32 | 0.09

At any temperature below Tc, the liquid state is bounded below by the
liquid-gas coexistence density, and using Tables II-IV of the same reference as a guide,
we take the values for three example temperatures as

model                 | temperature | density | pressure
-----                   -------   -----------   --------
full, with LRC        | 0.8 | 0.793 | 0.005
cut (but not shifted) | 0.8 | 0.765 | 0.008
cut-and-shifted       ! 0.8 | 0.730 | 0.013
full, with LRC        | 0.9 | 0.746 | 0.012
cut (but not shifted) | 0.9 | 0.714 | 0.020
cut-and-shifted       ! 0.9 | 0.665 | 0.030
full, with LRC        | 1.0 | 0.695 | 0.026
cut (but not shifted) | 1.0 | 0.652 | 0.036
cut-and-shifted       ! 1.0 | 0.578 | 0.062

For Lennard-Jones, the default state point (density=0.75 and temperature=1.0) 
lies in the liquid region of the phase diagram for all these variants of the model. 

For comparison with MD runs, a wide range of simulation data for LJ with 
rcut=2.5 (cut-and-shifted, denoted cs below) is summarized by
M Thol, G Rutkai, R Span, J Vrabec and R Lustig, Int J Thermophys, 36, 25 (2015)
who present data (in Supplementary Information) and a fitted approximate 
equation of state (which they also provide as a C++ program). 
Full LJ potential data are available from
J Kolafa and I Nezbeda, Fluid Phase Equil, 100, 1 (1994) and
M Mecke, A Muller, J Winkelmann, J Vrabec, J Fischer, R Span, W Wagner,
Int J Thermophys, 17, 391 (1996) and 19, 1493(E) (1998)
Below we use the Kolafa & Nezbeda formulae (denoted f) and estimate the results
for the cut-but-not-shifted potential (denoted c) using the same corrections as in the MC codes.
Below we compare with typical test runs from our programs using default parameters except where stated.
Note that E/N is the total internal energy per atom, including the ideal gas part.

Source           | density | temperature | E/N (cs) | P (cs)  | E/N (c) | P (c)  | E/N (f)  | P (f)   |
------           | ------- | ----------- | -------- | ------  | ------- | ------ | -------- | ------- |
Thol et al       |   0.75  |   1.00      | -2.9280  | 0.9909  |         |        |          |         |
Kolafa & Nezbeda |   0.75      1.00      |          |         | -3.3172 | 0.6951 | -3.7188  | 0.3939  |
-------          | ------- | ----------- | ---------| ------  | ------- | ------ | -------- | ------- |
bd nvt lj        |   0.75  |   1.00      | -2.93(1) | 0.98(2) |         |        | -3.73(1) | 0.38(2) |
