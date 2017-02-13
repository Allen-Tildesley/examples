# Brief Guide to Python Examples
This subdirectory contains Python versions of some of the example programs.
Python has some advantages over Fortran:
it is an interpreted language, rather than a compiled one,
which allows faster code development,
and it has a natural (and readable) programming style.
Variables may be created and redefined at runtime, as required.
Object-oriented programming fits very well with the Python design,
but it is well suited to other styles as well.
It is widely used as a vehicle to introduce students to scientific programming.
For an excellent introductory text, see

* _Learning Scientific Programming with Python,_ C Hill (Cambridge University Press, 2015).

Python has one major drawback, when compared with Fortran and other compiled languages:
it is very slow in execution.
To partially counter this,
the [NumPy](http://www.numpy.org/ "NumPy home page") package offers more efficient handling of array structures.
The syntax of NumPy is surprisingly close to that of Fortran in some respects,
especially for combining arrays.
On top of this, NumPy and the [SciPy](https://www.scipy.org/ "SciPy home page") package
provide a huge variety of scientific libraries,
and our Python examples make extensive use of these.
Nonetheless,
any code which cannot be handled in vectorized form by NumPy will still run slowly.
Several strategies are available to address these issues
(for example, [Cython](http://cython.org/ "Cython home page"),
[SWIG](http://www.swig.org/ "SWIG home page"),
and F2PY which is part of NumPy),
but we do not attempt to follow these here.
As always,
the main aim of these examples is to illustrate ideas in the text,
not to provide programs for practical use.
The reader should feel free to experiment with ways to make the programs run faster!

In the past few years, the community has been making the transition from Python 2 to Python 3.
There are some incompatibilities between the two,
and since a choice had to be made,
we have settled on __Python 3__ for these examples.
We indicate this by the string
```
#!/usr/bin/env python3
```
at the top of each source file.
On many systems,
this will allow the files to be run directly by typing their name;
otherwise it will be necessary to type, for example, `python3 t_tensor.py`,
or just `python t_tensor.py`,
depending on your particular installation of Python.
The examples _will not work_ with Python 2!
They have been tested with Python 3.6.0.
For an introduction to the differences between Python 2 and Python 3,
see the [What's New in Python 3.0](https://docs.python.org/3/whatsnew/3.0.html "What's New in Python 3.0") page.
The most obvious changes are

1. `print(a)` is a function; `print a` will return an error.
2. An expression like `1/2` will return a float; if you want a truncated integer, use `1//2`.
3. Various methods, such as `dict.keys()`, return views,
and various functions, such as `zip`, return iterators, instead of lists.

Anyone coming from a Fortran background should note the use of indentation in Python
to indicate the range of conditional constructs and loops;
there is no `end if` statement!
Fortran experts
should also be aware that indices for arrays
(and other entities) follow the C-convention of numbering from 0 upwards.
Arrays cannot have negative indices; rather, the notation `a[-1]` refers to
numbering backwards from the end of the array.
A more subtle point is that, in the slice notation `a[i:j]`,
the upper index `j` is _excluded_, so for example `a[1:3]` consists of the elements
`a[1]` and `a[2]`,
whereas in Fortran the analogous notation `a(i:j)` refers to elements `i` through `j` _inclusive_.
Yet more subtlety lies in the distinction between
an assignment statement which makes a fresh copy of an object,
such as an array or array slice,
and an assignment which merely gives a name to a view of that object.
In the second case,
no new memory locations are used, and
changing the value associated with one of the objects will affect both of them:
they are really the same object.
This is a powerful, and memory-efficient, approach, but can be confusing!

## Data Input
In the Fortran examples we use a `namelist` to input a few parameters from standard input,
but Python does not have this.
Instead,
to provide a keyword based syntax,
we input these values using the widespread [JSON](http://www.json.org/ "JSON home page") format.
Typical input for a very simple example might be
```
{ "nblock":20, "nstep":10000, "dt":0.002 }
```
and the `"key":value` pairs may be set out on different lines if you wish.
The appearance is very similar to a Python dictionary,
and indeed the data is loaded and parsed into a dictionary for further processing.
Note carefully that we use a colon `:` rather than an equals sign to separate the
key from the value,
and that the keys should be enclosed in double quotes `"..."`.
To avoid some fairly basic exceptions later in the program,
we usually type-check the data
(Python enthusiasts may disapprove of this)
so integer values must not have a decimal point,
while floating-point values must have one
(and at least one digit following the point, for example `1.0`, otherwise JSON raises an exception).
String values need to be enclosed, like the keys, in double quotes.
The variables which may be set in this way are typically considered one by one in our programs:
those whose names correspond to keys supplied in the input file are given the input values,
while the others are given values taken from another `defaults` dictionary.

## Test programs for potentials, forces and torques
Two program files are provided: `test_pot_atom.py` and `test_pot_linear.py`,
for pair potentials between, respectively, atoms and linear molecules.
These load, at runtime, a module containing a function to calculate
the necessary potential, forces and torques.
The aim is to demonstrate the numerical testing of the analytical derivatives
which go into the forces and torques:
small displacements and rotations are applied in order to do this.
The test is performed for a randomly selected configuration.
Some parameters are used to prevent serious overlap,
which might produce numerical overflow,
while keeping the particles close enough together to give non-zero results.
The values of these parameters may be adjusted via the input file in individual cases.
To run the programs without any tweaking,
simply give (through standard input in the usual way)
a record containing the bare minimum information,
namely a string which identifies
the model of interest, for example `{"model":"at"}`.
The supplied examples (with their identifying strings) are, for `test_pot_atom.py`:

* `"at"` `test_pot_at.py` the Axilrod-Teller three-body potential
* `"bend"` `test_pot_bend.py` the angle-bending part of a polymer chain potential
* `"twist"` `test_pot_twist.py` the angle-torsion part of a polymer chain potential

and for `test_pot_linear.py`

* `"dd"` `test_pot_dd.py` the dipole-dipole potential
* `"dq"` `test_pot_dq.py` the dipole-quadrupole and quadrupole-dipole potential
* `"gb"` `test_pot_gb.py` the Gay-Berne potential
* `"qq"` `test_pot_qq.py` the quadrupole-quadrupole potential

## T-tensor program
The program `t_tensor.py` compares the calculation of multipole energies by two methods:
using explicit formulae based on trigonometric functions of the Euler angles,
and via the Cartesian T-tensors.
Two linear molecules are placed in random positions and orientations,
within a specified range of separations,
and some of the contributions to the electrostatic energies and forces are calculated.
The program may be run using an empty record `{}`,
so as to take the program defaults,
or various parameters may be specified.
Several of the tensor manipulations are neatly expressed using NumPy library functions
such as `outer` (outer product) and `einsum` (Einstein summation).

* How easy would it be to add quadrupole-quadrupole energy, quadrupole-dipole forces,
quadrupole-quadrupole forces, and all the torques, calculated both ways??

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

In the current implementation,
the direct method makes little or no use of NumPy's efficient array manipulation functions,
and so is very slow.
For this reason,
the default run length is much shorter than for the Fortran example.
An alternative, much faster, direct method is also provided at the end of the program.
This uses the NumPy `correlate` library function.
The selected mode `'full'` corresponds to aligning the data array `v` with itself,
at all possible offsets that result in some overlap of values,
before computing the products (for each offset) and summing.
In this case, the resulting array is symmetric in time (offset),
and only the values from the mid-point onwards are required;
these must be normalized in the standard way
and identical results to the slow direct method (with `origin_interval=1`)
and FFT method,
are obtained.

## Pair distribution function
The program `pair_distribution.py` reads in a set of configurations and calculates
the pair correlation function _g(r)_.
We limit the number of configurations to a maximum of 1000 (numbered from 000 to 999)
simply so as to use a fixed naming scheme for the input configurations;
in a practical application, a trajectory file would be used instead.
The sum over all pairs is performed using a vectorized approach,
so as to take advantage of NumPy routines.
The coordinate array `r` is compared with a cyclically-shifted copy,
to give _N_ separation vectors `rij` which are processed all at once.
Only `n//2` shifts are needed to cover every distinct `ij` pair.
There is a slight subtlety on the last shift, if `n` is even:
both `ij` and `ji` pairs appear,
and so the usual incrementing factor 2 is replaced by a factor 1.
The actual histogramming is conveniently performed
by the built-in NumPy `histogram` routine.

We have tested this on a set of 500 configurations
of _N_=256 Lennard-Jones atoms,
cut (but not shifted) at _R_<sub>c</sub>=2.5&sigma;,
at the usual state point &rho;=0.75, _T_=1.0.
The interval between configurations was 100 MC sweeps.
Using the default resolution of 0.02&sigma;,
identical results were obtained as for the Fortran example.

![alt text](gr.png "g(r) test results")

## Error calculation
The program `error_calc.py` is a self-contained illustration of the effects of
correlations on the estimation of errors for a time series.
We produce the series using a generalized Langevin equation,
in the same manner as for the correlation function program (see above).
Since the correlation time of the GLE is exactly known,
we can predict the effects, and compare with the empirical estimates
obtained by different methods.
The program contains extensive comments to explain what is being calculated at each stage.

## FFT program
The aim of `fft3dwrap.py` is to illustrate the way a standard Fast Fourier Transform
library routine is wrapped in a user program.
We numerically transform a 3D Gaussian function,
and compare with the analytically, exactly, known result,
User input defines the number of grid points and the box size;
sensible defaults are provided.
The library that we use for this example is the built-in NumPy one.

## Hit-and-miss and sample-mean
The two programs `hit_and_miss.py` and `sample_mean.py` illustrate two very simple
Monte Carlo methods to estimate the volume of a 3D object.
They are both described in detail at the start of Chapter 4.
No user input is required.
For the built-in values defining the geometry, the exact result is 5/3.

## Quantum simulation programs
The program `qmc_walk_sho.py` solves the diffusion equation in imaginary time
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
with an empty input record `{}`;
these defaults include setting `et` initially to the exact ground state energy.
Other values such as `{et:0.6}` may be supplied through standard input in the usual way.
This type of simulation is sensitive to the initial value,
and quite noisy:
possible improvements are discussed in general terms in the text.
