# Examples
This software was written in 2016/17
by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>
and Dominic J. Tildesley <d.tildesley7@gmail.com> ("the authors"),
to accompany the book
[Computer Simulation of Liquids](https://global.oup.com/academic/product/computer-simulation-of-liquids-9780198803201),
second edition, 2017 ("the text"),
published by Oxford University Press ("the publishers").

## Licence
Creative Commons CC0 Public Domain Dedication.
To the extent possible under law, the authors have dedicated all copyright and related
and neighboring rights to this software to the PUBLIC domain worldwide.
This software is distributed without any warranty.
You should have received a copy of the
[CC0 Public Domain Dedication](./COPYING.txt)
along with this software.
If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.

## Disclaimer
The authors and publishers make no warranties about the software, and disclaim liability
for all uses of the software, to the fullest extent permitted by applicable law.
The authors and publishers do not recommend use of this software for any purpose.
It is made freely available, solely to clarify points made in the text.
When using or citing the software, you should not imply endorsement by the authors or publishers.

## Language
The programs contain some explanatory comments,
and are written, in the main, in Fortran 2003/2008.
This has some advantages: a built-in syntax for array operations,
a straightforward approach to modular programming, and a basic simplicity.
It is also a compiled language, which means that it is quite efficient,
and widely used, so it is easy to find compilers which are
optimized for different machine architectures.
The common tools for parallelizing scientific codes (OpenMP and MPI)
are compatible with Fortran.
The [User Guide](./GUIDE.md) contains some notes to assist in running the programs,
and some typical results.

We hope that those who are used to other program languages
will find little difficulty in converting these examples;
also we point out the provisions, in current Fortran standards,
for interoperability with C codes.

The [python-examples subdirectory](./python_examples) contains Python versions
of several of these same examples, also with an accompanying
[User Guide](./python_examples/GUIDE.md).

## Building the codes
On some computing platforms,
the supplied `SConstruct` and `SConscript` files will build all the working examples,
using SCons, an Open Source software construction tool based on Python.
The homepage for SCons is at <http://www.scons.org/>.
The `SConstruct` file may need to be edited
(for example, to point to the correct location of libraries
such as `fftw3` and `lapack` on your system).
Then, simply type `scons` to build each full example program in its own directory.
A few examples consist of individual routines or modules,
rather than working programs,
so there is no need to build them.

The build process for the Fortran examples has been tested using SCons v2.5.1
(older versions might not work properly).
If you don't like using SCons, or can't get it to work,
it is not difficult to compile the programs using other methods.
Bear in mind that, with Fortran, it is usually essential to compile any
modules that are used by the main program, before compiling the main program itself.
Take a look at the `SConstruct` file in any case,
as it shows the file dependencies for each example.
Also bear in mind that several alternative module files
(e.g. `md_lj_module.f90`, `md_lj_ll_module.f90`, `md_lj_omp_module.f90`)
contain modules with the same name (`md_module` in this case),
since they act as drop-in replacements for each other.
To avoid confusion during compilation due to intermediate files with the same name
(e.g. `md_module.mod`)
it is advisable to __compile each example in its own build directory__
(which is what the `SConstruct` file is configured to do)
or to delete all intermediate files before each individual compilation.

We have used gfortran v7.2 (and before that, v6.3) for testing,
but have attempted to stick to code which conforms to the Fortran 2008 standard.
In gfortran v6, calling the intrinsic `random_seed()` function would
generate the same sequence of random numbers every time.
In a few examples it is important to generate different sequences each time the program is run,
and for this purpose a subroutine `init_random_seed` was provided in the file `maths_module.f90`;
this routine, however, included some GNU extensions.
With gfortran v7, the random number generator has been changed,
and the intrinsic `random_seed()` function behaves in the desired fashion.
Therefore, the GNU-specific code has been transferred to a separate file, `gnu_v6_init_random_seed.f90`,
which is not included in the build process,
but may be of interest to those still using gfortran v6.
The current `init_random_seed` routine simply calls `random_seed()`.
__You should check the behaviour of the random number generator on your own system.__
Note that, by default, we do not select any optimization options in compilation.
If you are using a different compiler,
then the compiler and linker options in the `SConstruct` file will most likely need changing.

The above, general, advice should help you to build the codes on your system.
Unfortunately, due to the enormous variety of computing platforms and compilers,
__we cannot offer more specific advice on the build process.__

The Python versions do not require building, they are simply run through the Python interpreter.
