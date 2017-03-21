# Examples
This software was written in 2016/17
by Michael P. Allen <m.p.allen@warwick.ac.uk>/<m.p.allen@bristol.ac.uk>
and Dominic J. Tildesley <dominic.tildesley@epfl.ch> ("the authors"),
to accompany the book "Computer Simulation of Liquids", second edition, 2017 ("the text"),
published by Oxford University Press ("the publishers").

## Licence
Creative Commons CC0 Public Domain Dedication.
To the extent possible under law, the authors have dedicated all copyright and related
and neighboring rights to this software to the PUBLIC domain worldwide.
This software is distributed without any warranty.
You should have received a copy of the CC0 Public Domain Dedication along with this software.
If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.

## Disclaimer
The authors and publishers make no warranties about the software, and disclaim liability
for all uses of the software, to the fullest extent permitted by applicable law.
The authors and publishers do not recommend use of this software for any purpose.
It is made freely available, solely to clarify points made in the text.
When using or citing the software, you should not imply endorsement by the authors or publishers.

## Warning
These programs are currently in development. They are __not yet ready.__

## Language
The programs contain some explanatory comments, and
are written, in the main, in Fortran 2003/2008.
This has some advantages:
a built-in syntax for array operations,
a straightforward approach to modular programming,
and a basic simplicity.
It is also a compiled language,
which means that it is quite efficient,
and widely used,
so it is easy to find compilers which are
optimised for different machine architectures.
The common tools for parallelizing scientific codes
(OpenMP and MPI)
are compatible with Fortran.
The [User Guide](./GUIDE.md) contains some notes to assist in running the program,
and some typical results.

We hope that those who are used to other program languages
will find little difficulty in converting these examples;
also we point out the provisions,
in current Fortran standards,
for interoperability with C codes.

The [python-examples subdirectory](./python_examples) contains Python versions
of several of these same examples, also with an accompanying
[User Guide](./python_examples/GUIDE.md).

## Building the codes
The supplied SConstruct and SConscript files
will build all the working examples,
using SCons,
an Open Source software construction tool based on Python.
The homepage for SCons is at <http://www.scons.org/>.
Simply type `scons` to build each full example program in its own directory.
A few examples consist of individual routines or modules,
so there is no need to build them.

The build process for the Fortran examples has been tested using SCons v2.5.1
(older versions might not work properly).
If you don't like using SCons, or can't get it to work,
it is not difficult to compile the programs using other methods.
Bear in mind that, with Fortran, it is usually essential to compile any
modules that are used by the main program, before compiling the main program itself.

The Python versions do not require building, they are simply run through the Python interpreter.
