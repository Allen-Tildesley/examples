# Examples
These programs accompany the book
"Computer Simulation of Liquids" by M. P. Allen and D. J. Tildesley,
second edition, to be published by Oxford University Press.
They are made freely available, in the hope that they will be useful.
The intention is to clarify points made in the text,
rather than to provide a piece of code
suitable for direct use in a research application.

## Warning
These programs are currently in development. They are not yet ready.

## Disclaimer
We ascribe no commercial value to the programs themselves.
Although a few complete programs are provided,
our aim has been to offer building blocks rather than black boxes.
As far as we are aware, the programs work correctly,
but we can accept no responsibility for the consequences of any
errors;
we should be grateful to hear from you if you find any.
You should always check out a routine for your particular application.

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

We hope that those who are used to other program languages
will find little difficulty in converting these examples;
also we point out the provisions,
in current Fortran standards,
for interoperability with C codes.

## Building the codes
The supplied SConstruct and SConscript files
will build all the working examples,
using SCons,
an Open Source software construction tool based on Python.
The homepage for SCons is at http://www.scons.org/.
Simply type `scons` to build each full example program in its own directory.
A few examples consist of individual routines or modules,
so there is no need to build them.

