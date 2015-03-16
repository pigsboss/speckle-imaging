#Build Fortran Executables.

# Introduction #

FITSIO and BLAS libraries are required in order to build Fortran executables.


# Details #

Install the required libraries through the following steps:
  * Download CFITSIO sources from [CFITSIO Homepage](http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html).
  * Compile and link the CFITSIO sources, then copy the run-time library `libcfitsio.a` to the default path (for example, `/usr/lib`).
  * Download BLAS sources from [Netlib Repository](http://www.netlib.org/blas/blas.tgz).
  * Modify the `make.inc` file.
  * Run `make`, then copy the run-time library `libblas.a` to the default path too.