#Interface for calling CFITSIO from MATLAB.

# Introduction #

MFITSIO is the interface for calling CFITSIO procedures from MATLAB.


# Details #

Install MFITSIO:
  * Install CFITSIO to the system. Copy the run-time library `libcfitsio.a` to `/usr/lib` and header files contained in `cfitsio/include` to `/usr/local/include`.
  * Download MFITSIO from [MFITSIO Webpage](http://public.lanl.gov/eads/mfitsio).
  * Modify the `Makefile` or `.bashrc` so that `mex` can run smoothly.
  * Run `make`, and open MATLAB, then add the path where the `.m` files and `.mexglx` files are contained into MATLAB list of paths.

# Problems #

The MFITSIO project has no longer been updated since 2006. Now there are several problems:
  * There is no C compiler shipped with MATLAB for 64-bit Windows.
  * The pre-compiled binaries are DLLs which are no longer supported by MATLAB. The library compiled by `mex` and supported by MATLAB is `.mexwin32` or `.mexwin64` for Windows.