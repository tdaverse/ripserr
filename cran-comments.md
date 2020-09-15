This is an patch to ripserr due to ERRORs in tests on CRAN. The deadline to submit this update and avoid removal of ripserr from CRAN was 28 Sep 2020 (today is 15 Sep 2020).

The ERROR was caught by AddressSanitizer and was a memory leak in the underlying C++ code on which ripserr relies. The memory leak has been fixed and checked by using valgrind.

There was also an ERROR on a CRAN check (mac_oldrel) because ripserr assumed objects of class "matrix" are also of class "array"; this issue has been fixed and the error no longer occurs.

## Test environments

* local OS X install, R 4.0.2
* Ubuntu 16.04 (on Travis CI), R 4.0.0
* Windows Server 2012 R2 (on AppVeyor CI), R 4.0.2
* win-builder, r-devel

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs.