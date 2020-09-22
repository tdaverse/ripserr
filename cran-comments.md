This is the second resubmission for a patch to ripserr due to ERRORs in tests on CRAN. The deadline to submit this update and avoid removal of ripserr from CRAN was 28 Sep 2020 (first resubmission 15 Sep 2020; second resubmission today 22 Sep 2020).

The ERROR that persisted in AddressSanitizer occurred because the C++ std::vector reserve method was used to ensure capacity; objects were then assigned to the vector without actually initializing the reserved capacity. Thus, AddressSanitizer caught a container overflow error. This has been fixed by replacing the use of the std::vector reserve method with std::vector resize, which allocates + initializes the vector's elements. This should resolve the issue brought up by AddressSanitizer.

The robustness of ripserr's underlying C++ code has been improved by CRAN's rigorous checks. Thank you for your patience with these resubmissions and your help + guidance in finding the issue.

## From first resubmission

The ERROR was a memory leak in the underlying C++ code on which ripserr relies. The memory leak has been fixed and checked by using valgrind.

There was also an ERROR on a CRAN check (mac_oldrel) because ripserr assumed objects of class "matrix" are also of class "array"; this issue has been fixed and the error no longer occurs.

## Test environments

The following environments were rechecked prior to this resubmission.

* local OS X install, R 4.0.2
* Ubuntu 16.04 (on Travis CI), R 4.0.0
* Windows Server 2012 R2 (on AppVeyor CI), R 4.0.2
* win-builder, r-devel

## R CMD check results

There were no ERRORs, WARNINGs, or NOTEs.
