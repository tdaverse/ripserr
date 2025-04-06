## Submission comments for version 0.3.0 (2025 April 6)

This is a minor version in anticipation of an upcoming major version. The goal is to incorporate, stabilize, and make available a variety of contributions made over several years before introducing fundamental changes that will need to be reconciled with the rest of the package.
Comments from previous submissions are retained below for reference.

### Test environments

The submission version is being tested using the following tools:

* Local OS X installs, R v4.2.3 and v4.4.2
  * `devtools::check()`
  * `devtools::check(env_vars = c('_R_CHECK_SUGGESTS_ONLY_' = "true"))`
  * `devtools::check(env_vars = c('_R_CHECK_DEPENDS_ONLY_' = "true"), vignettes = FALSE)`
  * `devtools::check(remote = TRUE, manual = TRUE)`
* Win-Builder
  * `devtools::check_win_oldrelease()`
  * `devtools::check_win_release()`
  * `devtools::check_win_devel()`

The results are as follows.
Alerts (ERRORs, WARNINGs, or NOTEs) that appeared in previous checks are not mentioned in the results of subsequent checks.
Only tests that produced alerts not produced by previous alerts are reported and discussed.

#### Local OS X install, R v4.2.3

Some checks produced no alerts.

The `R_CHECK_DEPENDS_ONLY_` check produced the following NOTE:

```
── R CMD check results ────────────────────────────────────── ripserr 0.3.0 ────
Duration: 36.2s

❯ checking for future file timestamps ... NOTE
  unable to verify current time

0 errors ✔ | 0 warnings ✔ | 1 note ✖
```

This was likely due to a weak Internet connection.

The `remote = TRUE` check was not completed on this machine due to GitHub authentication problems.

#### Local OS X install, R v4.4.2

All checks produced no ERRORs or WARNINGs.

The last check, `devtools::check(remote = TRUE, manual = TRUE)`, produced the following NOTE, reflecting the transfer of maintenance:

```
❯ checking CRAN incoming feasibility ... [4s/18s] NOTE
  Maintainer: ‘Jason Cory Brunson <cornelioid@gmail.com>’
  
  New maintainer:
    Jason Cory Brunson <cornelioid@gmail.com>
  Old maintainer(s):
    Raoul Wadhwa <raoulwadhwa@gmail.com>
```

It also produced many notes regarding the HTML manual, beginning with the following:

```
❯ checking HTML version of manual ... NOTE
  Found the following HTML validation problems:
  PHom.html:15:44 (PHom.Rd:5): Error: <main> is not recognized!
  PHom.html:15:44 (PHom.Rd:5): Warning: discarding unexpected <main>
```

I don't know how to resolve these.
They do not arise from the Win-Builder checks and are not flagged for the previous CRAN version, so i presume that they do not pose a problem for CRAN.

#### Win-Builder

Some Win-Builder checks produced no additional alerts.

The current release check produced the following extra component to the "incoming feasibility" NOTE:

```
Found the following (possibly) invalid URLs:
  URL: https://www.ibge.gov.br/geociencias/organizacao-do-territorio/estrutura-territorial/15761-areas-dos-municipios.html?edicao=30133&t=acesso-ao-produto
    From: inst/doc/modeling-disease-using-PH.html
    Status: Error
    Message: Timeout was reached [www.ibge.gov.br]: Operation timed out after 60008 milliseconds with 0 bytes received
```

This URL is valid; it is in fact an update to a previous URL that became invalid since we acquired the data from it. See the discussion i initiated here: <https://stat.ethz.ch/pipermail/r-package-devel/2025q2/011578.html>

## Submission comments for archived package

### From fourth resubmission (15 Oct 2020)

The lines of C++ code causing the AddressSanitizer ERROR were modified so that each use of new is paired with an appropriate use of delete to prevent memory leaks.

### From third resubmission (1 Oct 2020)

The lines of C++ code causing the AddressSanitizer ERROR were modified so that random access (accessing memory that was not formally allocated, only reserved) is not used. Instead, the std::vector STL container is cleared and then new elements are added using push_back, which should avoid raising a container-overflow error with AddressSanitizer.

## From previous submissions to correct archived version

### From second resubmission (22 Sep 2020)

This is the second resubmission for a patch to ripserr due to ERRORs in tests on CRAN. The deadline to submit this update and avoid removal of ripserr from CRAN was 28 Sep 2020 (first resubmission 15 Sep 2020; second resubmission today 22 Sep 2020).

The ERROR that persisted in AddressSanitizer occurred because the C++ std::vector reserve method was used to ensure capacity; objects were then assigned to the vector without actually initializing the reserved capacity. Thus, AddressSanitizer caught a container overflow error. This has been fixed by replacing the use of the std::vector reserve method with std::vector resize, which allocates + initializes the vector's elements. This should resolve the issue brought up by AddressSanitizer.

The robustness of ripserr's underlying C++ code has been improved by CRAN's rigorous checks. Thank you for your patience with these resubmissions and your help + guidance in finding the issue.

### From first resubmission (15 Sep 2020)

The ERROR was a memory leak in the underlying C++ code on which ripserr relies. The memory leak has been fixed and checked by using valgrind.

There was also an ERROR on a CRAN check (mac_oldrel) because ripserr assumed objects of class "matrix" are also of class "array"; this issue has been fixed and the error no longer occurs.

## Test environments

The following environments were rechecked prior to this resubmission.

* local OS X install, R 4.0.2
* Ubuntu 16.04 (on Travis CI), R 4.0.0
* Windows Server 2012 R2 (on AppVeyor CI), R 4.0.2
* win-builder, r-devel

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE: submission for a package that was archived.
