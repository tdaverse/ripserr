## Submission comments for version 1.0.0 (2025 June ??)

This is the first major release increment, due to complete replacement of the outdated Ripser C++ library by its most recent version.
To limit the number of CRAN submissions, v0.4.0 (black hole data sets) is subsumed into this release.

### Test environments

The submission version is being tested using the following tools:

* local OS X install, R v4.2.3
  * `devtools::check()`
  * `devtools::check(env_vars = c('_R_CHECK_SUGGESTS_ONLY_' = "true"))`
* local OS X install, R v4.4.2
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

#### local OS X installs, R v4.2.3 and v4.4.2

There were no ERRORs or WARNINGs. An occasional NOTE was presumably due to Internet connection speeds.

#### Win-Builder

The old release check produced the following extra component to the "incoming feasibility" NOTE:

```
Found the following (possibly) invalid URLs:
  URL: http://www.ipeadata.gov.br/Default.aspx
    From: man/case_predictors.Rd
          inst/doc/modeling-disease-using-PH.html
    Status: Error
    Message: Timeout was reached [www.ipeadata.gov.br]: Failed to connect to www.ipeadata.gov.br port 80 after 42105 ms: Could not connect to server
```

This is a different URL than was flagged in a previous submission, but it is also valid, though it triggers security flags on some browsers.
