# next version

## output class

The new `return_class` parameter allows the user to specify whether to output persistence data in the legacy `'PHom'` class or the `'persistence'` class from {phutil}.
It defaults to `'PHom'` but will switch to `'persistence'` when the `'PHom'` class is deprecated in a future version.

# ripserr 1.0.0

This major version replaces an outdated version of the Ripser C++ library with its current version.

# ripserr 0.4.0

This minor edition includes the `blackhole` data set along with its accompanying documentation and R script.

# ripserr 0.3.0

This minor version subsumes several contributions over 4+ years, including invisibly returning persistence data on `print()`, an additional vignette, and debugging & correcting of documentation.

# ripserr 0.2.0

## Changes

* `vietoris_rips` and `cubical` are now S3 generics
* `vietoris_rips` accepts objects of class `data.frame`, `matrix`, `numeric`, `ts`, and `dist`
* `cubical` accepts objects of class `array` and `matrix` (for r-oldrel, in which not all `matrix` objects are also of class `array`)
* `vietoris_rips` and `cubical` return objects of class `PHom` and `data.frame` (instead of either `data.frame` or `matrix` previously)

# ripserr 0.1.1

## Logistics

* Added Code of Conduct (CODE_OF_CONDUCT.md) from Contributor Covenant
* Added Contribution Guide (CONTRIBUTING.md) from ggalluvial
* Removed empty vignette (`cubical`), will add in next minor update (likely 0.2.x)

## Bugs

* Fixed memory leak in Cubical Ripser C++ code
* Fixed citation issue in vignette (`vietoris-rips`)
* Removed TDAstats from Suggests to avoid circular dependency in the future
* Fixed issue in old R versions where 2-d matrices are not automatically also of type "array" (cubical() checks)

# ripserr 0.1.0

* Initial release
* Ports Ripser for Vietoris-Rips complex functionality
* Ports Cubical Ripser for cubical complex functionality
