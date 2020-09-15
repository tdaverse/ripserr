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