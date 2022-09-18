# Changelog

<!--- 
CHANGELOG STYLE GUIDE

Use the following change categories for each release (in that order):
###Changed : for changes in existing functionality
###Added   : for new functionality
###Removed : for removed functionality
###Fixed   : for bug fixes

Describe changes in imperative, e.g. "modify exception handling for ..."

Prefix breaking changes with **Breaking:** (and list them before other changes)

Syntax for links to issues (or pull requests):
([#1](https://github.com/phi-hein/VHoppAS/issues/1))
--->

## [1.7.0] - 19.09.2022

### Changed:
- **Breaking:** Rename the command line parameter "-sim" to "-job" (for consistency with new feature, see below)
- Allow the omission of parameter units in input files

### Added:
- Parallelize repetitions: Use the "-job" command line parameter in combination with "ParallelizeReps = y" to select a certain SimID and RepID instead of a SimID with all repetitions

### Fixed:
- Use the notation "Dsigma" (= "_conductivity diffusion coefficient_") for the diffusion coefficient calculated from conductivity via the Nernst-Einstein equation (required for Haven ratio)

## [1.6.0] - 12.09.2022

_Initial release._

<!--- List of links to releases: --->
[1.7.0]: https://github.com/phi-hein/VHoppAS/releases/tag/v1.7.0
[1.6.0]: https://github.com/phi-hein/VHoppAS/releases/tag/v1.6.0
