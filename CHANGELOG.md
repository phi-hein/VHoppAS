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

## [1.9.0] - June 2023

### Changed:
- **Breaking:** Rename several output values
- **Breaking:** Add `<Data>` tag to the DOS file specification (to signal the beginning of DOS values) ([#11])
- **Breaking:** Change default value of UseYZVariance setting from no to yes

### Added:
- **Breaking:** Result files contain charge diffusion coefficient

### Fixed:
- **Breaking:** Respect that electrons have a spin that allows two electrons with opposite spin "at the same place" (note: spin-degeneracy is assumed) ([#10])
- Realistic example parameters ([#5])
- Avoid potential reading errors for values with partially equal names
- Take oscillating electrons into account for calculating the mobile ratio in state energy histograms 

## 1.8.0 - January 2023

### Added:
- **Breaking:** Fermi level-related statistics ([#4])
    - Count the number of electrons above the Fermi level and number of electron holes below the Fermi level
    - Analyze the number of hops below, across and above the Fermi level and calculate their contribution to the conductivity

## 1.7.1 - October 2022

### Changed:
- KMC optimization with regard to oscillations: Avoid looping through all electrons or paths ([#7])

### Added:
- Upper boundary for attempt time and a lower boundary for temperature

## 1.7.0 - September 2022

### Changed:
- **Breaking:** Rename the command line parameter `-sim` to `-job` (for consistency with new feature, see below)
- Allow the omission of parameter units in input files

### Added:
- Parallelize repetitions: Use the `-job` command line parameter in combination with `ParallelizeReps = y` to select a certain SimID and RepID instead of a SimID with all repetitions

### Fixed:
- Use the notation `Dsigma` (= _charge diffusion coefficient_) for the diffusion coefficient calculated from conductivity via the Nernst-Einstein equation (required for Haven ratio)

## 1.6.0 - September 2022

_Initial release._

<!--- List of links to releases: --->
[1.9.0]: https://github.com/phi-hein/VHoppAS/releases/tag/v1.9.0

<!--- List of links to pull requests and issues: --->
[#4]: https://github.com/phi-hein/VHoppAS/issues/4
[#5]: https://github.com/phi-hein/VHoppAS/issues/5
[#7]: https://github.com/phi-hein/VHoppAS/issues/7
[#10]: https://github.com/phi-hein/VHoppAS/issues/10
[#11]: https://github.com/phi-hein/VHoppAS/issues/11
