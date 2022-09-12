# VHoppAS: Variable-range Hopping in Amorphous Solids
DOS-based kinetic Monte-Carlo simulation of 3D variable-range hopping for conductivity and self-diffusion coefficients of electrons in amorphous compounds.

## Features
- Randomized 3D positions and energies of localized states generated from an arbitrary user-supplied DOS
- Flexible cut-off definitions for DOS range, hop distance and hop energy (incl. auto-adjustment)
- Initial electron distribution: Fermi-Dirac or step function
- Multi-stage simulations with equilibration and repetitions for improved statistics
- Electron transport with and without electric field
- Calculation of Haven ratio and partial entropy of electrons
- Facile parameter studies, e.g. of temperature, chemical potential or electric field (optimized for parallel submission to cluster systems)
- Optional fit of the electron distribution for effective temperature and effective chemical potential (e.g. at high electric field strength)
- Detection of oscillations
- Extensive statistics (incl. detailed histograms)
- _and more_

## Documentation
- [User manual](USER_MANUAL.md): installation instructions and usage guide with examples
- [Technical notes](TECHNICAL_NOTES.md): theoretical background of the simulations and implementation details
- [Developer manual](DEV_MANUAL.md): guides for working on this project
- [Changelog](CHANGELOG.md)

## License information
This code is &copy; P. Hein, 2022, and it is made available under the GPL license enclosed with the software.

Over and above the legal restrictions imposed by this license, if you use this software or modified variants of it for an academic publication then you are obliged to provide proper attribution. This can be to this code directly,

- P. Hein. VHoppAS: Variable-range Hopping in Amorphous Solids, _\<version\>_ (2022). github.com/phi-hein/VHoppAS.

or to the paper that describes it

- _to be published_.

or (ideally) both.
