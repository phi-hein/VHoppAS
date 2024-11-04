# VHoppAS: Variable-range Hopping in Amorphous Solids
Kinetic Monte-Carlo simulation of 3D variable-range hopping of localized electrons in amorphous compounds, which yields the conductivity (and diffusion coefficients) based on a user-supplied density of states (e.g. from an ab-initio DOS).

This code was developed by Philipp Hein at the Institute of Physical Chemistry, RWTH Aachen University, Germany, under the supervision of Prof. Manfred Martin.

DOI (scientific article): [10.1002/aelm.202400407](https://doi.org/10.1002/aelm.202400407)  
DOI (citable archive): [10.5281/zenodo.10974845](https://doi.org/10.5281/zenodo.10974845)

## Features
- Randomized 3D positions and energies of localized states generated from an arbitrary user-supplied DOS
- Flexible cut-off definitions for DOS energy range, hop distance and hop energy (incl. auto-adjustment)
- Initial electron distribution: Fermi-Dirac or step function
- Multi-stage simulations with equilibration and repetitions for improved statistics
- Electron transport with and without electric field
- Validation of the Nernst-Einstein relation
- Calculation of the partial entropy of electrons
- Facile parameter studies, e.g. of temperature, chemical potential or electric field (optimized for parallel submission to cluster systems)
- Optional fit of the electron distribution for effective temperature and effective chemical potential (e.g. at high electric field strength)
- Detection of oscillations
- Extensive statistics (incl. detailed histograms)
- _and more_

## Documentation
- [User manual](USER_MANUAL.md): Installation instructions and usage guide with examples
- [Technical notes](TECHNICAL_NOTES.md): Theoretical background of the simulations and implementation details
- [Developer manual](DEV_MANUAL.md): Guides for working on this project
- [Changelog](CHANGELOG.md)

## License information
This code is &copy; P. Hein, 2023, and it is made available under the GPL license enclosed with the software.

Over and above the legal restrictions imposed by this license, if you use this software or modified variants of it, for example for an academic publication, then you are obliged to provide proper attribution by citing the following paper:

"Variable-Range Hopping Conduction in Amorphous, Non-Stoichiometric Gallium Oxide"  
P. Hein, T. Romstadt, F. Draber, J. Ryu, T. Böger, A. Falkenstein, M. Kim and M. Martin  
_Advanced Electronic Materials_ 2024, 2400407  
DOI: [10.1002/aelm.202400407](https://doi.org/10.1002/aelm.202400407)  
