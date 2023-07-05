# User manual
This document describes how the VHoppAS program can be build on the target system (where the simulations should run) and how to set up and execute simulations.

## How to build the VHoppAS program from source
Ensure that the following two components are installed on the target system (on cluster systems this may involve loading of the correct modules):

- **C++ compiler** that supports at least the **C++17** standard (e.g. GCC with version 8 or above; can be tested with `g++ --version`) 

- **CMake** (version 3.10 or above; can be tested with `cmake --version`) that is linked to the C++ compiler (e.g. by default or through the `CXX` environment variable)

&rarr; installation via package manager: `sudo apt update && sudo apt install g++ cmake`  
(see also: [gcc install](https://gcc.gnu.org/install/), [gcc binaries](https://gcc.gnu.org/install/binaries.html) or [cmake install](https://cmake.org/install/))

Download the source code of VHoppAS from [Github](https://github.com/phi-hein/VHoppAS) to a folder on the target system (unpack if necessary) and execute the following commands in this folder (which contains the `CMakeLists.txt` file):
- `mkdir build_Release`
- `cd build_Release`
- `cmake -DCMAKE_BUILD_TYPE=Release ..`
- `cmake --build .`

## Input files and simulation parameters
Each simulation requires two text files as input: 
1. a **DOS file** that contains the density of states of the investigated material
2. an **input file** that contains the simulation parameters

The structure and content of these files will be described below. It is recommended to create a parent directory for all simulations, which contains the `VHoppAS` executable and the DOS file(s). Either copy the `VHoppAS` executable from the `build_Release` folder or set up a symbolic link. For each group of simulations (i.e. for each ProjectID, see below) create a sub-folder that contains the respective input file. In this way, different simulations can use the same DOS file(s). Example directory structure:

<pre><code>- KMC:              (<i>Parent directory</i>)
   - VHoppAS        (<i>Executable</i>)
   - DOS.txt        (<i>DOS file</i>)
   - Project001:    (<i>Simulation sub-folder</i>)
      - Input.txt   (<i>Input file</i>)
   - Project002:    (<i>Another simulation sub-folder</i>)
      - Input.txt   (<i>Another input file</i>)
   ...
</code></pre> 

Example DOS and input files, that may serve as a starting point to create own files, can be generated from the VHoppAS executable with the command: `./VHoppAS -example`

### DOS file
The DOS file contains the density of states that has been obtained through ab-initio methods or an analytical model. Similar to for example DFT software, the employed DOS is a density of single-electron states, which means that each "state" can accommodate one electron. The specified DOS is considered as the density of localized states over which variable-range hopping can take place, which means that the DOS file should not contain values beyond the mobility edges. 

**DOS file layout** (e.g. `DOS.txt` file):
```
<DOS>
Type = PiecewiseLinear
RefTemp(K) = 300

<Data>
E(eV) DOS(1/cm3eV)
-1.000   3.99649E21
-0.975   3.99110E21
-0.950   4.70461E21
...
```
The DOS file must begin with the `<DOS>` tag, `Type = PiecewiseLinear` (which is the only type implemented so far) and the specification of the reference temperature `RefTemp(K) = ...` (alternatively `RefTemp = ...` but still in unit Kelvin). The meaning of `RefTemp` is explained further below. The `<Data>` tag signals the beginning of the DOS value pairs, which means that additional comments can be written between `RefTemp(K) = ...` and `<Data>`. The DOS must be specified as value pairs where the first value is the energy **relative to the charge-neutral Fermi level** (in unit eV) and the second value is the density of states (in unit 1/cm<sup>3</sup>eV). The line `E(eV) DOS(1/cm3eV)` is just explanatory for this and can be omitted. Note that the raw DOS from ab-initio methods might need to be transformed to these units and to the relative energy axis. **Important**: Zero on the `E` axis has to be at the Fermi level of the charge-neutral material, i.e. the Fermi level for the natural amount of valence electrons. It serves as reference for the actual Fermi level in the KMC-VRH simulation, which may differ to simulate a charged material with positive or negative charge density (`ChemPot` parameter in input file). `RefTemp` is the temperature, for which the Fermi level that is used as reference for the energy axis was determined (e.g. `RefTemp(K) = 0` when the 0 K-Fermi level from basic DFT-DOS calculations is used). In VHoppAS, the DOS is assumed to be temperature-independent, but the charge-neutral Fermi level (and therefore the zero position on the energy axis) can be adjusted automatically from the reference temperature to the simulation temperature (`EFTAdjust` option in input file). Closing tags for `<DOS>` and `<Data>` are not necessary at the end of the DOS file.

**Comment on spin-degeneracy**: 
VHoppAS is based on the assumption of spin-degeneracy, i.e. there is no energy difference between spin-up and spin-down states. The DOS is the sum of spin-up and spin-down DOS, which are equal. All input and output values refer to both spin-types together (even though internally only one spin-type is simulated and taken as representative for the other). Materials without spin-degeneracy can be simulated with the following workaround:  
Run a simulation with the spin-up DOS times 2 as input and another simulation with the spin-down DOS times 2 (and equal other parameters). The movement of spin-up and spin-down electrons is then simulated separately (spin-flip can be neglected in the steady-state). Using the doubled DOS in each case is necessary to have the correct cell volume. The arithmetic average (and not the sum) of the output values of both simulations then yields the desired results.

### Input file
The input file is structured into at least two sections: a block of general project settings (enclosed by `<MC-Project>...</MC-Project>`) and a block of simulation parameters (enclosed by `<Parameters>...</Parameters>`). Simulations are organized as projects in order to group related simulations, thereby allowing to successively add more simulations or to generate result tables afterwards (see `./VHoppAS -collect` command below). A project could be for example a temperature series, i.e. different simulations with varied temperature. Such multiple simulations could be specified as different single-type input files, but it is easier to just use one multi-type input file containing the optional third section (beginning with `<VariedParameters>` and without closing tag) that defines different parameter sets as a table. This table may contain any of the simulation parameters as columns (separated by one or more spaces, not tabs) and each line of the table (together with the non-varied parameters in the `<Parameters>`-section) corresponds to a different simulation. 

**SimID**: When the input file is processed, an identification number called `SimID` gets assigned to each line in the `<VariedParameters>` table, starting with `SimID = 1` for the first line, `SimID = 2` for the second line and so on, in order to identify the different parameter sets. In single-type input files without a `<VariedParameters>` section the single parameter set gets `SimID = 1`. The SimIDs are important for simulating only a selected parameter set (`-job` command line parameter, see below) and for matching result files to the specified parameter sets. Although not recommended, it is possible to explicitly specify the SimID (>= 1) in the `<Parameters>` block of single-type input files (`SimID = ...`). In multi-type input files (= those with a `<VariedParameters>` table) any exclicit `SimID` specification is ignored.

**RepID**: Another identification number is the `RepID`. Each parameter set (= each SimID) can be simulated multiple times through the `Repetitions` parameter (see further below). The `RepID` is the identification number for individual repetitions. It is 1-based like the SimID, such that `RepID = 1` is the first repetition (= first run), `RepID = 2` is the second and so on. The `RepID` is not present in the input file, but it is written to the result files (together with the SimID) to identify and discern the results of different repetitions.

**Input file layout** (e.g. `Input.txt` file; descriptions of settings and parameters are below; see also [technical notes](TECHNICAL_NOTES.md) for the meaning of quantities):
```
<MC-Project>
ProjectID = 1
Name = Test_TempFieldVar
DOS-File = "../DOS.txt"   
Output-File = "Result.txt"
Verbosity(0-2) = 2
EFTAdjust(y/n) = y
InitialFDDistrib(y/n) = y
TeffFit(y/n) = n
EnforceECount(y/n) = y
CutoffAutoAdjust(y/n) = y
DistAdjustPercentage = 10
EdiffAdjustPercentage = 25
OnlyCompareSimID(y/n) = n
UseYZVariance(y/n) = n
ParallelizeReps(y/n) = n
Description = Testing the simultaneous variation of temperature and electric field strength. 
</MC-Project>

<Parameters>
Repetitions = 3
Emin(eV) = -0.5
Emax(eV) = 0.5
ChemPot(eV) = 0.0
States = 120000
MinPaths = 60
DistCutoff(nm) = 2.14
EdiffCutoff(eV) = 0.75
PreHops = 20000000
EqHops = 300000000
HopLimit = 300000000
Seed = -61522861958
AttTime(s) = 3.25e-16
LocRadius(nm) = 0.15
</Parameters>

<VariedParameters>
Temp    GradPhi
K       V/cm 
273.15  0.0
320     60000
320     80000
...
```

The settings, parameters and table columns can be specified in any order (within each section), except for `Description` that (if present) has to be the last in the `<MC-Project>`-block. Units (incl. the brackets) and the unit line of the table may be omitted, but it is not possible to use other units. The content of the input file is case-sensitive.

**List of input settings**:  
The following settings belong into the `<MC-Project>` block. See also [technical notes](TECHNICAL_NOTES.md) for further clarifications.
- `ProjectID` (mandatory):  
_Type: Unsigned integer; Unit: None; Allowed values: Unrestricted._  
Identification number of the project, e.g. to discern output files that belong to this input file from those who do not belong to it.

- `Name` (optional, default = "Default"):  
_Type: Text; Allowed characters: Unrestricted (no line break)._  
Name of the project.

- `DOS-File` (mandatory):  
_Type: Path; Allowed values: Relative or absolute path to existing file (no line break)._  
Path to the DOS file that should be used. Must be enclosed by double-quotes (") when the path contains white-space characters. If relative path, then relative to working directory.

- `Output-File` (optional, default = "DefaultOutput.txt"):  
_Type: Path; Allowed values: Valid relative or absolute path (directory must exist; no line break)._  
Base name and path of the output file(s). Must be enclosed by double-quotes (") when the path contains white-space characters. If relative path, then relative to working directory.

- `Verbosity` (optional, default = `1`):  
_Allowed values: `0`, `1` or `2`._  
Determines how much information is reported during and after the simulation: `0` = Minimum, `1` = Medium and `2` = Maximum. Has no effect on result file(s), except that when `Verbosity = 0` no convergence tables are recorded (reduces memory requirements slightly).

- `EFTAdjust` (optional, default = `y`):  
_Type: Switch; Allowed values: `yes` (`y`|`yes`|`t`|`true`) or `no` (`n`|`no`|`f`|`false`)._  
Determines whether the reference Fermi level (zero on energy axis) will be automatically adjusted to the simulation temperature (= `yes`) or not (= `no`).

- `InitialFDDistrib` (optional, default = `y`):  
_Type: Switch; Allowed values: `yes` (`y`|`yes`|`t`|`true`) or `no` (`n`|`no`|`f`|`false`)._  
Determines whether at the beginning of the simulations the electrons will be distributed according to the Fermi-Dirac distribution (= `yes`) or according to a step function from lowest to highest state energy (= `no`).

- `TeffFit` (optional, default = `y`):  
_Type: Switch; Allowed values: `yes` (`y`|`yes`|`t`|`true`) or `no` (`n`|`no`|`f`|`false`)._  
Determines whether the steady-state electron distribution will be fitted with a Fermi-Dirac function to determine the effective temperature and effective chemical potential (= `yes`), which is important to obtain correct results beyond the weak-field regime, or whether the input temperature and chemical potential are used to calculate the effective charge carrier density (= `no`).

- `EnforceECount` (optional, default = `y`):  
_Type: Switch; Allowed values: `yes` (`y`|`yes`|`t`|`true`) or `no` (`n`|`no`|`f`|`false`)._  
Determines whether the total number of electrons should exactly correspond to the expectation from the Fermi-Dirac distribution (= `yes`) or whether stochastic deviations from the ideal electron count are allowed (= `no`). These deviations normally arise when the Fermi-Dirac probability in combination with random numbers is used to populate the localized states with electrons (i.e. when `InitialFDDistrib = y`). This option is irrelevant when `InitialFDDistrib = n`, because then always the expected number of electrons is used (as if `EnforceECount = y`).

- `CutoffAutoAdjust` (optional, default = `n`):  
_Type: Switch; Allowed values: `yes` (`y`|`yes`|`t`|`true`) or `no` (`n`|`no`|`f`|`false`)._  
Determines whether the spatial and energetic cut-offs will be automatically adjusted after the equilibration (= `yes`) or not (= `no`). If `yes`, then at the end of the equilibration new cut-offs will be determined based on the maximum hop distance and the maximum state energy difference that occured during the equilibration, in combination with the `DistAdjustPercentage` and `EdiffAdjustPercentage` settings (see their descriptions). Discarding unused paths in this way can speed up the following main simulation. This option has no effect when no equilibration hops are specified (if `EqHops` $`= 0`$).

- `DistAdjustPercentage` (optional, default $`= 0.0`$):  
_Type: Real number; Unit: None; Allowed values: $`\geq 0.0`$._  
Percentage for the automatic adjustment of the spatial path cut-off after the equilibration (if `CutoffAutoAdjust = y`). The new spatial cut-off will be the maximum hop distance (that occurred during equilibration) times ($`1.0 +`$ `DistAdjustPercentage`$`/100`$), but only if the new value is smaller than the original cut-off and if each state would still have at least two paths.

- `EdiffAdjustPercentage` (optional, default $`= 0.0`$):  
_Type: Real number; Unit: None; Allowed values: $`\geq 0.0`$._  
Percentage for the automatic adjustment of the energetic path cut-off after the equilibration (if `CutoffAutoAdjust = y`). The new energetic cut-off will be the maximum state energy difference (of hops that occurred during equilibration) times ($`1.0 +`$ `EdiffAdjustPercentage`$`/100`$), but only if the new value is smaller than the original cut-off and if each state would still have at least two paths.

- `OnlyCompareSimID` (optional, default = `n`):  
_Type: Switch; Allowed values: `yes` (`y`|`yes`|`t`|`true`) or `no` (`n`|`no`|`f`|`false`)._  
Determines whether different simulation results should be compared only based on their `ProjectID` and `SimID` (= `yes`) or based on all simulation parameters (= `no`). This is relevant when already finished simulations (i.e. result files) are compared to those specified in an input file, in order to determine which simulations (or repetitions) are still missing and will be simulated, or in order to determine which results belong together when collecting averaged results in a table (`./VhoppAS -collect` command). This option can be used to treat simulations that differ slightly in the input parameters as still belonging together (e.g. when more repetitions with different hop counts are added). If `OnlyCompareSimID = yes` is used, then it must be ensured that only result files that belong together are in the output directory.

- `UseYZVariance` (optional, default = `y`):  
_Type: Switch; Allowed values: `yes` (`y`|`yes`|`t`|`true`) or `no` (`n`|`no`|`f`|`false`)._  
Determines whether the tracer diffusion coefficient $`D_{yz}`$ (perpendicular to the electric field) should be calculated from the $y$- and $z$-displacement variances (= `yes`) or from the mean squared displacements in $y$- and $z$-direction (= `no`). `UseYZVariance = yes` is necessary especially in the high-field regime, because spurious imbalances in the random structure (for example slightly more or longer paths in $`({+}x,{+}y)`$-direction than in $`({+}x,{-}y)`$-direction) can lead to an artificial small flux perpendicular to the high electric field.

- `ParallelizeReps` (optional, default = `n`):  
_Type: Switch; Allowed values: `yes` (`y`|`yes`|`t`|`true`) or `no` (`n`|`no`|`f`|`false`)._  
Controls the meaning of the `-job` command line parameter of VHoppAS. If `ParallelizeReps = no` then `-job` allows to select a certain parameter set (= a certain line in the `<VariedParameters>` table) for which all repetitions will be simulated sequentially. If `ParallelizeReps = yes` then `-job` allows to select a certain individual repetition that will solely be simulated (see description of `-job` further below).

- `Description` (optional, default = ""):  
_Type: Text; Allowed characters: Unrestricted (incl. line breaks = multiple lines are allowed as description)._  
Project description, e.g. for comments towards the purpose of the project. If specified, this setting has to be the last one in the `<MC-Project>` block.

**List of input parameters**:  
Each of the following parameters can be either present in the `<Parameters>`-block or in the `<VariedParameters>`-table. See also [technical notes](TECHNICAL_NOTES.md) for further definitions.
- `Repetitions` (optional, default $`= 1`$):  
_Type: Unsigned integer; Unit: None; Allowed values: $1$ - $100$._  
Specifies how often each parameter set shall be simulated. Averaging over several repetitions with different random structure (see `Seed` parameter) is important for meaningful results with error bars.

- `Emin` (mandatory):  
_Type: Real number; Unit: eV; Allowed values: Within energy axis in the DOS-File and lower than `Emax` parameter._  
Lower boundary $`E_{\text{min}}`$ of the DOS range to be used (= lowest allowed state energy; on the energy axis of the supplied DOS-File, i.e. relative to the charge-neutral Fermi level).

- `Emax` (mandatory):  
_Type: Real number; Unit: eV; Allowed values: Within energy axis in the DOS-File and higher than `Emin` parameter._  
Upper boundary $`E_{\text{max}}`$ of the DOS range to be used (= highest allowed state energy; on the energy axis of the supplied DOS-File, i.e. relative to the charge-neutral Fermi level).

- `ChemPot` (mandatory):  
_Type: Real number; Unit: eV; Allowed values: Unrestricted (= can even be outside DOS range), but must lead to at least one electron and at least one unoccupied state in the cell._  
Chemical potential $\mu$ of the electrons, i.e. the position of the used Fermi level on the DOS energy axis (= relative to the charge-neutral Fermi level). This parameter determines the number of electrons initially placed in the simulation cell. Under weak-field conditions, it corresponds to the actual (relative) chemical potential of the electrons, while in the high-field regime the effective chemical potential of electrons (= steady-state Fermi level) can differ due to a higher effective temperature of the electrons. `ChemPot` $`= 0`$ leads to simulation of the charge-neutral material, while `ChemPot` $`> 0`$ shifts the Fermi level to higher energy (= more electrons = negatively charged material) and `ChemPot` $`< 0`$ shifts the Fermi level to lower energy (= less electrons = positively charged material).

- `States` (mandatory):  
_Type: Unsigned integer; Unit: None; Allowed values: $200$ - $20000000$ (actual maximum depends on available memory)._  
Number of localized states $`N_{\text{st}}`$ (determines the cell volume).

- `MinPaths` (either this or `DistCutoff` must be specified and non-zero):  
_Type: Unsigned integer; Unit: None; Allowed values: $2$ - theoretical maximum for given cell volume (actual maximum depends on available memory)._  
Minimum number of hopping paths per localized state. The spatial cutoff gets determined by finding the radius that yields at least `MinPaths` paths for each state. This parameter is ignored, when `DistCutoff` is specified and non-zero.

- `DistCutoff` (either this or `MinPaths` must be specified and non-zero):  
_Type: Real number; Unit: nm; Allowed values: $0.0$ - half cell size (actual maximum depends on available memory)._  
Spatial path cutoff $`r_{\text{cut}}`$ (= maximum state-to-state distance considered for hopping). Is increased automatically if it does not lead to at least two paths for each state. If `DistCutoff` is specified and non-zero, then the `MinPaths` parameter is ignored.

- `EdiffCutoff` (optional, default $`= 0.0`$):  
_Type: Real number; Unit: eV; Allowed values: $0.0$ - width of the DOS range (Emax minus Emin)._  
Energetic path cutoff $`E_{\text{cut}}`$ (= maximum state energy difference considered for hopping). This constraint is applied after the spatial cutoff (and may lead to states with less paths than `MinPaths`). It is increased automatically if it does not leave at least two paths for each state. The value `EdiffCutoff` $`= 0`$ corresponds to no energetic cutoff (as if `EdiffCutoff` $`= E_{\text{max}} - E_{\text{min}}`$).

- `PreHops` (optional, default $`= 0`$):  
_Type: Unsigned integer; Unit: None; Allowed values: $0$ - $2000000000000$._  
Number of hops $`N_{\text{pre}}`$ during pre-equilibration (to reach steady-state electron distribution). No pre-equilibration if `PreHops` $`= 0`$.

- `EqHops` (optional, default $`= 0`$):  
_Type: Unsigned integer; Unit: None; Allowed values: $0$ - $2000000000000$._  
Number of hops $`N_{\text{eq}}`$ during equilibration (to reach proper displacement distributions and linear relation between displacement averages and time). No equilibration if `EqHops` $`= 0`$.

- `HopLimit` (mandatory):  
_Type: Unsigned integer; Unit: None; Allowed values: $200$ - $2000000000000$._  
Number of hops $`N_{\text{sim}}`$ during the main simulation (after pre-equilibration and equilibration).

- `Seed` (mandatory):  
_Type: Signed integer; Unit: None; Allowed values: Unrestricted._  
Seed for the initialization of the Mersenne-Twister random number generator (MT19337). If `Seed` $`> 0`$, then the random number generator gets initialized with the following actual seed: `Seed` $+$ `SimID`$*1000$ $+$ `RepID` (in order to ensure different seed for different parameter sets and different repetitions). If `Seed` $`\leq 0`$ then the absolute value of `Seed` is used in the equation above and additionally the current system time (clock ticks since year 1970) gets added to the seed.

- `AttTime` (mandatory):  
_Type: Real number; Unit: s; Allowed values: $`> 0.0`$ and $`< 1.0e100`$._  
Hop attempt time $`\nu_0^{-1}`$ (= inverse attempt frequency).

- `LocRadius` (mandatory):  
_Type: Real number; Unit: nm; Allowed values: $`> 0.0`$._  
Localization radius $\alpha$ (= decay length of localized wavefunctions).

- `Temp` (mandatory):  
_Type: Real number; Unit: K; Allowed values: $`\geq 0.01`$._  
Simulation temperature $T$ (= phonon temperature in the Miller-Abrahams rate equation; = electron temperature in the weak-field regime).

- `GradPhi` (mandatory):  
_Type: Real number; Unit: V/cm; Allowed values: Unrestricted._  
Electric potential gradient $F$ in x-direction. Its absolute value is the electric field strength. For `GradPhi` $`> 0`$, the electric field points in $`{-}x`$-direction and electrons flow in $`{+}x`$-direction (opposite for `GradPhi` $`< 0`$). For `GradPhi` $`= 0`$ there is no electric field in the simulation, i.e. only diffusion. 

When VHoppAS is started, it analyzes whether there are already finished simulations (i.e. result files) in the output directory that correspond to simulations specified in the input file. It then only conducts those simulations and repetitions that are not already finished. Only result files whose filename begins similar as specified in the `Output-File` setting are considered for this. Result files correspond to a certain parameter set when the `ProjectID`, the `SimID` and all simulation parameters (except `Repetitions` and `Seed`; also not the settings in the `<MC-Project>` block) are equal (see also `OnlyCompareSimID` setting). While finding these corresponding result files is the basis for aggregating result tables (`./VHoppAS -collect` command), it is also helpful in several other scenarios. For example, when VHoppAS was terminated (e.g. because runtime limit on a cluster was reached) after some but not all simulations or repetitions were finished, then it can be restarted with the same input file to simulate only what is missing. When all simulations or repetitions are finished, it is also possible to add more lines at the end of the `<VariedParameters>`-table or to increase the number of `Repetitions` when more simulations or repetitions are needed.

## Run simulations from the command line
In order to start VHoppAS simulations, a relative or absolute path to the input file has to be supplied as command line argument after `-input` (separated by white-space; if relative then relative to working directory). If the input file path contains white-spaces, it must be enclosed in double-quotes ("). It is recommended to start VHoppAS while in the directory that contains the input file (= working directory), because this facilitates the `DOS-File` and `Output-File` settings. For example, while in the Project001 directory of the example directory layout explained above, enter the command
```
../VhoppAS -input Input.txt
```
to sequentially run all simulations and repetitions specified in the input file `Input.txt`. It is recommended to redirect the quite substantial and informative logging output to a log file (output redirection with `&>`). Instead of running the simulations, it is also possible to just validate whether the input file has valid parameters by adding the `-validate` command line argument:
```
../VhoppAS -validate -input Input.txt
```
This command re-writes the input file, thereby adding omitted units and aligning the parameter table.

The `-job` command line argument can be added to execute only a specific simulation or repetition of the input file (instead of all). It has to be followed by an identification number `<ID>`, whose meaning depends on the `ParallelizeReps` setting in the input file. Example:
```
../VhoppAS -job 3 -input Input.txt
```
If `ParallelizeReps = no` (default), then `-job` allows to select a certain line in the `<VariedParameters>`-table by using the respective `SimID` as `<ID>` (first line = `1`, second line = `2` and so on). All repetitions of the selected parameter set are then simulated sequentially. If `ParallelizeReps = yes`, then `-job` allows to select a certain repetition. The `<ID>` is then the number of the repetition in the continous 1-based list of all repetitions in the input file, for example when the `<VariedParameters>`-table contains three or more lines (parameter sets) with `Repetitions = 5` for each, then `-job 12` would select the second repetition (`RepID = 2`) of the third parameter set (`SimID = 3`). Only this singular repetition would then be simulated. The `-job` functionality can be used to run simulations or repetitions in parallel (while VHoppAS itself is not parallelized), as demonstrated for example in the script further below.

VHoppAS can also aggregate the results of multi-type input files into tables. When `-job` is not used, i.e. when all simulations and repetitions of an input file are simulated sequentially, then this is done automatically at the end. Otherwise, when `-job` is used, then the different result files can be aggregated with this command:
```
../VhoppAS -collect -input Input.txt
```
This command instructs VHoppAS to find all corresponding result files in the output directory and match them to the parameter sets defined in the input file (as described above). It then creates two files, containing the results in tables: a summary file ("_Summary" appended to the `Output-File` name) that contains all individual results and a mean file ("_Mean" appended to the `Output-File` name) that for each `SimID` contains the average and standard deviation of the respective repetitions (the "sample standard deviation" is used; for some quantities the minimum and/or maximum are given instead of average and std. dev.; column prefix `d` for std. dev., column prefix `min` for minimum and column prefix `max` for maximum).

## Run simulations on a computer cluster with a job scheduler
VHoppAS itself is not parallelized, such that each VHoppAS execution is a single-core calculation. It is therefore recommended to use the `-job` functionality to distribute the different simulations (or even repetitions) of a project on multiple cores. For example, on a computer cluster with a SLURM job manager the simulations can be started as array jobs with the SimIDs supplied to the `#SBATCH --array` option as task IDs (e.g. `#SBATCH --array=1-44` to start individual calculations for all SimIDs from `1` to `44`). 

**Example job script** for a SLURM cluster (e.g. `job_script.sh` file; shebang, email-address and the loading of compiler runtime libraries must be adapted by user; assumes that VHoppAS executable resides in the parent directory; must be submitted from the directory where the input file resides):
<pre><code>#!&lt;<i>shebang</i>&gt;

# Submit to job queue by command: sbatch job_script.sh
# (submit from the working directory that contains the input file)

# specify mail notification
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=&lt;<i>email-address</i>&gt;

# specify runtime limit (d-hh:mm:ss)
#SBATCH --time=1-00:00:00

# ask for memory
#SBATCH --mem-per-cpu=3800M

# name the job
#SBATCH --job-name=hpc_001

# define sub-jobs (selected simulation IDs)
#SBATCH --array=1-5,8

# declare the merged STDOUT/STDERR file
#SBATCH --output=Result(%a)-%A.out

### begin of executable commands
# print job info
date +"%a %d.%m.%Y, %T"
echo "Starting job: ${SLURM_JOB_NAME}"
echo "Job array ID: ${SLURM_ARRAY_JOB_ID}"
echo "Job array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Directory: ${SLURM_SUBMIT_DIR}"

# load runtime libraries of used compiler
<i>module purge</i>
<i>module load GCC/11.3.0</i>
echo

# execute simulation
cd ${SLURM_SUBMIT_DIR}
../VHoppAS -job ${SLURM_ARRAY_TASK_ID} -input Input.txt

exitcode=$?
echo
echo "Exit code: ${exitcode}"
date +"%a %d.%m.%Y, %T"
exit ${exitcode}
</code></pre>
Afterwards, the results can be aggregated with `../VhoppAS -collect -input Input.txt`. When some tasks did not finish (e.g. due to runtime limit), the missing job IDs (either SimIDs or serialized repetition numbers depending on `ParallelizeReps`) are written to the summary and mean files (`Incomplete JobIDs: ...`) and can be re-supplied to `#SBATCH --array` to run the missing tasks.

## Output files and simulation results
For each individual simulation (= repetition), VHoppAS yields an individual result file (`Output-File` name with "(`SimID`.`RepID`)" appended) that contains the input settings (`<MC-Project>...</MC-Project>` section), the input parameters (`<Parameters>...</Parameters>` section, incl. the `SimID`) and the result values (`<Results>...</Results>` section, incl. the `RepID`). If `Verbosity` is not `0`, then the result file also contains convergence tables for the equilibration and the main simulation, showing the development of several quantities over the course of the equilibration and the main simulation (`<EquilibrationConvergence>` and `<SimulationConvergence>` sections). Additionally, each result file contains different histogram tables with additional statistics, as listed further below. 

**Result file layout** (e.g. `Result(31.3).txt` file):
<pre><code>&lt;MC-Project&gt;
... (<i>project settings from input file</i>)
&lt;/MC-Project&gt;

&lt;Parameters&gt;
SimID = 31
... (<i>parameter set from input file</i>)
&lt;/Parameters&gt;

&lt;Results&gt;
RepID = 3
... (<i>result quantities, see list below</i>)
&lt;/Results&gt;

&lt;EquilibrationConvergence&gt;
... (<i>development of quantities over the course of the equilibration</i>)
&lt;/EquilibrationConvergence&gt;

&lt;SimulationConvergence&gt;
... (<i>development of quantities during the main simulation, relative to the final values of the equilibration</i>)
&lt;/SimulationConvergence&gt;

... (<i>histograms, see list further below</i>)
</code></pre>

**List of result quantities**:  
See also [technical notes](TECHNICAL_NOTES.md) for further definitions.
- `Conductivity` (Unit: S/m):  
Conductivity $\sigma$ of the electrons.

- `Mobility` (Unit: cm<sup>2</sup>/Vs):  
Drift mobility $u$ of the electrons.

- `TracerDCoeff`, `TracerDCoeffP`, `TracerDCoeffT` (Unit: cm<sup>2</sup>/s):  
Tracer diffusion coefficient $D$, parallel tracer diffusion coefficient $`D_x`$ (`P` suffix; parallel to electric field) and transverse tracer diffusion coefficient $`D_{yz}`$ (`T` suffix; perpendicular to electric field). Could also be called self-diffusion coefficients (due to homogeneous KMC cell) and correspond to the respective chemical diffusion coefficients (due to the effective charge carrier formalism).

- `ChargeDCoeff` (Unit: cm<sup>2</sup>/s):  
Charge diffusion coefficient $`D_{\sigma}`$, i.e. the diffusion coefficient calculated from the conductivity or mobility with the Nernst-Einstein relation (using $`T_{\mathrm{eff}}`$ and $`\mu_{\mathrm{eff}}`$).

- `HavenRatio`, `HavenRatioP`, `HavenRatioT` (Unit: None):  
Haven ratio $`H_{\mathrm{R}}=D/D_{\sigma}`$, parallel Haven ratio $`H_{\mathrm{R},x}=D_x/D_{\sigma}`$ (`P` suffix) and transverse Haven ratio $`H_{\mathrm{R},yz}=D_{yz}/D_{\sigma}`$ (`T` suffix).

- `PartialEntropy` (Unit: eV/K):  
Partial entropy $`s_{\mathrm{e}}`$ of the electrons (for the specific random structure = calculated from the state energies $`E_i`$).

- `EffChemPot` (Unit: eV):  
Effective chemical potential $`\mu_{\mathrm{eff}}`$ of the electrons, as obtained from fitting the final electron distribution with the Fermi-Dirac function (if `TeffFit = y`). Is set to the input chemical potential (`ChemPot` parameter) if fitting is disabled (`TeffFit = n`).

- `EffTemp` (Unit: K):  
Effective temperature $`T_{\mathrm{eff}}`$ of the electrons, as obtained from fitting the final electron distribution with the Fermi-Dirac function (if `TeffFit = y`). Is set to the simulation temperature (`Temp` parameter) if fitting is disabled (`TeffFit = n`).

- `EffCarriers` (Unit: None):  
Number of effective charge carriers $`N_{\mathrm{eff}}`$ (for the specific random structure = calculated from the state energies $`E_i`$).

- `EffCarrierDensity` (Unit: 1/cm<sup>3</sup>):  
Effective charge carrier density $`n_{\mathrm{eff}}=N_{\mathrm{eff}}/V_{\mathrm{cell}}`$ (for the specific random structure = calculated from the state energies $`E_i`$ and the cell volume $`V_{\mathrm{cell}}`$).

- `Electrons`, `MobileElectrons`, `ZeroHopElectrons`, `OscElectrons` (Unit: None):  
Total numbers of electrons, of mobile electrons (= electrons with non-oscillating hops), of stationary electrons (= electrons with zero hops) and of oscillating electrons (= electrons with only oscillating hops) at the end of the simulation.

- `Time` (Unit: s):  
Simulated timespan $`\tau_{\mathrm{sim}}`$ (relative to the equilibration).

- `MeanDisp`, `MeanDispx`, `MeanDispy`, `MeanDispz` (Unit: nm):  
Mean displacement $`\langle R \rangle_{\mathrm{sim}}`$, mean `x`-displacement $`\langle X \rangle_{\mathrm{sim}}`$, mean `y`-displacement $`\langle Y \rangle_{\mathrm{sim}}`$ and mean `z`-displacement $`\langle Z \rangle_{\mathrm{sim}}`$ of effective charge carriers (relative to the equilibration).

- `MeanDispVarx`, `MeanDispVary`, `MeanDispVarz` (Unit: nm<sup>2</sup>):  
Variance of `x`-displacement $`\langle (X-\langle X \rangle_{\mathrm{mob}})^2 \rangle_{\mathrm{sim}}`$, variance of `y`-displacement $`\langle (Y-\langle Y \rangle_{\mathrm{mob}})^2 \rangle_{\mathrm{sim}}`$ and variance of `z`-displacement $`\langle (Z-\langle Z \rangle_{\mathrm{mob}})^2 \rangle_{\mathrm{sim}}`$ of effective charge carriers (relative to the equilibration).

- `MeanSqDispx`, `MeanSqDispy`, `MeanSqDispz` (Unit: nm<sup>2</sup>):  
Mean squared `x`-displacement $`\langle X^2 \rangle_{\mathrm{sim}}`$, mean squared `y`-displacement $`\langle Y^2 \rangle_{\mathrm{sim}}`$ and mean squared `z`-displacement $`\langle Z^2 \rangle_{\mathrm{sim}}`$ of effective charge carriers (relative to the equilibration).

- `NonOscHops`, `NonOscRatio` (Unit: None):  
Total number of non-oscillating hops during the main simulation and the ratio between the total number of non-oscillating hops and the total number of hops during the main simulation.

- `NonOscXDirRatio` (Unit: None):  
Ratio between the number of non-oscillating hops in $+x$-direction (= from state $i$ to state $j$ with $`x_{ij} = x_j - x_i \geq 0`$) and the total number of non-oscillating hops during the main simulation.

- `MeanFieldContrib` (Unit: eV):  
Average of the electric field contribution ($`-\mathrm{e}x_{ij}F`$ term with $F$ = electric potential gradient) of all hops that occurred during the main simulation.

- `TotalEnergy` (Unit: eV):  
Sum of the state energies $`E_i`$ of all occupied states (at the end of the simulation = final electron distribution).

- `CellSize` (Unit: nm):  
Edge length of the cubic KMC cell (= 3rd root of the cell volume).

- `MaxPathDist`, `MaxUsedPathDist` (Unit: nm):  
Actual spatial cut-off for hopping paths (= result of the automatic cut-off adjustment after the equilibration or equal to `DistCutoff` if no equilibration or no auto-adjust) and the maximum state-to-state distance of used paths (= maximum $`r_{ij}`$ of all hops that took place).

- `MaxPathEdiff`, `MaxUsedPathEdiff` (Unit: eV):  
Actual energetic cut-off for hopping paths (= result of the automatic cut-off adjustment after the equilibration or equal to `EDiffCutoff` if no equilibration or no auto-adjust) and the maximum absolute state energy difference of used paths (= maximum $`|E_{ij}|`$ of all hops that took place).

- `MaxPaths`, `MeanPaths` (Unit: None):  
Maximum number of paths per state (= number of paths of the localized state with the highest number of paths) and average number of paths per state.

- `MinUsedStateE`, `MaxUsedStateE` (Unit: eV):  
Minimum and maximum state energy $`E_i`$ of used states ("used" states are those that were part of a hop, i.e. those that have a non-zero number of outgoing or incoming hops).

- `ElectronsAboveEf`, `HolesBelowEf` (Unit: None):  
Number of occupied states above the Fermi level (= above $`\mu_{\mathrm{eff}}`$) and number of unoccupied states below the Fermi level (= below $`\mu_{\mathrm{eff}}`$) at the end of the simulation.

**List of histograms**:  
Histograms are generated as tables, where the first three columns correspond to the lower boundary, upper boundary and center of the bins. The ensuing columns then refer to these bins. The shown quantities should be rather self-explanatory (see sourcecode for details). All data in the histograms refers to equilibration and main simulation together (since the relevant variables are not reset inbetween for technical reasons). Each result file contains the following histograms:
- `Histogram:StateEnergies`:  
Distributions vs. **state energy** ($`E_i`$). For example which states were used, how many outgoing or incoming hops from which state energy and final energetic distribution of different electron categories (mobile, stationary, effective, ...). Note that `EffCarriers(calc)` is the correct effective charge carrier distribution, while `EffCarriers(a-c)` only serve as comparison to demonstrate that it is not possible to identify certain electrons as effective charge carriers.

- `Histogram:PathTimes`:  
Distributions vs. **hop time expectancy value for all possible paths** ($`\tau_{ij}`$ = inverse of Miller-Abrahams rate; not the randomized time). For example which paths were used and how many hops occured with a certain time value.

- `Histogram:HopTimes`:  
Distributions vs. **hop time expectancy value for the used paths** ($`\tau_{ij}`$ = inverse of Miller-Abrahams rate; not the randomized time). For example how many hops for which used paths.

- `Histogram:Displacements`:  
Distributions vs. individual **displacement** ($`R_k`$). For example how many electrons have a certain displacement.

- `Histogram:xDisplacements`:  
Distributions vs. individual **displacement in $x$-direction** ($`X_k`$), i.e. in direction of the electric potential gradient. For example how many electrons have a certain displacement in $x$-direction.

- `Histogram:yDisplacements`:  
Distributions vs. individual **displacement in $y$-direction** ($`Y_k`$). For example how many electrons have a certain displacement in $y$-direction.

- `Histogram:zDisplacements`:  
Distributions vs. individual **displacement in $z$-direction** ($`Z_k`$). For example how many electrons have a certain displacement in $z$-direction.

- `Histogram:HopCounts`:  
Distributions vs. **number of hops**. For example how many states have a certain number of outgoing hops or how many electrons made a certain number of hops.

- `Histogram:OscHopCounts`:  
Distributions vs. **number of oscillations** (= back-and-forth hops). For example how many states have a certain number of outgoing oscillations or how many electrons have a certain number of oscillations.

- `Histogram:NonOscHopCounts`:  
Distributions vs. **number of non-oscillating hops**. For example how many states have a certain number of outgoing non-oscillating hops or how many electrons made a certain number of non-oscillating hops.

- `Histogram:StateEnergyDifferences`:  
Distributions vs. **state energy difference** ($`E_{ij}`$ in the Miller-Abrahams equation; does not contain the contribution from electric field). For example how many paths exist with a certain energy difference between the partaking states or how many hops occurred with a certain state energy difference.

- `Histogram:FieldEnergyContributions`:  
Distributions vs. **energetic contribution of the electric field** ($`-\mathrm{e}x_{ij}F`$ term in the Miller-Abrahams equation). For example how many paths or how many hops occurred with a certain field contribution.

- `Histogram:Distances`:  
Distributions vs. **hop distance** ($`r_{ij}`$ in the Miller-Abrahams equation). For example how many hops occurred with a certain spatial distance between the states.

- `Histogram:NextTransitionTimes`:  
Distributions vs. each electron's **time until the next hop would occur** ($`\tau_{\text{min},k}`$). For example how many electrons have a certain $`\tau_{\text{min},k}`$ or how large are the $`\tau_{\text{min},k}`$ of stationary electrons at the end of the simulation.

- `Histogram:RelativeOccupationTimes`:  
Distributions vs. **fraction of time that a state has been occupied**. For example how many states are occupied all the time or how many states have never been occupied by an electron.