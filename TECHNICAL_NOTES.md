# Technical notes

This document provides a step-by-step description of the KMC-VRH simulations implemented in VHoppAS, including the underlying maths. For more explanations see our respective paper and for more details see the sourcecode (especially the `TEngine.cpp` file, which contains most of the actual simulation code).

## Fundamental concepts
The simulations are based on the [Miller-Abrahams rate equation](https://doi.org/10.1103/PhysRev.120.745 "A. Miller and E. Abrahams, Physical Review, 1960, 120, 745–755.")
```math
\nu_{ij} = \left\{
\begin{array}{ll}
\nu_0\exp\left(-\frac{2|\vec{r}_{ij}|}{\alpha}-\frac{E_{ij}-\mathrm{e}\vec{r}_{ij}\vec{F}}{k_{\text{B}}T}\right) & \text{for } E_{ij}-\mathrm{e}\vec{r}_{ij}\vec{F} > 0 \\
\nu_0\exp\left(-\frac{2|\vec{r}_{ij}|}{\alpha}\right) & \, \text{for } E_{ij}-\mathrm{e}\vec{r}_{ij}\vec{F} \leq 0 \\
\end{array}\right.
```
which describes the expected hopping rate $\nu_{ij}$ from an occupied localized state $i$ to an unoccupied localized state $j$, depending on the distance $`\vec{r}_{ij}=\vec{r}_{j}-\vec{r}_{i}`$ and the energy difference $`E_{ij}=E_{j}-E_{i}`$ between the states. The involved parameters are the attempt frequency $\nu_0$, the localization radius $\alpha$, the electric potential gradient $\vec{F}=\vec{\nabla}\varphi$ and the temperature $T$ ($\mathrm{e}$ is the elementary charge, $k_{\text{B}}$ is the Boltzmann constant).

Furthermore, the simulations are based on the effective charge carrier formalism. In general, only a fraction of the electrons in a system, i.e. those close to the Fermi level, can contribute to charge transport at a time (e.g. due to the Pauli exclusion principle). Since not the same electrons act as charge carriers all the time (e.g. due to trapping/de-trapping in deeper states) and since electrons are indistinguishable, the electron movement observed in the simulation must be re-interpreted as the movement of effective charge carriers, whose amount is given by the [effective charge carrier density](https://doi.org/10.1063/1.4871757 "V. Palenskis, AIP Advances, 2014, 4, 047119.")
```math
n_{\text{eff}} = \int_{\infty}^{-\infty} g(E)f_{\text{FD}}(E)\left[1-f_{\text{FD}}(E)\right]\,\mathrm{d}E
```
where $g(E)$ is the density of states (DOS) and $f_{\text{FD}}(E)$ is the Fermi-Dirac distribution.

## Spatial and energetic axes
The simulated KMC cell that contains the localized states is cubic with periodic boundary conditions in the three directions $x$, $y$ and $z$. The electric potential gradient is oriented in $x$-direction according to 
```math
\vec{F} = \vec{\nabla}\varphi = \left(
\begin{array}{ccc}
F \\ 0 \\ 0 
\end{array}
\right)
```
such that $|F|$ is the electric field strength and for $F > 0$ the electrons flow in $x$-direction ($F$ is the `GradPhi` input value). The cell is considered homogeneous in the sense of spatially constant input parameters, for example constant temperature ($\vec{\nabla}T=0$), constant electric field ($\vec{\nabla}^2\varphi=0$) and position-independent DOS ($\vec{\nabla}g(E)=0$).

The energy axis $E$ is supplied with the DOS (in the DOS input file). $E=0$ there corresponds to the position of the reference Fermi level, i.e. the Fermi level of the charge-neutral material at a certain reference temperature $T_{\text{ref}}$. Depending on the `EFTAdjust` input setting, this reference Fermi level gets adjusted to the simulation temperature $T$ by finding $\Delta{}E$ (bisection search algorithm) such that
```math
\int_{-\infty}^{\infty}g(E)\left[\exp\left(\frac{E}{k_{\text{B}}T_{\text{ref}}}\right)+1\right]^{-1}\,\mathrm{d}E
= \int_{-\infty}^{\infty}g(E)\left[\exp\left(\frac{E+\Delta{}E}{k_{\text{B}}T}\right)+1\right]^{-1}\,\mathrm{d}E
```
```math
\quad(\text{or}\quad
\int_{-\infty}^{0}g(E)\,\mathrm{d}E
= \int_{-\infty}^{\infty}g(E)\left[\exp\left(\frac{E+\Delta{}E}{k_{\text{B}}T}\right)+1\right]^{-1}\,\mathrm{d}E
\quad\text{if}\quad
T_{\text{ref}} = 0\text{ K})
```
and shifting the energy axis according to $E \rightarrow E + \Delta{}E$, such that $E=0$ then corresponds to the charge-neutral Fermi level at the simulation temperature $T$ (under the assumption of a temperature-independent DOS). The actual Fermi level in the simulations may differ from this reference Fermi level in order to simulate the material in a charged state (= more or less electrons than normal). Note that in the KMC cell the Fermi level (= the electrochemical potential of electrons) is not spatially constant due to the linearly changing electric potential $\varphi(x)$, but it still has a constant position on the energy axis $E$ of the DOS since the state energies are similarly shifted by the electric potential. $E$ here always corresponds to energies without the contribution from the electric potential. In contrast to the Fermi level and the electric potential, the chemical potential of the electrons is constant throughout the KMC cell. A certain "position of the Fermi level on the axis $E$" therefore corresponds to a certain chemical potential of the electrons. It is then possible to define a relative chemical potential $\mu$ (relative to the chemical potential of the charge-neutral material) such that $E=\mu$ is the position of the Fermi level throughout the KMC cell (`ChemPot` input value). $\mu \neq 0$ then leads to simulations with a charged material.

Since only the region around the Fermi level is relevant for hopping (due to the simulaneous presence of occupied and unoccupied states), only localized states in a chosen energy range $`[E_{\text{min}},E_{\text{max}}]`$ on the axis $E$ are regarded (`Emin` and `Emax` input values; in most cases symmetric around `ChemPot`). For a given total number of localized states $N_{\text{st}}$ (`States` input value) the volume of the simulated KMC cell is then
```math
V_{\text{cell}} = N_{\text{st}}\left[\int_{E_{\text{min}}}^{E_{\text{max}}}g(E)\,\mathrm{d}E\right]^{-1}
```
Consequently, the edge length of the cubic cell is $`V_{\text{cell}}^{1/3}`$ (`CellSize` output value).

## Relative internal quantities
The edge length of the cell is used as spatial scaling factor $`f_{\text{sp}}=V_{\text{cell}}^{1/3}`$ for all internal quantities related to spatial dimensions. The relative spatial axes $\hat{x}$, $\hat{y}$ and $\hat{z}$ that confine the cell then range from $0$ to $1$, such that the position of a localized state $i$ is given by
```math
\vec{r}_i = f_{\text{sp}}\hat{r}_i
\quad\text{with}\quad
\hat{r}_i = \left(
\begin{array}{ccc}
\hat{x}_i \\ \hat{y}_i \\ \hat{z}_i 
\end{array}
\right)
\quad\text{and}\quad
\hat{x}_i \in [0,1], \hat{y}_i \in [0,1], \hat{z}_i \in [0,1]
```
For a hopping path from state $i$ to state $j$ this yields $`\vec{r}_{ij}=\vec{r}_j-\vec{r}_i=f_{\text{sp}}(\hat{r}_j-\hat{r}_i)=f_{\text{sp}}\hat{r}_{ij}`$. Note that the periodic boundary conditions must be considered for these state-to-state vectors, meaning that $`\vec{r}_{ij}`$ and $`\hat{r}_{ij}`$ correspond to the shortest route from state $i$ to state $j$, which can be either within the cell or through any of the periodic boundaries. 

Similarly, an energetic scaling factor $`f_{\text{en}}=E_{\text{max}}-E_{\text{min}}`$ is introduced for all internal quantities related to energies. The relative energy axis $\hat{E}$ ranges from $-1$ to $0$, such that the energy of a localized state $i$ is given by
```math
E_i = E_{\text{max}} + f_{\text{en}}\hat{E}_i
\quad\text{with}\quad
\hat{E}_i \in [-1,0]
```
For a hopping path from state $i$ to state $j$ this yields $`E_{ij}=E_j-E_i=f_{\text{en}}(\hat{E}_j-\hat{E}_i)=f_{\text{en}}\hat{E}_{ij}`$. The position of the Fermi level on the relative energy axis $\hat{E}$ is $`\hat{\mu}=(\mu-E_{\text{max}})/f_{\text{en}}`$.

The relative values of further input quantities are:
- Inverse localization radius ($\alpha$ is the `LocRadius` input value): $`\hat{\alpha}^{-1} = f_{\text{sp}}/\alpha`$
- Electric potential gradient ($F$ is the `GradPhi` input value): $`\hat{F} = \mathrm{e}F f_{\text{sp}}/f_{\text{en}}`$
- Thermal energy ($T$ is the `Temp` input value): $`\hat{E}_{\text{th}} = k_{\text{B}}T/f_{\text{en}}`$ 
- Spatial cutoff ($r_{\text{cut}}$ is the `DistCutoff` input value): $`\hat{r}_{\text{cut}} = r_{\text{cut}}/f_{\text{sp}}`$
- Energetic cutoff ($E_{\text{cut}}$ is the `EdiffCutoff` input value): $`\hat{E}_{\text{cut}} = E_{\text{cut}}/f_{\text{en}}`$

**Time-to-hop expectancy value** with relative quantities (= inverse Miller-Abrahams rate equation; $\nu_0^{-1}$ is the `AttTime` input value):
```math
\tau_{ij} = \nu_{ij}^{-1} = \left\{
\begin{array}{ll}
\nu_0^{-1}\exp\left(2|\hat{r}_{ij}|\hat{\alpha}^{-1}\right)\exp\left(\frac{\hat{E}_{ij}-\hat{x}_{ij}\hat{F}}{\hat{E}_{\text{th}}}\right) & \text{for } \hat{E}_{ij}-\hat{x}_{ij}\hat{F} > 0 \\
\nu_0^{-1}\exp\left(2|\hat{r}_{ij}|\hat{\alpha}^{-1}\right) & \, \text{for } \hat{E}_{ij}-\hat{x}_{ij}\hat{F} \leq 0 \\
\end{array}\right.
```

**Fermi-Dirac distribution function** with relative quantities:
```math
f_{\text{FD}}(\hat{E}) = \left[\exp\left(\frac{\hat{E}-\hat{\mu}}{\hat{E}_{\text{th}}}\right)+1\right]^{-1}
```

## Preparations before the simulation
Note: VHoppAS is based on the assumption of spin-degeneracy (= localized states do not differ in energy depending on the electron spin + spin-up and spin-down electrons behave equally), such that only one spin-type (= half of the states and half of the electrons) need to be simulated and count twice for the results. This is the reason for the factors $1/2$ and $2$ that occur occasionally in the following sections. It is convention that all input and output values refer to both spin-types together. See also the respective comment in the [user manual](USER_MANUAL.md).

### 1. Generate the localized states
Create $`\hat{N}_{\text{st}} = N_{\text{st}}/2`$ localized states with
- uniform random position:  
$\hat{x}_i$, $\hat{y}_i$, $\hat{z}_i$ = uniform random numbers in the interval $[0,1)$.  
(`std::uniform_real_distribution` of the C++ standard library)
- DOS-weighted random energy:  
$\hat{E}_i$ = random number in the interval $[-1,0)$ with the normalized DOS $g(\hat{E})$ as distribution function.  
(`std::piecewise_linear_distribution` of the C++ standard library)

### 2. Find all hopping paths within the cutoffs
If no spatial cutoff was explicitly specified (if `DistCutoff` $=0$ or not defined in input file), then determine it from the `MinPaths` input value:
- Find the spatial cutoff $`\hat{r}_{\text{cut}}`$ for which all localized states $i$ have at least `MinPaths` paths to other states $j$ with $`|\hat{r}_{ij}| < \hat{r}_{\text{cut}}`$.

If a spatial cutoff was explicitly specified (if `DistCutoff` $\neq 0$), then validate it to ensure at least two paths per state:
- Find the minimum spatial cutoff $`\hat{r}_{\text{cut}}^{*}`$ for which all localized states $i$ would have at least two paths to other states $j$ with $`|\hat{r}_{ij}| < \hat{r}_{\text{cut}}^{*}`$. If $`\hat{r}_{\text{cut}} < \hat{r}_{\text{cut}}^{*}`$ then use $`\hat{r}_{\text{cut}}^{*}`$ as new $`\hat{r}_{\text{cut}}`$.

The resulting spatial cutoff must not exceed half the cell size ($`\hat{r}_{\text{cut}}`$ must be smaller than $0.5$) to avoid having two different paths between the same two states (through the periodic boundaries).

If no energetic cutoff was explicitly specified (if `EdiffCutoff` $=0$ or not defined in input file), then set it to the full range:
- Set $`\hat{E}_{\text{cut}}=1`$, which is equivalent to no energetic cutoff.

If an energetic cutoff was explicitly specified (if `EdiffCutoff` $\neq 0$), then validate it to ensure at least two paths per state:
- Find the minimum energetic cutoff $`\hat{E}_{\text{cut}}^{*}`$ for which all localized states $i$ would have at least two paths to other states $j$ with both $`|\hat{r}_{ij}| < \hat{r}_{\text{cut}}`$ and $`|\hat{E}_{ij}| < \hat{E}_{\text{cut}}^{*}`$. If $`\hat{E}_{\text{cut}} < \hat{E}_{\text{cut}}^{*}`$ then use $`\hat{E}_{\text{cut}}^{*}`$ as new $`\hat{E}_{\text{cut}}`$.

Finally, for each localized state $i$ create all outgoing hopping paths which satisfy both $`|\hat{r}_{ij}| < \hat{r}_{\text{cut}}`$ and $`|\hat{E}_{ij}| < \hat{E}_{\text{cut}}`$.

### 3. Calculate the time-to-hop expectancy values
For each outgoing path of each state, calculate the time-to-hop expectancy value $\tau_{ij}$ (see equation above), which corresponds to the time that it takes **on average** until an electron on state $i$ would hop to the unoccupied state $j$. The $\tau_{ij}$ values do not change during the simulation because the positions and energies of the states do not change.

### 4. Populate the localized states with electrons
Expected number of electrons in the cell (rounded to integer; for $`f_{\text{FD}}(\hat{E}_i)`$ see the Fermi-Dirac equation further above):
```math
\hat{N}_{\text{e}}^{\text{Exp}} = \sum_i^{\hat{N}_{\text{st}}}f_{\text{FD}}(\hat{E}_i)
```

If `InitialFDDistrib = yes`, then place on each state $i$ an electron if a uniform random number $RN \in [0,1)$ is smaller than the respective Fermi-Dirac probability $`f_{\text{FD}}(\hat{E}_i)`$. The resulting number of Fermi-Dirac-distributed electrons is called $`\hat{N}_{\text{e}}`$. If $`\hat{N}_{\text{e}}`$ differs from $`\hat{N}_{\text{e}}^{\text{Exp}}`$ and if `EnforceECount = yes`, then re-evaluate randomly chosen states with a new $RN$ until the number of electrons $`\hat{N}_{\text{e}}`$ matches the expectation value $`\hat{N}_{\text{e}}^{\text{Exp}}`$.

If `InitialFDDistrib = no`, then distribute $`\hat{N}_{\text{e}} = \hat{N}_{\text{e}}^{\text{Exp}}`$ electrons among the localized states from lowest to highest state energy, which leads to a step function as initial electron distribution.

Since each internal electron counts twice, for example the `Electrons` output value corresponds to $`N_{\text{e}} = 2\hat{N}_{\text{e}}`$. For comparison, the expectation value for $N_{\text{e}}$ is also calculated from the DOS by numerical integration (rounded to integer):
```math
N_{\text{e}}^{\text{DOS}} = V_{\text{cell}}\int_{-1}^{0}g(\hat{E})f_{\text{FD}}(\hat{E})\,\mathrm{d}\hat{E}
```

### 5. Initialize the electrons with randomized hop times
For each electron $k$:
- Retrieve the list of its possible hopping paths from the state it occupies.
- For each path calculate a randomized time-to-hop $`\tau_{ij}^{*} = -\ln(RN)\tau_{ij}`$ with a uniform random number $RN \in [0,1)$. This equation for $`\tau_{ij}^{*}`$ originates from the assumption that hops are Poisson-distributed rare events with expectation value $\tau_{ij}$. If the respective target state $j$ is occupied by another electron, then set $`\tau_{ij}^{*} = \infty`$ (i.e. the hop will never happen). 
- Find the lowest and the 2nd lowest $`\tau_{ij}^{*}`$ of this electron (= $`\tau_{\text{min},k}`$ and $`\tau_{\text{2nd},k}`$).

Rationale: $`\tau_{\text{min},k}`$ is the path that the electron $k$ would take next (after the time $`\tau_{\text{min},k}`$ has passed), unless another electron hops to the respective target state beforehand and blocks this path. Simultaneously determining $`\tau_{\text{2nd},k}`$ is useful for example for the following situation. If an electron is oscillating between two states $i$ and $j$ and say the state $i$ is energetically favorable ($`\hat{E}_i \ll \hat{\mu}`$), then there is a high chance that the state $i$ is also the target state of the $`\tau_{\text{min},k}`$ of several surrounding electrons. Each time the occupation of $i$ switches due to the oscillating electron, a new $`\tau_{\text{min},k}`$ of these surrounding electrons would need to be searched. This can be avoided when their $`\tau_{\text{2nd},k}`$ is known: When $`\tau_{\text{min},k}`$ to state $i$ gets blocked, then $`\tau_{\text{2nd},k}`$ becomes the new $`\tau_{\text{min},k}`$ (and $`\tau_{\text{2nd},k}`$ is set to "undefined"). When state $i$ becomes unoccupied again (the oscillating electron leaves), then most likely the respective time-to-hop $\tau^{*}$ is lower than $\tau_{\text{min},k}$ (because $`\hat{E}_i`$ is so low). In this case $\tau_{\text{min},k}$ becomes $\tau_{\text{2nd},k}$ (as it was before) and the path to $i$ becomes $\tau_{\text{min},k}$ again, thereby restoring the original lowest and second lowest paths for the surrounding electron without the need to search them anew. Note that this is only one of the cases where $\tau_{\text{2nd},k}$ is useful (see also the KMC algorithm below).

## KMC algorithm
The KMC algorithm simulates the movement of the electrons based on random numbers. The most important quantities that are recorded during the algorithm are:
- Simulated timespan $\tau$: This variable records the progress of time during the movement of the electrons.
- Electron displacement vector $(\hat{X}_k, \hat{Y}_k, \hat{Z}_k)$: These variables record each electrons' individual movement in the direction of the three spatial axes $\hat{x}$, $\hat{y}$ and $\hat{z}$.

These values are initialized to zero at the beginning of the simulation. Several more statistics are recorded during the simulation (e.g. how many hops each electron makes or how often a certain path is used), but discussing them here would be too lengthy (see `KMCLoop` function in `TEngine.cpp` for this).

### Basic principle
The basic principle of the here employed KMC algorithm is the following:
1. Find the electron that hops next, i.e. the electron with the overall lowest time $`\tau_{\text{hop}} = \min(\tau_{\text{min},k})`$.
2. Execute the respective hop: Transfer the electron to the target state and add the movement to its displacement vector.
3. Advance the time: Add $\tau_{\text{hop}}$ to the simulated timespan $\tau$ and subtract $\tau_{\text{hop}}$ from all other $\tau_{ij}^{*}$ (incl. all $\tau_{\text{min},k}$).
4. Assign new $\tau_{ij}^{*}$ to the electron that hopped and to previously or newly blocked paths (incl. update the $\tau_{\text{min},k}$ of all affected electrons).
5. Repeat from 1. until the specified hop count is reached.

The actual algorithm is more complicated due to optimizations and the reader is referred to the heavily annotated `KMCLoop` function in the `TEngine.cpp` file for all details. Here, only the use of second lowest times shall be explained in some detail. Additional to the $\tau_{\text{min},k}$ and $\tau_{\text{2nd},k}$ that were introduced for each electron $k$ further above, also the global variables $\tau_{\text{min}}$ and $\tau_{\text{2nd}}$ are introduced that refer to the electron that would hop next and to the electron that would hop after it (more specifically $\tau_{\text{min}}$ and $\tau_{\text{2nd}}$ are the $\tau_{\text{min},k}$ of these electrons). Actually, each $\tau_{\text{min},k}$ and $\tau_{\text{2nd},k}$ is accompanied by an index that links it to the respective path, and $\tau_{\text{min}}$ and $\tau_{\text{2nd}}$ are accompanied by an index that links them to the respective electrons. In this way, for each of these time values also the electron and the specific path that is meant is known (here we just use the time values and do not explicitly introduce these indices for ease of reading). In consequence, each of these time values may either point to a valid path, which can be either non-blocked ($<\infty$) or blocked ($=\infty$), or it may be "undefined" (point to no electron or path). An exception to this is that the $\tau_{\text{min},k}$ must never be "undefined". At the beginning of the simulation, both $\tau_{\text{min}}$ and $\tau_{\text{2nd}}$ are "undefined" (while $\tau_{\text{min},k}$ and $\tau_{\text{2nd},k}$ are already initialized as explained further above). During the KMC loop, the variable $\tau_{\text{hop}}$ (see algorithm above) refers to the electron that currently hops (and the path it takes) and it is ensured that this is always a valid non-blocked path.

### Usage of $\tau_{\text{min}}$ and $\tau_{\text{2nd}}$
$\tau_{\text{min}}$ acts as a supplier for $\tau_{\text{hop}}$ (and $\tau_{\text{2nd}}$ as a supplier for $\tau_{\text{min}}$) at the beginning of each KMC loop iteration:
- If $\tau_{\text{min}}$ is not already known (= if $\tau_{\text{min}}$ is "undefined"), iterate over all electrons $k$ to find the lowest and 2nd lowest $\tau_{\text{min},k}$, and set $\tau_{\text{min}}$ and $\tau_{\text{2nd}}$ accordingly.
- Use $\tau_{\text{min}}$ as $\tau_{\text{hop}}$, use $\tau_{\text{2nd}}$ as new $\tau_{\text{min}}$ (= next hop candidate) and set $\tau_{\text{2nd}}$ to "undefined".

It is ensured in the algorithm that at the beginning of the KMC loop $\tau_{\text{min}}$ never refers to a blocked path.

During the KMC loop iteration, in principle the following rules apply whenever any $\tau_{\text{min},k}$ changes, i.e. whenever any electron $k$ gets an updated lowest time-to-hop (how this happens is discussed in the next section). These rules only apply when $\tau_{\text{min}}$ is not "undefined" (because then there would be nothing to compare to):
- If the new $\tau_{\text{min},k}<\infty$ (= valid non-blocked path):  
    (special rules apply when $\tau_{\text{min},k}$ previously was either $\tau_{\text{min}}$ or $\tau_{\text{2nd}}$, not shown here)
    - If $\tau_{\text{min},k}<\tau_{\text{min}}$ (= new global lowest time-to-hop):
        - Use $\tau_{\text{min}}$ as new $\tau_{\text{2nd}}$ and set $\tau_{\text{min},k}$ as new $\tau_{\text{min}}$.
    - Else if $\tau_{\text{2nd}}$ is not "undefined" and $\tau_{\text{min},k}<\tau_{\text{2nd}}$ (= new global 2nd lowest time-to-hop):
        - Set $\tau_{\text{min},k}$ as new $\tau_{\text{2nd}}$.
- If the new $\tau_{\text{min},k}=\infty$ (= all paths of the electron are blocked):
    - If $\tau_{\text{min},k}$ previously was $\tau_{\text{min}}$ (= global lowest time-to-hop became blocked):
        - Use $\tau_{\text{2nd}}$ as new $\tau_{\text{min}}$ and set $\tau_{\text{2nd}}$ to "undefined".
    - If $\tau_{\text{min},k}$ previously was $\tau_{\text{2nd}}$ (= global 2nd lowest time-to-hop became blocked):
        - Set $\tau_{\text{2nd}}$ to "undefined".

Three constraints must always be respected additional to the above rules:
- $\tau_{\text{min}}$ must refer to a different path than $\tau_{\text{hop}}$ (not necessarily a different electron)
- $\tau_{\text{min}}$ and $\tau_{\text{2nd}}$ must refer to different electrons
- $\tau_{\text{min}}$ and $\tau_{\text{2nd}}$ must lead to different target states (otherwise $\tau_{\text{min}}$ could block $\tau_{\text{2nd}}$)

### Keeping $\tau_{\text{min},k}$ and $\tau_{\text{2nd},k}$ up to date
During each KMC loop iteration, there are three situations that influence the $\tau_{\text{min},k}$ and $\tau_{\text{2nd},k}$:

Firstly, when the electron that currently hops was transferred to the target state, it gets new $`\tau_{ij}^{*} = -\ln(RN)\tau_{ij}`$ for all paths outgoing from this state (or $`\tau_{ij}^{*}=\infty`$ for those that are blocked). The lowest and the 2nd lowest of these new $\tau_{ij}^{*}$ become $\tau_{\text{min},k}$ and $\tau_{\text{2nd},k}$, respectively, of this electron.

Secondly, the hop that takes place switches the initial state from occupied to unoccupied, thereby enabling all paths to this state. All surrounding electrons that can reach this newly empty state each get a new $`\tau_{ij}^{*} = -\ln(RN)\tau_{ij} < \infty`$ for the respective path (it was $`\tau_{ij}^{*}=\infty`$ for the blocked path before). Each new $\tau_{ij}^{*}$ is then a potential candidate for $\tau_{\text{min},k}$ or $\tau_{\text{2nd},k}$ of the respective electron:
- If $\tau_{ij}^{*}<\tau_{\text{min},k}$ (= new lowest time-to-hop):
    - Use $\tau_{\text{min},k}$ as new $\tau_{\text{2nd},k}$ and set $\tau_{ij}^{*}$ as new $\tau_{\text{min},k}$.
- Else if $\tau_{\text{2nd},k}$ is not "undefined" and $\tau_{ij}^{*}<\tau_{\text{2nd},k}$ (= new 2nd lowest time-to-hop):
    - Set $\tau_{ij}^{*}$ as new $\tau_{\text{2nd},k}$.

Thirdly, the hop that takes place switches the target state from unoccupied to occupied, thereby disabling all paths to this state. All surrounding electrons that could reach this newly blocked state each get a new $`\tau_{ij}^{*} = \infty`$ for the respective path (it was $`\tau_{ij}^{*}<\infty`$ for the available path before):
- If this $\tau_{ij}^{*}$ previously was $\tau_{\text{min},k}$ (= lowest time-to-hop became blocked):
    - If $\tau_{\text{2nd},k}$ is not "undefined":
        - Use $\tau_{\text{2nd},k}$ as new $\tau_{\text{min},k}$ and set $\tau_{\text{2nd},k}$ to "undefined".
    - Else (if $\tau_{\text{2nd},k}$ was "undefined"):
        - Iterate through all $\tau_{ij}^{*}$ of this electron to find $\tau_{\text{min},k}$ and $\tau_{\text{2nd},k}$.
- If this $\tau_{ij}^{*}$ previously was $\tau_{\text{2nd},k}$ (= 2nd lowest time-to-hop became blocked):
    - Set $\tau_{\text{2nd},k}$ to "undefined".

Since it does not affect the order of the $\tau_{\text{min},k}$, changing $\tau_{\text{min},k}$ and $\tau_{\text{2nd},k}$ when the time advances is avoided (contrary to what is said in the algorithm principle above), by storing for each electron the time it hopped last. This avoids the subtraction of $\tau_{\text{hop}}$ from all $\tau_{ij}^{*}$ (incl. the $\tau_{\text{min},k}$ and $\tau_{\text{2nd},k}$) and it avoids unnecessary updates of $\tau_{\text{min}}$ and $\tau_{\text{2nd}}$.

## Pre-equilibration, equilibration and the main simulation
Each simulation is structured into three sequential steps, which do not differ in the algorithm but in the way the resulting statistics are used:

1. **Pre-equilibration** (optional):  
The KMC algorithm introduced above is executed until $`\hat{N}_{\text{pre}} = N_{\text{pre}}/2`$ hops are reached ($N_{\text{pre}}$ is the `PreHops` input value). Afterwards, all recorded statistics are discarded, e.g. the simulated timespan and the electron displacements are reset to zero. The only purpose of this pre-equilibration step is that the electrons can redistribute and adopt a steady-state distribution.

2. **Equilibration** (optional):  
The KMC algorithm is executed for $`\hat{N}_{\text{eq}} = N_{\text{eq}}/2`$ more hops ($N_{\text{eq}}$ is the `EqHops` input value). Afterwards, the reached simulated timespan and all displacement averages are saved as reference values for the following main simulation (see next section about the calculation of the results), but no values are reset. The purpose of this equilibration step is that the electron displacements $\hat{X}_k$, $\hat{Y}_k$ and $\hat{Z}_k$ can reach continuously developing distributions, where the displacement averages increase linear with the simulated timespan.

3. **Main simulation**:  
The KMC algorithm is executed for $`\hat{N}_{\text{sim}} = N_{\text{sim}}/2`$ more hops ($N_{\text{sim}}$ is the `HopLimit` input value). The movement of the electrons during this period (relative to the end of the equilibration) determines the simulation results (see next section).

## Calculation of result quantities
### Effective charge carriers
Before any simulation results can be calculated from the simulated timespan $\tau$ and the electron displacement vectors $`(\hat{X}_k, \hat{Y}_k, \hat{Z}_k)`$, the number of effective charge carriers $N_{\text{eff}}$ must be determined (`EffCarriers` output value):
```math
N_{\text{eff}} = 2 \sum_i^{\hat{N}_{\text{st}}}f_{\text{FD}}(E_i,\mu_{\text{eff}},T_{\text{eff}})[1-f_{\text{FD}}(E_i,\mu_{\text{eff}},T_{\text{eff}})]
\quad\text{with}\quad
f_{\text{FD}}(E,\mu_{\text{eff}},T_{\text{eff}}) = \left[\exp\left(\frac{E-\mu_{\text{eff}}}{k_{\text{B}}T_{\text{eff}}}\right)+1\right]^{-1}
```
The effective charge carrier density is then $`n_{\text{eff}} = N_{\text{eff}}/V_{\text{cell}}`$ (`EffCarrierDensity` output value). For comparison, also the values expected from the DOS are calculated and reported in the log output:
```math
N_{\text{eff}}^{\text{DOS}} = V_{\text{cell}} n_{\text{eff}}^{\text{DOS}}
\quad\text{and}\quad
n_{\text{eff}}^{\text{DOS}} = \int_{\infty}^{-\infty} g(E)f_{\text{FD}}(E,\mu_{\text{eff}},T_{\text{eff}})\left[1-f_{\text{FD}}(E,\mu_{\text{eff}},T_{\text{eff}})\right]\,\mathrm{d}E
```
These equations require the knowledge of the effective temperature $T_{\text{eff}}$ (`EffTemp` output value) and the effective chemical potential $\mu_{\text{eff}}$ (`EffChemPot` output value) of the electrons (= the actual temperature and Fermi level position of the steady-state electron distribution). If `TeffFit = no`, then simply $T_{\text{eff}} = T$ (`Temp` input value) and $\mu_{\text{eff}} = \mu$ (`ChemPot` input value) are used. If `TeffFit = yes`, then the electron distribution is fitted with a Fermi-Dirac function by finding the minimum of
```math
\sum_i^{\hat{N}_{\text{st}}}\left(w_i - f_{\text{FD}}(E_i,\mu_{\text{eff}},T_{\text{eff}})\right)^2
```
with regard to $T_{\text{eff}}$ through the robust bisection search algorithm ($w_i = 1$ when state $i$ is occupied and $w_i = 0$ if the state is unoccupied; there is exactly one minimum such that convergence is guaranteed when starting with temperatures on opposite sides of the minimum). The search only has to be performed along one degree of freedom (here $T_{\text{eff}}$ was chosen) because $T_{\text{eff}}$ and $\mu_{\text{eff}}$ are not independent. Since the total electron density is constant, for each tested $T_{\text{eff}}$ within this fitting procedure the corresponding $\mu_{\text{eff}}$ is determined similar to the $T$-adjustment of the Fermi level discussed further above (= also via bisection search with guaranteed convergence).

The partial entropy $s_{\text{e}}$ of the electrons (`PartialEntropy` output value), which is related to the effective charge carrier density, is calculated in a comparable way (for a derivation of this quantity see the supplement of our paper):
```math
s_{\text{e}} = k_{\text{B}}\frac{N_{\text{therm}}}{N_{\text{eff}}}
\quad\text{with}\quad
N_{\text{therm}} = 2 \sum_i^{\hat{N}_{\text{st}}}f_{\text{FD}}(E_i,\mu_{\text{eff}},T_{\text{eff}})[1-f_{\text{FD}}(E_i,\mu_{\text{eff}},T_{\text{eff}})]\frac{E_i-\mu_{\text{eff}}}{k_{\text{B}}T_{\text{eff}}}
```
```math
s_{\text{e}}^{\text{DOS}} = k_{\text{B}}\frac{N_{\text{therm}}^{\text{DOS}}}{N_{\text{eff}}^{\text{DOS}}}
\quad\text{with}\quad
N_{\text{therm}}^{\text{DOS}} = V_{\text{cell}}\int_{\infty}^{-\infty} g(E)f_{\text{FD}}(E,\mu_{\text{eff}},T_{\text{eff}})\left[1-f_{\text{FD}}(E,\mu_{\text{eff}},T_{\text{eff}})\right]\frac{E-\mu_{\text{eff}}}{k_{\text{B}}T_{\text{eff}}}\,\mathrm{d}E
```

### Electron and hop categories
Additional to the total hop counter, also for each electron it is counted how often it hops and for each path it is counted how often it was used (= individual hop counters). Thereby, the hops are distinguished into "oscillatory" and "non-oscillatory" hops as explained in the following:

Each time an electron hops from one state to another and then back to the previous state, this is an **oscillation**. For example, a certain electron $k$ may hop from state $i$ to state $j$ and the next hop this electron makes takes it back to state $i$. It is irrelevant whether other electrons moved inbetween these two hops, i.e. the two hops do not have to be immediately after each other in the algorithm (just for the certain electron). The two hops "$`i \rightarrow j`$" and "$`j \rightarrow i`$" together count as an oscillation and each of these hops is termed "oscillatory hop". This means that these two hops not only count as $+2$ for the individial hop counter of the electron (and the total hop counter), but also as $+2$ for the individual oscillatory hop counter of the electron (and the total oscillatory hop counter). The individual oscillatory hop counters of the two paths "$`i \rightarrow j`$" and "$`j \rightarrow i`$" are each increased by $+1$ (like their individual hop counters). Important: The counting of these oscillations is based on pairs, meaning that when in this example the electron would now again hop from $i$ to $j$ this would not count as oscillatory hop even though it previously hopped from $j$ to $i$ (because this was part of the previous oscillation), unless it would then also hop back to $i$ afterwards, thereby completing the next oscillation. In summary, each hop is
- **oscillatory** if it is part of an oscillation
- **non-oscillatory** if it not part of an oscillation

Only the non-oscillatory hops contribute to each electrons' displacement. The oscillatory hops are still an inherent part of the physics of the system and cannot be neglected in general, for example because they control which paths are blocked or temporarily non-blocked.

The electrons are also distinguished into categories based on their behavior. Each electron can be either "mobile", "stationary" or "oscillating", according to the following definitions:
- **mobile**: The electron has a non-zero number of non-oscillatory hops.
- **stationary**: The electron has zero hops (= did not move at all).
- **oscillating**: The electron has only oscillatory hops (= non-zero number of hops but no non-oscillatory hops).

Counting the number of mobile electrons in the cell yields $`\hat{N}_{\text{mob}}`$. Since each internal electron counts twice (to get the value for spin-up and spin-down electrons together), their total number is $`N_{\text{mob}} = 2 \hat{N}_{\text{mob}}`$ (`MobileElectrons` output value). Similarly, the numbers of stationary electrons (`ZeroHopElectrons` output value) and oscillatory electrons (`OscElectrons` output value) are determined. 

In the log output and in the histograms appear further sub-categories of the mobile electrons:
- "Zero-displacement electrons" (`ZeroDispElectrons`): Mobile (!) electrons that still have zero displacement, which means they returned to their initial position through a series of non-oscillatory hops.
- "Electrons with a single hop": Mobile electrons with just a single hop (which is necessarily non-oscillatory).
- "Electrons with a single hop except osc.": Mobile electrons that have a non-zero number of oscillatory hops and just a single non-oscillatory hop.
- "Electrons with only non-osc. hops": Mobile electrons that have only non-oscillatory hops (and no oscillations).

### Displacement averages
The individual displacements $`(\hat{X}_k, \hat{Y}_k, \hat{Z}_k)`$ of the electrons $k$ are aggregated into the following displacement averages. Therein, each internal electron counts twice (for the two spin-types) and the movement of the electrons is attributed to the movement of effective charge carriers. Stationary and oscillating electrons are excluded since they always have no displacement and do not contribute to transport.

**Mean displacements of mobile electrons**:
```math
\langle X\rangle_{\text{mob}} = \frac{2f_{\text{sp}}}{N_{\text{mob}}}\sum_k^{\hat{N}_{\text{mob}}}\hat{X}_k
\;,\quad
\langle Y\rangle_{\text{mob}} = \frac{2f_{\text{sp}}}{N_{\text{mob}}}\sum_k^{\hat{N}_{\text{mob}}}\hat{Y}_k
\quad\text{and}\quad
\langle Z\rangle_{\text{mob}} = \frac{2f_{\text{sp}}}{N_{\text{mob}}}\sum_k^{\hat{N}_{\text{mob}}}\hat{Z}_k
```

**Mean displacements of effective charge carriers**:
```math
\langle R\rangle = \frac{2f_{\text{sp}}}{N_{\text{eff}}}\sum_k^{\hat{N}_{\text{mob}}}\left(\hat{X}_k^2 + \hat{Y}_k^2 + \hat{Z}_k^2\right)^{1/2}
\;,\quad
\langle X\rangle = \frac{2f_{\text{sp}}}{N_{\text{eff}}}\sum_k^{\hat{N}_{\text{mob}}}\hat{X}_k
\;,\quad
\langle Y\rangle = \frac{2f_{\text{sp}}}{N_{\text{eff}}}\sum_k^{\hat{N}_{\text{mob}}}\hat{Y}_k
\quad\text{and}\quad
\langle Z\rangle = \frac{2f_{\text{sp}}}{N_{\text{eff}}}\sum_k^{\hat{N}_{\text{mob}}}\hat{Z}_k
```

**Mean squared displacements of effective charge carriers**:
```math
\langle X^2\rangle = \frac{2f_{\text{sp}}^2}{N_{\text{eff}}}\sum_k^{\hat{N}_{\text{mob}}}\hat{X}_k^2
\;,\quad
\langle Y^2\rangle = \frac{2f_{\text{sp}}^2}{N_{\text{eff}}}\sum_k^{\hat{N}_{\text{mob}}}\hat{Y}_k^2
\quad\text{and}\quad
\langle Z^2\rangle = \frac{2f_{\text{sp}}^2}{N_{\text{eff}}}\sum_k^{\hat{N}_{\text{mob}}}\hat{Z}_k^2
```

**Displacement variances of effective charge carriers**:
```math
\langle(X-\langle X\rangle_{\text{mob}})^2\rangle 
= \langle X^2\rangle - \frac{N_{\text{eff}}}{N_{\text{mob}}}\langle X\rangle^2
\;,\quad
\langle(Y-\langle Y\rangle_{\text{mob}})^2\rangle 
= \langle Y^2\rangle - \frac{N_{\text{eff}}}{N_{\text{mob}}}\langle Y\rangle^2
\quad\text{and}\quad
\langle(Z-\langle Z\rangle_{\text{mob}})^2\rangle 
= \langle Z^2\rangle - \frac{N_{\text{eff}}}{N_{\text{mob}}}\langle Z\rangle^2
```

When an equilibration step is used (`EqHops` $\neq 0$), then the displacement averages of the main simulation are relative to the values at the end of the equilibration (if no equilibration then they just correspond to the above equations). For this, the above averages are calculated at the end of the equilibration (here the notation `eq` is used for these values) and used as reference for the results of the main simulation (here the notation `sim` is used for relative values). If `TeffFit = yes` then the number of effective charge carriers might differ slightly between the end of the equilibration ($N_{\text{eff}}^{\text{eq}}$) and the end of the main simulation ($N_{\text{eff}}$), because the electron distribution and therefore the fitted $T_{\text{eff}}$ and $\mu_{\text{eff}}$ might differ slightly. This is corrected for by adjusting the averages from the equilibration with the factor $N_{\text{eff}}^{\text{eq}}/N_{\text{eff}}$.

**Mean displacements of effective charge carriers** (`MeanDisp`, `MeanDispx`, `MeanDispy` and `MeanDispz` output values):
```math
\langle R\rangle_{\text{sim}} = \langle R\rangle - \frac{N_{\text{eff}}^{\text{eq}}}{N_{\text{eff}}}\langle R\rangle_{\text{eq}}
\;,\quad
\langle X\rangle_{\text{sim}} = \langle X\rangle - \frac{N_{\text{eff}}^{\text{eq}}}{N_{\text{eff}}}\langle X\rangle_{\text{eq}}
\;,\quad
\langle Y\rangle_{\text{sim}} = \langle Y\rangle - \frac{N_{\text{eff}}^{\text{eq}}}{N_{\text{eff}}}\langle Y\rangle_{\text{eq}}
\quad\text{and}\quad
\langle Z\rangle_{\text{sim}} = \langle Z\rangle - \frac{N_{\text{eff}}^{\text{eq}}}{N_{\text{eff}}}\langle Z\rangle_{\text{eq}}
```

**Mean squared displacements of effective charge carriers** (`MeanSqDispx`, `MeanSqDispy` and `MeanSqDispz` output values):
```math
\langle X^2\rangle_{\text{sim}} = \langle X^2\rangle - \frac{N_{\text{eff}}^{\text{eq}}}{N_{\text{eff}}}\langle X^2\rangle_{\text{eq}}
\;,\quad
\langle Y^2\rangle_{\text{sim}} = \langle Y^2\rangle - \frac{N_{\text{eff}}^{\text{eq}}}{N_{\text{eff}}}\langle Y^2\rangle_{\text{eq}}
\quad\text{and}\quad
\langle Z^2\rangle_{\text{sim}} = \langle Z^2\rangle - \frac{N_{\text{eff}}^{\text{eq}}}{N_{\text{eff}}}\langle Z^2\rangle_{\text{eq}}
```

**Displacement variances of effective charge carriers** (`MeanDispVarx`, `MeanDispVary` and `MeanDispVarz` output values):
```math
\langle(X-\langle X\rangle_{\text{mob}})^2\rangle_{\text{sim}} = \langle(X-\langle X\rangle_{\text{mob}})^2\rangle - \frac{N_{\text{eff}}^{\text{eq}}}{N_{\text{eff}}}\langle(X-\langle X\rangle_{\text{mob}})^2\rangle_{\text{eq}}
\;,\quad
\langle(Y-\langle Y\rangle_{\text{mob}})^2\rangle_{\text{sim}} = \langle(Y-\langle Y\rangle_{\text{mob}})^2\rangle - \frac{N_{\text{eff}}^{\text{eq}}}{N_{\text{eff}}}\langle(Y-\langle Y\rangle_{\text{mob}})^2\rangle_{\text{eq}}
```
```math
\quad\text{and}\quad
\langle(Z-\langle Z\rangle_{\text{mob}})^2\rangle_{\text{sim}} = \langle(Z-\langle Z\rangle_{\text{mob}})^2\rangle - \frac{N_{\text{eff}}^{\text{eq}}}{N_{\text{eff}}}\langle(Z-\langle Z\rangle_{\text{mob}})^2\rangle_{\text{eq}}
```

Similarly, the **simulated timespan** must be related to the equilibration (`Time` output value): $`\tau_{\text{sim}} = \tau - \tau_{\text{eq}}`$

### Electron transport properties
Based on the aforementioned quantities, the following major result values can be calculated:

**Conductivity** (`Conductivity` output value):
```math
\sigma = \mathrm{e}n_{\text{eff}}\frac{\langle X \rangle_{\text{sim}}}{F \tau_{\text{sim}}}
```
If there is no electric field (`GradPhi` $=0$), then $\sigma$ is set to zero.

**Drift mobility** (`Mobility` output value):
```math
u = -\frac{\langle X \rangle_{\text{sim}}}{F \tau_{\text{sim}}} = -\frac{\sigma}{\mathrm{e}n_{\text{eff}}}
```
If there is no electric field (`GradPhi` $=0$), then $u$ is set to zero.

**Tracer diffusion coefficients** (`TracerDCoeff`, `TracerDCoeffP` and `TracerDCoeffT` output values): 
```math
D = \frac{D_x + 2 D_{yz}}{3}
\;,\quad
D_x = \frac{\langle(X-\langle X\rangle_{\text{mob}})^2\rangle_{\text{sim}}}{2\tau_{\text{sim}}}
\quad\text{and}\quad
D_{yz} 
= \frac{\langle Y^2 \rangle_{\text{sim}} + \langle Z^2 \rangle_{\text{sim}}}{4\tau_{\text{sim}}}
```
As explained in the supplement of our paper, these tracer diffusion coefficients equal the chemical diffusion coefficients, because they contain the thermodynamic factor through the use of effective charge carriers. If `UseYZVariance = yes`, then $D_{yz}$ is calculated with the variances $`\langle(Y-\langle Y\rangle_{\text{mob}})^2\rangle_{\text{sim}}`$ and $`\langle(Z-\langle Z\rangle_{\text{mob}})^2\rangle_{\text{sim}}`$ instead of the mean squared displacements. If there is no electric field (`GradPhi` $=0$), then both $D_x$ and $D_{yz}$ are calculated with mean squared displacements, i.e. with $`\langle X^2 \rangle_{\text{sim}}`$ and $`\langle Y^2 \rangle_{\text{sim}} + \langle Z^2 \rangle_{\text{sim}}`$, respectively (irrespective of `UseYZVariance`).

**Charge diffusion coefficient** (`ChargeDCoeff` output value):  
```math
D_{\sigma} = -\frac{k_{\text{B}}T_{\text{eff}}u}{\mathrm{e}}
```
If there is no electric field (`GradPhi` $=0$), then $D_{\sigma}$ is set to zero.

**Haven ratios** (`HavenRatio`, `HavenRatioP` and `HavenRatioT` output values):
```math
H_{\mathrm{R}}=\frac{D}{D_{\sigma}}
\;,\quad
H_{\mathrm{R},x}=\frac{D_x}{D_{\sigma}}
\quad\text{and}\quad
H_{\mathrm{R},yz}=\frac{D_{yz}}{D_{\sigma}}
```
From a theoretical point of view, the Haven ratios should equal unity (because the tracer diffusion coefficients equal the chemical diffusion coefficients), such that these values should scatter around unity and serve more as a consistency check. Deviations might occur for example in the case of anormalous diffusion (= non-linear relation between displacement variance and time), e.g. for $D_x$ at high field strength. If there is no electric field (`GradPhi` $=0$), then $H_{\mathrm{R}}$, $H_{\mathrm{R},x}$ and $H_{\mathrm{R},yz}$ are set to zero.

All these output values are not only determined at the end of the main simulation but also in 2.5 % intervals over the course of equilibration and main simulation (if `Verbosity` $\geq 1$). This everytime includes a fitting of the electron distribution if `TeffFit = yes`. For the equilibration, the displacement averages and the simulated timespan are of course not relative but the total values (e.g. $\langle X \rangle$ and $\tau$ used to calculate $\sigma$). In the output file, the `<Results>` block refers to the final values at the end of the main simulation, while the convergence tables (and the log output) contain also the intermediate values. If `TeffFit = yes`, then the repetitive fitting typically leads to slightly differing numbers of effective charge carriers for each intermediate result. This causes an additional noise for the intermediate quantities based on $N_{\text{eff}}$. While the progress table in the log output contains these variations (since it is written live during the simulation), in both convergence tables of the output file all affected quantities are corrected to the final $N_{\text{eff}}$ at the end of the main simulation (by multiplication with the factor $`N_{\text{eff},i}/N_{\text{eff}}`$ for each intermediate result $i$). This means for example that the last line of the equilibration progress table in the log output contains the $\langle X\rangle_{\text{eq}}$ value, while the corresponding last line of the `<EquilibrationConvergence>` table in the output file contains the corrected value $`\langle X\rangle_{\text{eq}}N_{\text{eff}}^{\text{eq}}/N_{\text{eff}}`$.
