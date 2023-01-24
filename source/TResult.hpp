#ifndef MC_TResult_H_
#define MC_TResult_H_

#include <string>
#include <vector>
#include <ostream>
#include <memory>
#include <array>

#include "THistogram.hpp"

namespace MC
{

struct TProgress
{
	// Percentage of hop limit
	double m_Percentage;

	// Simulated hops
	std::uint64_t m_TotalHops;

	// Non-oscillatory hops
	std::uint64_t m_NonOscHops;

	// Ratio of non-oscillatory hops divided by total hops
	double m_NonOscHopRatio;

	// Simulated timespan (in s)
	double m_TotalTime;

	// Drift conductivity (in S/m = A/Vm)
	double m_DriftConductivity;

	// Drift mobility (in cm2/Vs; based on effective carrier density)
	double m_DriftMobility;

	// Diffusion coefficient (average of all directions; in cm2/s; based on effective carrier density)
	double m_DiffusionCoefficient;

	// Diffusion coefficient parallel to electric field (in cm2/s; based on effective carrier density)
	double m_DiffusionCoefficientParallel;

	// Diffusion coefficient perpendicular to electric field (in cm2/s; based on effective carrier density)
	double m_DiffusionCoefficientTransverse;

	// Partial entropy of electrons (in eV/K)
	double m_PartialEntropy;

	// Sum of occupied state energies (in eV)
	double m_TotalEnergy;

	// Effective chemical potential (in eV)
	double m_EffChemPot;

	// Effective temperature (in K)
	double m_EffTemp;

	// Effective carriers
	double m_EffCarriers;

	// Density of effective carriers (in 1/cm3)
	double m_EffCarrierDensity;

	// Electrons with > 0 non-oscillatory hops
	std::uint32_t m_MobileElectrons;

	// Electrons with zero hops
	std::uint32_t m_ZeroHopElectrons;

	// Electrons with only oscillatory hops
	std::uint32_t m_OscElectrons;

	// Mean displacements (in nm; per effective carrier)
	double m_MeanDisp;
	double m_MeanDisp_x;
	double m_MeanDisp_y;
	double m_MeanDisp_z;

	// Mean displacement variances (in nm2; per effective carrier)
	double m_MeanDispVariance_x;
	double m_MeanDispVariance_y;
	double m_MeanDispVariance_z;

	// Mean squared displacements (in nm2; per effective carrier)
	double m_MeanSquaredDisp_x;
	double m_MeanSquaredDisp_y;
	double m_MeanSquaredDisp_z;

	// Maximum distance of used paths (in nm)
    double m_MaxUsedPathDist;

	// Maximum absolute state energy difference of used paths (in eV)
    double m_MaxUsedPathEdiff;

	// Minimum and maximum used state energy (in eV)
	double m_MinUsedStateEnergy;
	double m_MaxUsedStateEnergy;

	// Electron-hole counts
	std::uint32_t m_ElectronsAboveEf;
	std::uint32_t m_HolesBelowEf;
};

class TResult 
{
public:
    // Result descriptors and units (without white-spaces)
	static const std::array<std::string,2> s_RepID;
	static const std::array<std::string,2> s_DriftConductivity;
	static const std::array<std::string,2> s_DriftMobility;
	static const std::array<std::string,2> s_DiffusionCoefficient;
	static const std::array<std::string,2> s_DiffusionCoefficientParallel;
	static const std::array<std::string,2> s_DiffusionCoefficientTransverse;
	static const std::array<std::string,2> s_HavenRatio;
	static const std::array<std::string,2> s_HavenRatioParallel;
	static const std::array<std::string,2> s_HavenRatioTransverse;
	static const std::array<std::string,2> s_PartialEntropy;
	static const std::array<std::string,2> s_EffChemPot;
	static const std::array<std::string,2> s_EffTemp;
	static const std::array<std::string,2> s_EffCarriers;
	static const std::array<std::string,2> s_EffCarrierDensity;
	static const std::array<std::string,2> s_ElectronCount;
	static const std::array<std::string,2> s_MobileElectrons;
	static const std::array<std::string,2> s_ZeroHopElectrons;
	static const std::array<std::string,2> s_OscElectrons;
	static const std::array<std::string,2> s_TotalTime;
	static const std::array<std::string,2> s_MeanDisp;
	static const std::array<std::string,2> s_MeanDisp_x;
	static const std::array<std::string,2> s_MeanDisp_y;
	static const std::array<std::string,2> s_MeanDisp_z;
	static const std::array<std::string,2> s_MeanDispVariance_x;
	static const std::array<std::string,2> s_MeanDispVariance_y;
	static const std::array<std::string,2> s_MeanDispVariance_z;
	static const std::array<std::string,2> s_MeanSquaredDisp_x;
	static const std::array<std::string,2> s_MeanSquaredDisp_y;
	static const std::array<std::string,2> s_MeanSquaredDisp_z;
	static const std::array<std::string,2> s_NonOscHops;
	static const std::array<std::string,2> s_NonOscHopRatio;
	static const std::array<std::string,2> s_InXDirNonOscRatio;
	static const std::array<std::string,2> s_MeanFieldContribution;
	static const std::array<std::string,2> s_TotalEnergy;
	static const std::array<std::string,2> s_CellSize;
	static const std::array<std::string,2> s_MaxPathDist;
	static const std::array<std::string,2> s_MaxUsedPathDist;
	static const std::array<std::string,2> s_MaxPathEdiff;
	static const std::array<std::string,2> s_MaxUsedPathEdiff;
	static const std::array<std::string,2> s_MaxPathCount;
	static const std::array<std::string,2> s_MeanPathCount;
	static const std::array<std::string,2> s_MinUsedStateEnergy;
	static const std::array<std::string,2> s_MaxUsedStateEnergy;
	static const std::array<std::string,2> s_ElectronsAboveEf;
	static const std::array<std::string,2> s_HolesBelowEf;

    // Default constructor
    TResult();

	// Copy constructor (ignores histograms)
	static std::unique_ptr<TResult> ValueCopy(const TResult* const result);

	// Write results block
	void Write(std::ostream& o_str) const;

	// Write table header (as individual elements)
	static std::vector<std::array<std::string,2>> WriteTableHeader();

	// Write table line (as individual elements)
	std::vector<std::string> WriteTableLine() const;

	// Read results (true = successful read)
	bool Read(const std::string& str, bool raise_errors = true);

	// Write analysis for results of multiple simulations
    static void WriteMultiAnalysis(std::ostream& o_str, const std::vector<std::unique_ptr<const TResult>>& results);

	// Write analysis for results of multiple simulations
    static void WriteMultiAnalysis(std::ostream& o_str, const std::vector<std::vector<std::unique_ptr<const TResult>>>& results);

	// Write table header for mean and stddev of multiple simulations
	static std::vector<std::array<std::string,2>> WriteMultiTableHeader();

	// Write table line for mean and stddev of multiple simulations
	static std::vector<std::string> WriteMultiTableLine(const std::vector<std::unique_ptr<const TResult>>& results);

	// Save progress of equilibration or simulation
	void SaveProgress(double percentage, bool is_eq);

	// Write progress table header
	void WriteProgressHeader(std::ostream& o_str, std::uint64_t maxhops, bool has_field) const;

	// Write last progress line of equilibration or simulation
	void WriteProgressLine(std::ostream& o_str, std::uint64_t maxhops, bool has_field, bool is_eq) const;

	// Write convergence table of equilibration or simulation
	void WriteConvergenceTable(std::ostream& o_str, bool is_eq) const;

	// Simulated hops (not written to or loaded from file)
	std::uint64_t m_TotalHops;

	// Repetition-ID
	std::uint32_t m_RepID;

	// Drift conductivity (in S/m = A/Vm)
	double m_DriftConductivity;

	// Drift mobility (in cm2/Vs; based on effective carrier density)
	double m_DriftMobility;

	// Diffusion coefficient (average of all directions; in cm2/s; based on effective carrier density)
	double m_DiffusionCoefficient;

	// Diffusion coefficient parallel to electric field (in cm2/s; based on effective carrier density)
	double m_DiffusionCoefficientParallel;

	// Diffusion coefficient perpendicular to electric field (in cm2/s; based on effective carrier density)
	double m_DiffusionCoefficientTransverse;

	// Haven ratio (average of all directions)
	double m_HavenRatio;

	// Haven ratio parallel to electric field (Dx / Dsigma)
	double m_HavenRatioParallel;

	// Haven ratio transverse to electric field (Dyz / Dsigma = Haven ratio at weak field)
	double m_HavenRatioTransverse;

	// Partial entropy of electrons (in eV/K)
	double m_PartialEntropy;

	// Effective chemical potential (in eV)
	double m_EffChemPot;

	// Effective temperature (in K)
	double m_EffTemp;

	// Effective carriers
	double m_EffCarriers;

	// Density of effective carriers (in 1/cm3)
	double m_EffCarrierDensity;

	// Total number of electrons
	std::uint32_t m_ElectronCount;

	// Electrons with > 0 non-oscillatory hops
	std::uint32_t m_MobileElectrons;

	// Electrons with zero hops
	std::uint32_t m_ZeroHopElectrons;

	// Electrons with only oscillatory hops
	std::uint32_t m_OscElectrons;

	// Simulated timespan (in s)
	double m_TotalTime;

	// Mean displacements (in nm; per effective carrier)
	double m_MeanDisp;
	double m_MeanDisp_x;
	double m_MeanDisp_y;
	double m_MeanDisp_z;

	// Mean displacement variances (in nm2; per effective carrier)
	double m_MeanDispVariance_x;
	double m_MeanDispVariance_y;
	double m_MeanDispVariance_z;

	// Mean squared displacements (in nm2; per effective carrier)
	double m_MeanSquaredDisp_x;
	double m_MeanSquaredDisp_y;
	double m_MeanSquaredDisp_z;

	// Non-oscillatory hops
	std::uint64_t m_NonOscHops;

	// Ratio of non-oscillatory hops divided by total hops
	double m_NonOscHopRatio;

	// Ratio of hops in +x direction divided by total hops
	double m_InXDirNonOscRatio;

	// Mean energy contibution of electric field to executed hops (in eV)
	double m_MeanFieldContribution;

	// Sum of occupied state energies (in eV)
	double m_TotalEnergy;

	// Cell size (in nm)
	double m_CellSize;

	// Maximum distance of paths (in nm)
    double m_MaxPathDist;

	// Maximum distance of used paths (in nm)
    double m_MaxUsedPathDist;

	// Maximum absolute state energy difference of paths (in eV)
    double m_MaxPathEdiff;

	// Maximum absolute state energy difference of used paths (in eV)
    double m_MaxUsedPathEdiff;

	// Maximal number of paths per state
	std::uint32_t m_MaxPathCount;

	// Average number of paths per state
	double m_MeanPathCount;

	// Minimum and maximum used state energy (in eV)
	double m_MinUsedStateEnergy;
	double m_MaxUsedStateEnergy;

	// Electron-hole counts
	std::uint32_t m_ElectronsAboveEf;
	std::uint32_t m_HolesBelowEf;

	// Equilibration progress
	std::vector<TProgress> m_EqProgress;

	// Simulation progress
	std::vector<TProgress> m_Progress;

	// Histogram: State energies (in eV)
	THistogram m_HStateEnergy;

	// Histogram: Pre-calculated path times (in s; on logarithmic scale)
	THistogram m_HPathTime;

	// Histogram: Pre-calculated path times used in hops (in s; on logarithmic scale)
	THistogram m_HHopTime;

	// Histograms: Electron displacements (in nm)
	THistogram m_HDisp;
	THistogram m_HDisp_x;
	THistogram m_HDisp_y;
	THistogram m_HDisp_z;
	
	// Histograms: Electron hop count, osc. count and non-osc. count (on logarithmic scale)
	THistogram m_HHopCount;
	THistogram m_HOscHopCount;
	THistogram m_HNonOscHopCount;

	// Histogram: State energy difference of paths and hops (in eV)
	THistogram m_HStateEnergyDifference;

	// Histogram: Field energy contribution of paths and hops (in eV)
	THistogram m_HFieldEnergyContribution;

	// Histogram: Spatial distance of paths and hops (in nm)
	THistogram m_HDistance;

	// Histogram: Next transition time (current minimal path time for each electron incl. random factor; in s; on logarithmic scale)
	THistogram m_HNextTime;

	// Histogram: Occupied timespan relative to simulated timespan
	THistogram m_HRelOccTime;
};

} // MC namespace
#endif  // MC_TResult_H_