#ifndef MC_TEngineData_H_
#define MC_TEngineData_H_

#include <vector>
#include <random>
#include <memory>
#include <cstdint>

#include "Constants.hpp"
#include "TDOS.hpp"
#include "TParamSet.hpp"
#include "TResult.hpp"
#include "TLocalState.hpp"
#include "TElectron.hpp"

namespace MC
{

class TEngineData 
{
public:	
	// Default constructor
	TEngineData ();

	// Specify output verbosity
	void SetVerbosity(Verbosity vl);

	// Specify density of states
	void SetDOS(std::unique_ptr<TDOS> dos);

	// Validate parameters
	static void ValidateParameters(const TParamSet& params);

	// Specify parameters
	void SetParameters(std::unique_ptr<const TParamSet> params, std::uint32_t rep_id, std::uint64_t seed_increment = 0);

	// Return results object (removes histograms from this object, generate again if needed)
	void GetResults(std::unique_ptr<const TResult>& result);

	// Print diagnostics data (e.g. for debug)
	void PrintDiagnostics() const;

protected:
	// Verbosity level
	Verbosity m_VL;

    // Density of states
	std::unique_ptr<TDOS> m_DOS;

	// Simulation parameter set
	std::unique_ptr<const TParamSet> m_ParamSet;

    // Ready switch: true = DOS and parameters are properly defined
    bool m_ParametersReady;

	// Parameter: Thermal energy (kBT; scaled to relative energetic dimension)
	double m_ThermalEnergy;

	// Parameter: Electric potential energy gradient in +x-direction (scaled to relative spatial and energetic dimensions)
	double m_GradPhiEnergy;

	// Parameter: Inverse localization radius (scaled to relative spatial dimensions)
	double m_InvLocRadius;

	// SimData: Random number generator
	std::mt19937_64 m_RndGen;

	// SimData: Localized states
	std::vector<TLocalState> m_Structure;

	// SimData: Electrons
	std::vector<TElectron> m_Electrons;

	// SimData: Simulated timespan (in s)
	double m_TotalTime;

	// Ready switch: true = simulation finished properly
    bool m_SimulationReady;

	// Simulation results
	std::unique_ptr<TResult> m_SimResult;

	// Ready switch: true = simulation finished properly and results are calculated
    bool m_ResultsReady;
};

} // MC namespace
#endif  // MC_TEngineData_H_