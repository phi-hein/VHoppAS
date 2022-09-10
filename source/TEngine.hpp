#ifndef TEngine_H_
#define TEngine_H_

#include "TEngineData.hpp"
#include "TResult.hpp"

namespace MC
{

class TEngine : public TEngineData
{
public:	
	// Default constructor
	TEngine();

	// Fill the structure with local states
	void GenerateStructure();

	// Find hopping paths between local states and calculate hopping times
	void GeneratePaths();

	// Recalculate hopping times for all hopping paths (e.g. after changes of electric field or temperature)
	void UpdateHoppingTimes();

	// Fill the structure with electrons according to F-D distribution
	void GenerateElectrons();

	// Execute the simulation
	void RunSimulation();

	// Generate simulation results from simulation data
	void GenerateResults();

private:
	// Initialize the simulation data (e.g. randomized path times)
	void InitializeSimulation();

	// Execute KMC hops until hop counter equals the limit
	void KMCLoop(std::uint64_t& hop_counter, const std::uint64_t hop_limit) noexcept;

	// Generate statistics at end of equilibration (incl. optional auto-adjust of cut-offs)
	void CalculatePostEquilibrationStatistics();

	// Fit state occupation with Fermi-Dirac distribution (optional) and calculate effective charge carrier density
	void CalculateEffectiveCarrierDensity(bool fit_eff);

	// Extract result values from simulation data (optional: relative to reference result)
	void CalculateResultValues(const TResult* const ref_result);
};

} // MC namespace
#endif  // TEngine_H_