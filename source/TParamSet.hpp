#ifndef MC_TParamSet_H_
#define MC_TParamSet_H_

#include <cstdint>
#include <string>
#include <vector>
#include <array>

namespace MC
{

class TParamSet 
{
public:
    // Default constructor
    TParamSet();

	// Equality comparison
	bool operator==(const TParamSet& rhs) const;

	// Write all parameters
	std::string Write(bool is_output = true) const;

	// Write non-varied parameters
	std::string WriteConstant() const;

	// Write table header (description and unit of varied parameters as individual elements)
	std::vector<std::array<std::string,2>> WriteTableHeader(bool is_output = true) const;

	// Write table line (values of varied parameters as individual elements)
	std::vector<std::string> WriteTableLine(bool is_output = true) const;

	// Read parameters (true = successful read)
	bool Read(const std::string& general, const std::string& header = "", const std::string& line = "", bool raise_errors = true);

	// Identification number of this parameter set (for multiple simulations)
	std::uint32_t m_SimID;

	// Switch (copy from controller): 
    // true = adjust DOS energy axis to temperature
    // false = zero of DOS energy axis remains as in DOS file
	bool m_EFTAdjust;

	// Switch (copy from controller): 
    // true = distribute electrons initially according to Fermi-Dirac function
    // false = distribute electrons initially according to step function (lowest state energy first)
	bool m_InitialFDDistrib;

	// Switch (copy from controller): 
    // true = fit occupation for effective temperature and effective chemical potential
    // false = use input temperature and chemical potential for effective carriers
	bool m_TeffFit;

    // Switch (copy from controller): 
    // true = enforce that the number of electrons is exactly the number predicted from the F-D integral of random states
    // false = occupy all states through random number < F-D probability (leads to fluctuation of electron count)
	bool m_EnforceECount;

	// Switch (copy from controller): 
    // true = automatic adjustment of path cutoffs (spatial and energetic) after equilibration (used ranges + Constant::autocutoff_inc %)
    // false = no automatic adjustment of path cutoffs
    bool m_CutoffAutoAdjust;

	// Percentage for auto-adjust of distance cutoff (copy from controller; optional; 0.0 if not set)
    double m_DistCutoffAdjustPercentage;

    // Percentage for auto-adjust of energy difference cutoff (copy from controller; optional; 0.0 if not set)
    double m_EdiffCutoffAdjustPercentage;

	// Switch (copy from controller): 
    // true = compare parameter sets only based on SimID (ignoring differences in other parameters)
    // false = compare all parameters
    bool m_OnlyCompareSimID;

    // Switch (copy from controller): 
    // true = calculate Dyz from variances instead of mean squared displacements (only for non-zero field)
    // false = calculate Dyz from mean squared displacements
    bool m_UseYZVariance;

	// Number of repeated simulations (default: 1)
	std::uint32_t m_Repetitions;

    // Minimum possible state energy (in eV; corresponds to relative -1)
    double m_MinStateEnergy;

    // Maximum possible state energy (in eV; corresponds to relative 0)
    double m_MaxStateEnergy;

    // Chemical potential (in eV)
	double m_ChemPot;
    
    // Number of localized states
	std::uint32_t m_StateCount;

	// Minimum number of paths per state (when no predefined distance cutoff)
	std::uint32_t m_MinPathCount;

	// Path distance cutoff (in nm)
	double m_DistCutoff;

	// Path state energy difference cutoff (absolute; in eV)
	double m_EdiffCutoff;

	// Total number of pre-equilibration hops
	std::uint64_t m_PreHopLimit;

	// Total number of equilibration hops
	std::uint64_t m_EqHopLimit;

    // Total number of hops
	std::uint64_t m_HopLimit;

    // Random generator seed (positive: seed = value; negative or zero: seed = system time + abs(value))
    std::int64_t m_RndSeed;

	// Hopping-attempt time (in s; inverse attempt frequency)
	double m_AttemptTime;

	// Temperature (in K)
	double m_Temperature;

	// Electric potential gradient in +x-direction (in V/cm; equals negative electric field)
	double m_PhiGradient;

	// Localization radius (in nm)
	double m_LocRadius;

	// Parameter descriptors and units (without white-spaces)
	static const std::array<std::string,2> s_SimID;
	static const std::array<std::string,2> s_Repetitions;
	static const std::array<std::string,2> s_MinStateEnergy;
    static const std::array<std::string,2> s_MaxStateEnergy;
	static const std::array<std::string,2> s_ChemPot;
	static const std::array<std::string,2> s_StateCount;
	static const std::array<std::string,2> s_MinPathCount;
	static const std::array<std::string,2> s_DistCutoff;
	static const std::array<std::string,2> s_EdiffCutoff;
	static const std::array<std::string,2> s_PreHopLimit;
	static const std::array<std::string,2> s_EqHopLimit;
	static const std::array<std::string,2> s_HopLimit;
    static const std::array<std::string,2> s_RndSeed;
	static const std::array<std::string,2> s_AttemptTime;
	static const std::array<std::string,2> s_Temperature;
	static const std::array<std::string,2> s_PhiGradient;
	static const std::array<std::string,2> s_LocRadius;

	// Constant parameter flags (true = parameter is not varied for multiple simulations)
	bool c_Repetitions;
	bool c_MinStateEnergy;
    bool c_MaxStateEnergy;
	bool c_ChemPot;
	bool c_StateCount;
	bool c_MinPathCount;
	bool c_DistCutoff;
	bool c_EdiffCutoff;
	bool c_PreHopLimit;
	bool c_EqHopLimit;
	bool c_HopLimit;
    bool c_RndSeed;
	bool c_AttemptTime;
	bool c_Temperature;
	bool c_PhiGradient;
	bool c_LocRadius;
};

} // MC namespace
#endif  // MC_TParamSet_H_