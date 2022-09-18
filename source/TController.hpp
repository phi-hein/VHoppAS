#ifndef TController_H_
#define TController_H_

#include <string>
#include <vector>
#include <memory>
#include <filesystem>
#include <ostream>
#include <cstdint>

#include "Constants.hpp"
#include "TDOS.hpp"
#include "TParamSet.hpp"
#include "TResult.hpp"

namespace MC
{

class TController
{
public:
    // Parameter descriptors
    static const std::string s_ProjectID;
    static const std::string s_ProjectName;
    static const std::string s_DOSFile;
    static const std::string s_OutputFile;
    static const std::string s_VL;
    static const std::string s_EFTAdjust;
    static const std::string s_InitialFDDistrib;
    static const std::string s_TeffFit;
    static const std::string s_EnforceECount;
    static const std::string s_CutoffAutoAdjust;
    static const std::string s_DistCutoffAdjustPercentage;
    static const std::string s_EdiffCutoffAdjustPercentage;
    static const std::string s_OnlyCompareSimID;
    static const std::string s_UseYZVariance;
    static const std::string s_ParallelizeReps;
    static const std::string s_ProjectDescription;

	// Default constructor
	TController();

    // Create example input files in current working directory
    void GenerateExampleInputFiles();

    // Read input file (optional: select job, 0 = all)
    void ReadInputFile(const std::string& filename, const std::uint32_t job_id = 0);

	// Execute simulations
	void ExecuteSimulations();

    // Collect results of already finished simulations
    void CollectResults();

    // Write input file (based on m_ParamSets)
    void WriteInputFile(const std::filesystem::path& filename) const;

private:
    // Write controller parameters
    void WriteHeader(std::ostream& o_str) const;

    // Write output file for single simulation
    void WriteOutputFile(const std::filesystem::path& filename, const TParamSet& params, const TResult& result) const;

    // Write summary output file for all simulations (based on m_ParamSets and m_Results)
    void WriteOutputFileSummary(const std::filesystem::path& filename, const std::string& incomplete_job_ids) const;

    // Write output file with mean of repetitions (based on m_ParamSets and m_Results)
    void WriteOutputFileMean(const std::filesystem::path& filename, const std::string& incomplete_job_ids) const;

    // Ready-switch: true = has valid input
    bool m_IsReady;

    // Finished-switch: true = has completed all simulations
    bool m_IsFinished;

    // Input file path
    std::filesystem::path m_InputFile;

    // Selected parameter set (0 = all)
    std::uint32_t m_SelectedSimID;

    // Selected repetition (0 = all)
    std::uint32_t m_SelectedRepID;

    // Identification number of this project
    std::uint32_t m_ProjectID;

    // Project name
    std::string m_ProjectName;

    // DOS file path
    std::filesystem::path m_DOSFile;

    // Output file path
    std::filesystem::path m_OutputFile;

    // Verbosity
    Verbosity m_VL;

    // Switch: 
    // true = adjust DOS energy axis to temperature
    // false = zero of DOS energy axis remains as in DOS file
	bool m_EFTAdjust;

    // Switch: 
    // true = distribute electrons initially according to Fermi-Dirac function
    // false = distribute electrons initially according to step function (lowest state energy first)
	bool m_InitialFDDistrib;

	// Switch: 
    // true = fit occupation for effective temperature and effective chemical potential
    // false = use input temperature and chemical potential for effective carriers
	bool m_TeffFit;

    // Switch: 
    // true = enforce that the number of electrons is exactly the number predicted from the F-D integral of random states
    // false = occupy all states through random number < F-D probability (leads to fluctuation of electron count)
	bool m_EnforceECount;

    // Switch: 
    // true = automatic adjustment of path cutoffs (spatial and energetic) after equilibration (used ranges + Constant::autocutoff_inc %)
    // false = no automatic adjustment of path cutoffs
    bool m_CutoffAutoAdjust;

    // Percentage for auto-adjust of distance cutoff (optional; 0.0 if not set)
    double m_DistCutoffAdjustPercentage;

    // Percentage for auto-adjust of energy difference cutoff (optional; 0.0 if not set)
    double m_EdiffCutoffAdjustPercentage;

    // Switch: 
    // true = compare parameter sets only based on SimID (ignoring differences in other parameters)
    // false = compare all parameters
    bool m_OnlyCompareSimID;

    // Switch: 
    // true = calculate Dyz from variances instead of mean squared displacements (only for non-zero field)
    // false = calculate Dyz from mean squared displacements
    bool m_UseYZVariance;

    // Switch: 
    // true = run only a single repetition (job <id> selects from serialized [SimID][RepID] list)
    // false = run all repetitions of a certain SimID (job <id> selects from [SimID] list)
    bool m_ParallelizeReps;

    // Description
    std::string m_ProjectDescription;

    // DOS
    std::unique_ptr<TDOS> m_DOS;

    // Parameters
    std::vector<std::unique_ptr<TParamSet>> m_ParamSets;

    // Results
    std::vector<std::vector<std::unique_ptr<const TResult>>> m_Results;
};

} // MC namespace
#endif  // TController_H_