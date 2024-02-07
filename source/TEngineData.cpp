#include "TEngineData.hpp"

#include <iostream>
#include <stdexcept>
#include <chrono>
#include <string>

#include "CustomExceptions.hpp"

// Default constructor
MC::TEngineData::TEngineData()
    : m_VL(Verbosity::MAXIMUM), m_ParametersReady(false), m_ThermalEnergy(0.0), m_GradPhiEnergy(0.0), 
    m_InvLocRadius(0.0), m_RndGen(0), m_TotalTime(0.0), m_SimulationReady(false), m_ResultsReady(false)
{

}

// Specify output verbosity
void MC::TEngineData::SetVerbosity(Verbosity vl)
{
    m_VL = vl;
}

// Specify density of states
void MC::TEngineData::SetDOS(std::unique_ptr<TDOS> dos)
{
    if (!dos)
    {
        m_DOS.reset(nullptr);
        m_ParametersReady = false;
        m_SimulationReady = false;
        m_ResultsReady = false;
        return;
    }
    if (dos->HasDOS() == false) throw EX::TInvalidStatus("DOS is not properly defined.",__func__);

    m_DOS = std::move(dos);
    m_ParametersReady = false;
    m_SimulationReady = false;
    m_ResultsReady = false;
}

// Validate parameters
void MC::TEngineData::ValidateParameters(const TParamSet& params)
{
    if ((params.m_Repetitions < 1U) || (params.m_Repetitions > 100U))
    {
        throw EX::TInvalidInput("Number of repetitions out of range (1 - 100).");
    }
	if ((params.mI_StateCount < 100U) || (params.mI_StateCount > 10000000U))
    {
        throw EX::TInvalidInput("Number of states out of range (200 - 20000000).");
    }
    // Formula for max paths: The maximum cut-off distance is half of the cell size because otherwise
    // the periodic boundary conditions would allow two different paths to the same target state. The
    // respective sphere with radius d/2 (d = cell size) and volume Vsphere = 4/3*pi*(d/2)^3 has the 
    // same state density as the cell with volume Vcell = d^3, such that the number of states in the
    // sphere is Nsphere = Ncell * Vsphere / Vcell 
    // (minus 1 yields the approximate maximum number of hop destinations)
    const std::uint32_t max_paths = static_cast<std::uint32_t>(params.mI_StateCount*Constant::pi/6.0) - 1;
	if ((params.m_DistCutoff == 0.0) &&
        ((params.m_MinPathCount < 2U) || (params.m_MinPathCount > max_paths)))
    {
        throw EX::TInvalidInput("Minimum number of paths out of range (2 - " + std::to_string(max_paths) + ").");
    }
    if (params.m_DistCutoff < 0.0)
    {
        throw EX::TInvalidInput("Distance cutoff is negative.");
    }
    if (params.m_EdiffCutoff < 0.0)
    {
        throw EX::TInvalidInput("Absolute energy difference cutoff is negative.");
    }
    if (params.mI_PreHopLimit > 1000000000000U)
    {
        throw EX::TInvalidInput("Pre-equilibration hops out of range (0 - 2000000000000).");
    }
    if (params.mI_EqHopLimit > 1000000000000U)
    {
        throw EX::TInvalidInput("Equilibration hops out of range (0 - 2000000000000).");
    }
    if ((params.mI_HopLimit < 100U) || (params.mI_HopLimit > 1000000000000U))
    {
        throw EX::TInvalidInput("Hop limit out of range (200 - 2000000000000).");
    }
    if (params.m_AttemptTime <= 0.0)
    {
        throw EX::TInvalidInput("Attempt time zero or negative.");
    }
    if (params.m_AttemptTime > 1.0E100)
    {
        throw EX::TInvalidInput("Attempt time unrealistic large (> 1.0E100 s).");
    }
    if (params.m_Temperature <= 0.0)
    {
        throw EX::TInvalidInput("Temperature zero or negative.");
    }
    if (params.m_Temperature < 0.01)
    {
        throw EX::TInvalidInput("Temperature unrealistic low (< 0.01 K).");
    }
    if (params.m_PhiGradient < 0.0)
    {
        throw EX::TInvalidInput("Electric field strength is negative.");
    }
    if (params.m_LocRadius <= 0.0)
    {
        throw EX::TInvalidInput("Localization radius zero or negative.");
    }
}

// Specify parameters
void MC::TEngineData::SetParameters(std::unique_ptr<const TParamSet> params, std::uint32_t rep_id, std::uint64_t seed_increment)
{
    if (!m_DOS) throw EX::TInvalidStatus("DOS is not available.",__func__);
    
    if (!params)
    {
        if (m_VL >= Verbosity::MEDIUM) std::cout << "Clearing parameters." << std::endl;
        m_ParamSet.reset(nullptr);
        m_SimResult.reset(nullptr);
        m_DOS->DeleteParameters();
        m_ParametersReady = false;
        m_SimulationReady = false;
        m_ResultsReady = false;
        return;
    }

    if (m_DOS->HasDOS() == false) throw EX::TInvalidStatus("DOS is not properly defined.",__func__);

    if (m_VL >= Verbosity::MEDIUM) std::cout << "Applying parameters: " << std::flush;
    
    // Validate parameters
    ValidateParameters(*params);

    // Apply parameters to DOS
    std::string dos_output = "";
    m_DOS->ApplyParameters(*params,dos_output);

    // Validate cellsize-related parameters
    if (params->m_DistCutoff >= 0.5 * m_DOS->GetSpatialFactor())
    {
        throw EX::TInvalidInput("Distance cutoff exceeds half cell size (self-interaction).");
    }

    // Transfer parameters
    m_ParamSet = std::move(params);
    m_ThermalEnergy = Constant::kboltz*(m_ParamSet->m_Temperature)/(m_DOS->GetEnergyFactor());
    m_GradPhiEnergy = m_ParamSet->m_PhiGradient*1.0E-7*(m_DOS->GetSpatialFactor())/(m_DOS->GetEnergyFactor());
    m_InvLocRadius = m_DOS->GetSpatialFactor()/(m_ParamSet->m_LocRadius);
    std::string seed_str;
    if (m_ParamSet->m_RndSeed > 0)
    {
        seed_str = std::to_string(m_ParamSet->m_RndSeed) + " + " + std::to_string(seed_increment);
        m_RndGen.seed(static_cast<std::uint64_t>(m_ParamSet->m_RndSeed) + seed_increment);
    }
    else
    {
        auto clock_seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        seed_str = std::to_string(clock_seed) + " + " + std::to_string(-m_ParamSet->m_RndSeed) 
            + " + " + std::to_string(seed_increment);
        m_RndGen.seed(clock_seed + static_cast<std::uint64_t>(-m_ParamSet->m_RndSeed) + seed_increment);
    }

    // Reset structure and electrons
    m_Structure = std::vector<TLocalState>();
    m_Electrons = std::vector<TElectron>();

    // Create results object
    try
    {
        m_SimResult = std::make_unique<TResult>();
        m_SimResult->m_RepID = rep_id;
        m_SimResult->m_CellSize = m_DOS->GetSpatialFactor();
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot create results container.",__func__,e.what());
    }

    // Print parameter-related overview
    if (m_VL >= Verbosity::MEDIUM) 
    {
        std::cout << "done" << std::endl;
        if (dos_output != "") std::cout << "  " << dos_output << std::endl;
        std::cout << "  RNG-Seed: " << seed_str << std::endl;
        std::cout << "  Temperature: " << m_ParamSet->m_Temperature << " K, Localization radius: " 
            << m_ParamSet->m_LocRadius << " nm" << std::endl;
        std::cout << "  State energy range: " << m_ParamSet->m_MinStateEnergy << " eV to " << m_ParamSet->m_MaxStateEnergy 
            << " eV, Chemical potential: " << m_ParamSet->m_ChemPot << " eV" << std::endl;
        std::cout << "  Electric field strength: " << m_ParamSet->m_PhiGradient << " V/cm, Voltage drop across cell: " 
            << m_GradPhiEnergy * m_DOS->GetEnergyFactor() << " V" << std::endl;
        std::cout << "  Inverse hop attempt frequency: " << m_ParamSet->m_AttemptTime << " s" << std::endl;
        std::cout << "  Pre-equilibration hops: " << m_ParamSet->mO_PreHopLimit() 
            << ", Equilibration hops: " << m_ParamSet->mO_EqHopLimit() 
            << ", Simulation hops: " << m_ParamSet->mO_HopLimit() << std::endl;
    }

    m_ParametersReady = true;
    m_SimulationReady = false;
    m_ResultsReady = false;
}

// Return results object (removes histograms from this object, generate again if needed)
void MC::TEngineData::GetResults(std::unique_ptr<const TResult>& result)
{
    if (m_ResultsReady)
    {
        // Leave copy of non-histogram data (in case e.g. new structure is generated without new parameters)
        std::unique_ptr<TResult> temp = TResult::ValueCopy(m_SimResult.get());

        result = std::move(m_SimResult);
        m_SimResult = std::move(temp);
    }
    else
    {
        result = nullptr;
        throw EX::TInvalidStatus("No results are available.",__func__);
    }
    m_ResultsReady = false;
}

// Print diagnostics data (e.g. for debug)
void MC::TEngineData::PrintDiagnostics() const
{
    if (!m_DOS) 
        throw EX::TInvalidStatus("DOS is not available.",__func__); 
    if (m_DOS->IsReady() == false) 
        throw EX::TInvalidStatus("DOS is not ready.",__func__);   
    if (!m_ParamSet) 
        throw EX::TInvalidStatus("Parameters are not available.",__func__);   
    if (m_ParametersReady == false) 
        throw EX::TInvalidStatus("Parameters are not ready.",__func__);
    if ((m_Structure.empty()) || (m_Electrons.empty())) 
        throw EX::TInvalidStatus("Structure or electrons are not defined.",__func__);

    std::cout << " --- DIAGNOSTICS ---" << std::endl;

    std::cout << "Energy scaling factor: " << m_DOS->GetEnergyFactor() << " eV" << std::endl;
    std::cout << "Spatial scaling factor: " << m_DOS->GetSpatialFactor() << " nm" << std::endl;
    std::cout << "Thermal energy (kBT): " << m_ThermalEnergy << " <=> " 
        << m_ThermalEnergy * m_DOS->GetEnergyFactor() << " eV" << std::endl;
    	std::cout << "Electric potential energy gradient in +x-direction: " << m_GradPhiEnergy << " <=> " 
        << m_GradPhiEnergy * m_DOS->GetEnergyFactor() / m_DOS->GetSpatialFactor() << " eV/nm" << std::endl;
	std::cout << "Inverse localization radius: " << m_InvLocRadius << " <=> " 
        << m_InvLocRadius / m_DOS->GetSpatialFactor() << " 1/nm" << std::endl;

    std::cout << "Example state (first): " << m_Structure.front().Write() << " <=> " << 
        m_Structure.front().Write(m_DOS->GetSpatialFactor(), m_ParamSet->m_MaxStateEnergy, 
        m_DOS->GetEnergyFactor()) << std::endl;
    std::cout << "Example state (last): " << m_Structure.back().Write() << " <=> " << 
        m_Structure.back().Write(m_DOS->GetSpatialFactor(), m_ParamSet->m_MaxStateEnergy, 
        m_DOS->GetEnergyFactor()) << std::endl;

    std::cout << " --- DIAGNOSTICS FINISHED ---" << std::endl;
}