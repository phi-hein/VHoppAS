#include "TEngine.hpp"

#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <memory>
#include <cmath>
#include <array>
#include <chrono>

#include "TPath.hpp"
#include "Constants.hpp"
#include "GlobalFunctions.hpp"
#include "CustomExceptions.hpp"

// Default constructor
MC::TEngine::TEngine()
    : TEngineData()
{

}

// Fill the structure with local states
void MC::TEngine::GenerateStructure()
{
    if (m_VL >= Verbosity::MEDIUM) std::cout << "Generating random structure: " << std::flush;
    if (!m_DOS)
        throw EX::TInvalidStatus("DOS is not available.",__func__);
    if (!m_ParamSet)
        throw EX::TInvalidStatus("Parameters are not available.",__func__);
    if ((m_ParametersReady == false) || (m_DOS->IsReady() == false))
        throw EX::TInvalidStatus("Parameters are not properly defined.",__func__);

    m_Structure = std::vector<TLocalState>();
    try
    {
        m_Structure.reserve(m_ParamSet->m_StateCount);
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot create localized states.",__func__,e.what());
    }
    
    std::uniform_real_distribution<double> position_distribution (0.0,1.0);
    while (m_Structure.size() < m_ParamSet->m_StateCount)
    {
        // Random position in a unit cube
        double pos_x = position_distribution(m_RndGen);
        double pos_y = position_distribution(m_RndGen);
        double pos_z = position_distribution(m_RndGen);

        // Random state energy from DOS
        double energy = m_DOS->GetRandomEnergy(m_RndGen);

        m_Structure.emplace_back(pos_x,pos_y,pos_z,energy,m_ParamSet->m_StateCount);
    }

    // Reset electrons
    m_Electrons.clear();

    // Print structure-related overview
    if (m_VL >= Verbosity::MEDIUM) 
    {
        std::cout << "done" << std::endl;
        std::cout << "  Number of localized states: " << m_ParamSet->m_StateCount << std::endl;
        std::cout << "  Cell size: " << m_DOS->GetSpatialFactor() << " nm, Cubic cell volume: " 
            << 1.0E-21*pow(m_DOS->GetSpatialFactor(),3.0) << " cm3" << std::endl;
    }

    m_SimulationReady = false;
    m_ResultsReady = false;
}

// Find hopping paths between local states and calculate hopping times
void MC::TEngine::GeneratePaths()
{
    if (!m_DOS)
        throw EX::TInvalidStatus("DOS is not available.",__func__);
    if (!m_ParamSet)
        throw EX::TInvalidStatus("Parameters are not available.",__func__);
    if ((m_ParametersReady == false) || (m_DOS->IsReady() == false))
        throw EX::TInvalidStatus("Parameters are not properly defined.",__func__);
    if (m_Structure.empty())
        throw EX::TInvalidStatus("Structure is not defined.",__func__);
    if (!m_SimResult)
        throw EX::TInvalidStatus("Result object not available.",__func__);

    if (m_VL >= Verbosity::MEDIUM)
    {
        if (m_ParamSet->m_DistCutoff != 0.0)
            std::cout << "Validate cut-offs: " << std::flush;
        else
            std::cout << "Calculate cut-offs: " << std::flush;
    }

    // Determine cut-offs
    // - if distance cut-off is not defined:
    //     -> find required distance cut-off (in relative units: 0 to 1) for the minimal number of paths
    // - if distance cut-off is defined:
    //     -> find required distance cut-off (in relative units: 0 to 1) for at least two paths per state
    // - find required ediff cut-off (in relative units: 0 to 1) for at least two paths per state
    //   (ediff cut-off evaluated for paths that are within the distance cut-off)
    double dist_cutoff = 0.0;
    bool dist_cutoff_too_low = false;
    std::vector<double> min_dists;  
    if (m_ParamSet->m_DistCutoff != 0.0)
        min_dists = std::vector<double>(2, 0.0);
    else
        min_dists = std::vector<double>(m_ParamSet->m_MinPathCount, 0.0);    
    for (auto state = m_Structure.begin(); state < m_Structure.end(); ++state)
    {
        const double i_x = state->m_Pos_x;
        const double i_y = state->m_Pos_y;
        const double i_z = state->m_Pos_z;
        auto dist_end_it = min_dists.begin();
        auto dist_last_it = min_dists.begin();
        for (auto nn_state = m_Structure.begin(); nn_state < m_Structure.end(); ++nn_state)
        {
            if (state == nn_state) continue;

            double nn_x = nn_state->m_Pos_x;
            if (i_x - nn_x > 0.5) nn_x += 1.0;
            if (nn_x - i_x > 0.5) nn_x -= 1.0;
            double nn_y = nn_state->m_Pos_y;
            if (i_y - nn_y > 0.5) nn_y += 1.0;
            if (nn_y - i_y > 0.5) nn_y -= 1.0;
            double nn_z = nn_state->m_Pos_z;
            if (i_z - nn_z > 0.5) nn_z += 1.0;
            if (nn_z - i_z > 0.5) nn_z -= 1.0;
            const double nn_dist = sqrt((i_x-nn_x)*(i_x-nn_x) + (i_y-nn_y)*(i_y-nn_y) + (i_z-nn_z)*(i_z-nn_z));

            if (dist_end_it == min_dists.begin())
            {
                ++dist_end_it;
                *dist_last_it = nn_dist;
                continue;
            }
            if (dist_end_it != min_dists.end())
            {
                ++dist_end_it;
                ++dist_last_it;
                *dist_last_it = nn_dist;
            }
            else if (nn_dist < *dist_last_it)
            {
                *dist_last_it = nn_dist;
            }
            else
            {
                continue;
            }
            auto dist_sort_it = std::upper_bound(min_dists.begin(),dist_last_it,nn_dist);
            if (dist_sort_it != dist_last_it)
            {
                std::rotate(dist_sort_it,dist_last_it,dist_end_it);
            }
        }
        if (*dist_last_it > dist_cutoff)
        {
            dist_cutoff = *dist_last_it;
        }
    }
    dist_cutoff += 0.001*dist_cutoff;

    if (m_ParamSet->m_DistCutoff != 0.0)
    {
        if (m_ParamSet->m_DistCutoff / m_DOS->GetSpatialFactor() > dist_cutoff)
            dist_cutoff = m_ParamSet->m_DistCutoff / m_DOS->GetSpatialFactor();
        else
            dist_cutoff_too_low = true;
    }
    if (dist_cutoff >= 0.5)
    {
        throw EX::TInvalidInput("Cut-off radius exceeds half cell (self-interaction).");
    }

    double ediff_cutoff = 0.0;
    bool ediff_cutoff_too_low = false;
    if (m_ParamSet->m_EdiffCutoff != 0.0)
    {
        std::vector<double> min_ediffs (2, 0.0);
        for (auto state = m_Structure.begin(); state < m_Structure.end(); ++state)
        {
            const double i_x = state->m_Pos_x;
            const double i_y = state->m_Pos_y;
            const double i_z = state->m_Pos_z;
            auto ediff_end_it = min_ediffs.begin();
            auto ediff_last_it = min_ediffs.begin();
            for (auto nn_state = m_Structure.begin(); nn_state < m_Structure.end(); ++nn_state)
            {
                if (state == nn_state) continue;

                double nn_x = nn_state->m_Pos_x;
                if (i_x - nn_x > 0.5) nn_x += 1.0;
                if (nn_x - i_x > 0.5) nn_x -= 1.0;
                double nn_y = nn_state->m_Pos_y;
                if (i_y - nn_y > 0.5) nn_y += 1.0;
                if (nn_y - i_y > 0.5) nn_y -= 1.0;
                double nn_z = nn_state->m_Pos_z;
                if (i_z - nn_z > 0.5) nn_z += 1.0;
                if (nn_z - i_z > 0.5) nn_z -= 1.0;
                const double nn_dist = sqrt((i_x-nn_x)*(i_x-nn_x) + (i_y-nn_y)*(i_y-nn_y) + (i_z-nn_z)*(i_z-nn_z));

                if (nn_dist > dist_cutoff) continue;

                const double nn_ediff = fabs(state->m_Energy - nn_state->m_Energy);

                if (ediff_end_it == min_ediffs.begin())
                {
                    ++ediff_end_it;
                    *ediff_last_it = nn_ediff;
                    continue;
                }
                if (ediff_end_it != min_ediffs.end())
                {
                    ++ediff_end_it;
                    ++ediff_last_it;
                    *ediff_last_it = nn_ediff;
                }
                else if (nn_ediff < *ediff_last_it)
                {
                    *ediff_last_it = nn_ediff;
                }
                else
                {
                    continue;
                }
                auto ediff_sort_it = std::upper_bound(min_ediffs.begin(),ediff_last_it,nn_ediff);
                if (ediff_sort_it != ediff_last_it)
                {
                    std::rotate(ediff_sort_it,ediff_last_it,ediff_end_it);
                }
            }
            if (*ediff_last_it > ediff_cutoff)
            {
                ediff_cutoff = *ediff_last_it;
            }
        }
        ediff_cutoff += 0.001*ediff_cutoff;

        if (m_ParamSet->m_EdiffCutoff / m_DOS->GetEnergyFactor() > ediff_cutoff)
            ediff_cutoff = m_ParamSet->m_EdiffCutoff / m_DOS->GetEnergyFactor();
        else
            ediff_cutoff_too_low = true;
        if (ediff_cutoff > 1.0) ediff_cutoff = 1.0;
    }
    else ediff_cutoff = 1.0;

    m_SimResult->m_MaxPathDist = dist_cutoff * m_DOS->GetSpatialFactor();
    m_SimResult->m_MaxPathEdiff = ediff_cutoff * m_DOS->GetEnergyFactor();

    // Print cutoff-related overview
    if (m_VL >= Verbosity::MEDIUM) 
    {
        std::cout << "done" << std::endl;
        if (dist_cutoff_too_low)
            std::cout << "  ! Distance cut-off too low -> adjusted to ensure at least two paths per state." << std::endl;
        if (ediff_cutoff_too_low)
            std::cout << "  ! Energy difference cut-off too low -> adjusted to ensure at least two paths per state." << std::endl;
        std::cout << "  Path cut-off for distance: " << m_SimResult->m_MaxPathDist << " nm" << std::endl;
        std::cout << "  Path cut-off for state energy difference: " << m_SimResult->m_MaxPathEdiff << " eV" << std::endl;
        std::cout << "Generating hopping paths: " << std::flush;
    }

    // Delete all paths
    for (auto& state : m_Structure)
    {
        state.m_Paths = std::vector<TPath>();
    }

    // For each localized state find the states within the cut-offs
    std::uint64_t path_count = 0;
    GF::TMinMaxMean<std::uint32_t> paths_per_state;
    for (std::uint32_t i = 0; i < m_Structure.size(); ++i)
    {
        const double i_x = m_Structure[i].m_Pos_x;
        const double i_y = m_Structure[i].m_Pos_y;
        const double i_z = m_Structure[i].m_Pos_z;
        std::uint32_t k = i + 1U;
        while (k < m_Structure.size())
        {
            // Shortest neighbor position with periodic boundaries
            double nn_x = m_Structure[k].m_Pos_x;
            if (i_x - nn_x > 0.5) nn_x += 1.0;
            if (nn_x - i_x > 0.5) nn_x -= 1.0;
            double nn_y = m_Structure[k].m_Pos_y;
            if (i_y - nn_y > 0.5) nn_y += 1.0;
            if (nn_y - i_y > 0.5) nn_y -= 1.0;
            double nn_z = m_Structure[k].m_Pos_z;
            if (i_z - nn_z > 0.5) nn_z += 1.0;
            if (nn_z - i_z > 0.5) nn_z -= 1.0;
            double nn_dist = sqrt((i_x-nn_x)*(i_x-nn_x) + (i_y-nn_y)*(i_y-nn_y) + (i_z-nn_z)*(i_z-nn_z));
            double nn_ediff = fabs(m_Structure[i].m_Energy - m_Structure[k].m_Energy);

            if ((nn_dist <= dist_cutoff) && (nn_ediff <= ediff_cutoff))
            {
                // Add pair of paths (including reverse_id)
                std::uint32_t i_size = static_cast<std::uint32_t>(m_Structure[i].m_Paths.size());
                std::uint32_t k_size = static_cast<std::uint32_t>(m_Structure[k].m_Paths.size());
                try
                {
                    m_Structure[i].m_Paths.emplace_back(k,k_size,0.0);
                    m_Structure[k].m_Paths.emplace_back(i,i_size,0.0);
                }
                catch(const std::bad_alloc& e)
                {
                    throw EX::TOutOfMemory("Cannot create all paths.",__func__,e.what());
                }
                path_count += 2;
            }

            ++k;
        }
        paths_per_state.check(static_cast<std::uint32_t>(m_Structure[i].m_Paths.size()));
    }

    if (paths_per_state.min < 2) 
        throw EX::TInvalidStatus("States with less than two paths detected.",__func__);
    
    m_SimResult->m_MaxPathCount = paths_per_state.max;
    m_SimResult->m_MeanPathCount = paths_per_state.mean;

    // Print path-related overview
    if (m_VL >= Verbosity::MEDIUM) 
    {
        std::cout << "done" << std::endl;
        if (m_ParamSet->m_DistCutoff == 0.0)
            std::cout << "  Required number of paths per state: " << m_ParamSet->m_MinPathCount << std::endl;
        std::cout << "  Average number of paths per state: " << paths_per_state << std::endl;
        std::cout << "  Total number of paths: " << path_count << std::endl;
    }

    m_SimulationReady = false;
    m_ResultsReady = false;

    // Calculate hopping times for each path
    UpdateHoppingTimes();
}

// Recalculate hopping times for all hopping paths (e.g. after changes of electric field or temperature)
void MC::TEngine::UpdateHoppingTimes()
{
    if (m_VL >= Verbosity::MEDIUM) std::cout << "Calculate hopping times: " << std::flush;
    if (!m_DOS)
        throw EX::TInvalidStatus("DOS is not available.",__func__);
    if (!m_ParamSet)
        throw EX::TInvalidStatus("Parameters are not available.",__func__);
    if ((m_ParametersReady == false) || (m_DOS->IsReady() == false))
        throw EX::TInvalidStatus("Parameters are not properly defined.",__func__);
    if (m_Structure.empty())
        throw EX::TInvalidStatus("Structure is not defined.",__func__);
    for (auto& state : m_Structure)
    {
        if (state.m_Paths.empty())
            throw EX::TInvalidStatus("Paths are not properly defined.",__func__);
    }

    GF::TMinMaxMean<double> avg_nn_dist;
    GF::TMinMaxMean<double> avg_nn_ediff;
    GF::TMinMaxMean<double> avg_nn_ediff_field;
    for (auto& state : m_Structure)
    {
        double state_nn_dist = 0.0;
        double state_nn_ediff = 0.0;
        double state_nn_ediff_field = 0.0;
        for (std::size_t j = 0; j < state.m_Paths.size(); ++j)
        {
            TPath& path = state.m_Paths[j];

            // Calculate path length
            double nn_x = m_Structure[path.m_StateID].m_Pos_x;
            if (state.m_Pos_x - nn_x > 0.5) nn_x += 1.0;
            if (nn_x - state.m_Pos_x > 0.5) nn_x -= 1.0;
            double nn_y = m_Structure[path.m_StateID].m_Pos_y;
            if (state.m_Pos_y - nn_y > 0.5) nn_y += 1.0;
            if (nn_y - state.m_Pos_y > 0.5) nn_y -= 1.0;
            double nn_z = m_Structure[path.m_StateID].m_Pos_z;
            if (state.m_Pos_z - nn_z > 0.5) nn_z += 1.0;
            if (nn_z - state.m_Pos_z > 0.5) nn_z -= 1.0;
            double nn_dist = sqrt((state.m_Pos_x - nn_x)*(state.m_Pos_x - nn_x) + 
                (state.m_Pos_y - nn_y)*(state.m_Pos_y - nn_y) + (state.m_Pos_z - nn_z)*(state.m_Pos_z - nn_z));

            // Calculate time (inverse rate)
            double energy_difference = m_Structure[path.m_StateID].m_Energy - state.m_Energy
                - (nn_x - state.m_Pos_x)*m_GradPhiEnergy;
            if (energy_difference < 0.0)
            {
                path.m_Time = m_ParamSet->m_AttemptTime*exp(2.0*m_InvLocRadius*nn_dist);
            }
            else
            {
                path.m_Time = m_ParamSet->m_AttemptTime*exp(2.0*m_InvLocRadius*nn_dist + 
                    energy_difference/m_ThermalEnergy);
            }

            // Distance and energy statistics (restricted to nearest-neighbors)
            if ((j == 0) || (nn_dist < state_nn_dist))
            {
                state_nn_dist = nn_dist;
                state_nn_ediff = fabs(m_Structure[path.m_StateID].m_Energy - state.m_Energy);
                state_nn_ediff_field = fabs(energy_difference);
            }
        }
        
        avg_nn_dist.check(state_nn_dist * m_DOS->GetSpatialFactor());
        avg_nn_ediff.check(state_nn_ediff * m_DOS->GetEnergyFactor());
        avg_nn_ediff_field.check(state_nn_ediff_field * m_DOS->GetEnergyFactor());
    }

    // Print energy- and distance-related overview
    if (m_VL >= Verbosity::MEDIUM) 
    {
        std::cout << "done" << std::endl;
        std::cout << "  Average distance betw. nearest-neighbors: " << avg_nn_dist << " nm" << std::endl;
        std::cout << "  Average absolute energy difference betw. nearest-neigbors: " << avg_nn_ediff << " eV" << std::endl;
        if (m_ParamSet->m_PhiGradient != 0.0)
        {
            std::cout << "  Average abs. energy diff. betw. nearest-neigbors incl. field: " << avg_nn_ediff_field << " eV" << std::endl;
        }
    }

    m_SimulationReady = false;
    m_ResultsReady = false;
}

// Fill the structure with electrons according to F-D distribution
void MC::TEngine::GenerateElectrons()
{
    if (m_VL >= Verbosity::MEDIUM) std::cout << "Generate electrons: " << std::flush;
    if (!m_DOS)
        throw EX::TInvalidStatus("DOS is not available.",__func__);
    if (!m_ParamSet)
        throw EX::TInvalidStatus("Parameters are not available.",__func__);
    if ((m_ParametersReady == false) || (m_DOS->IsReady() == false))
        throw EX::TInvalidStatus("Parameters are not properly defined.",__func__);
    if (m_Structure.empty())
        throw EX::TInvalidStatus("Structure is not defined.",__func__);
    if (!m_SimResult)
        throw EX::TInvalidStatus("Result object not available.",__func__);

    // Calculate expected number of electrons based on random structure
    double expected_esum = 0.0;
    for (const auto& state : m_Structure)
    {
        expected_esum += m_DOS->GetOccupationProbability(state.m_Energy,m_ThermalEnergy);
    }
    const std::uint32_t expected_ecount = static_cast<std::uint32_t>(expected_esum + 0.5);

    // Reserve memory for electrons
    m_Electrons = std::vector<TElectron>();
    try
    {
        m_Electrons.reserve(expected_ecount);
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot create all electrons.",__func__,e.what());
    }

    std::int32_t enforced_electrons = 0;
    if (m_ParamSet->m_InitialFDDistrib)
    {
        // Place electrons in structure based on F-D probability
        std::uniform_real_distribution<double> probability (0.0,1.0);
        for (std::uint32_t i = 0; i < m_Structure.size(); ++i)
        {
            if (probability(m_RndGen) < m_DOS->GetOccupationProbability(m_Structure[i].m_Energy,m_ThermalEnergy))
            {
                m_Structure[i].m_ElectronID = static_cast<std::uint32_t>(m_Electrons.size());
                m_Electrons.emplace_back(i);
            }
            else
            {
                m_Structure[i].m_ElectronID = m_ParamSet->m_StateCount;
            }
        }

        // Enforce expected number of electrons (if enabled)
        if ((m_ParamSet->m_EnforceECount) && (m_Electrons.size() != expected_ecount))
        {
            if (expected_ecount == 0)
            {
                throw EX::TInvalidInput("Enforcing expected electron count leads to zero electrons.");
            }
            if (expected_ecount >= m_Structure.size())
            {
                m_Electrons = std::vector<TElectron>();
                throw EX::TInvalidInput("Enforcing expected electron count leads to zero empty states.");
            }

            // Generate list of empty or occupied state indices
            std::vector<std::size_t> state_list;
            if (m_Electrons.size() < expected_ecount)
            {
                for (std::size_t i = 0; i < m_Structure.size(); ++i)
                    if (m_Structure[i].m_ElectronID == m_ParamSet->m_StateCount)
                        state_list.push_back(i);
            }
            else
            {
                for (std::size_t i = 0; i < m_Structure.size(); ++i)
                    if (m_Structure[i].m_ElectronID != m_ParamSet->m_StateCount)
                        state_list.push_back(i);
            }

            // Randomly re-roll state occupation to reach expected electron count
            while (m_Electrons.size() != expected_ecount)
            {
                std::size_t i = 0;
                if (state_list.size() > 1)
                {
                    std::uniform_int_distribution<std::size_t> list_entry (0, state_list.size() - 1);
                    i = list_entry(m_RndGen);
                }

                if (probability(m_RndGen) < m_DOS->GetOccupationProbability(m_Structure[state_list[i]].m_Energy,m_ThermalEnergy))
                {
                    if (m_Electrons.size() < expected_ecount)
                    {
                        m_Structure[state_list[i]].m_ElectronID = static_cast<std::uint32_t>(m_Electrons.size());
                        m_Electrons.emplace_back(state_list[i]);
                        state_list.erase(state_list.begin() + i);
                        ++enforced_electrons;
                    }
                }
                else
                {
                    if (m_Electrons.size() > expected_ecount)
                    {
                        m_Electrons.erase(m_Electrons.begin() + m_Structure[state_list[i]].m_ElectronID);
                        for (auto it = m_Electrons.begin() + m_Structure[state_list[i]].m_ElectronID; it != m_Electrons.end(); ++it)
                        {
                            --(m_Structure[it->m_CurrentStateID].m_ElectronID);
                        }
                        m_Structure[state_list[i]].m_ElectronID = m_ParamSet->m_StateCount;
                        state_list.erase(state_list.begin() + i);
                        --enforced_electrons;
                    }
                }
            }
        }
    }
    // Place electrons in structure based on step function (lowest state energies first)
    else
    {
        // Validate expected electron count
        if (expected_ecount == 0)
        {
            throw EX::TInvalidInput("Parameters lead to expected zero electrons.");
        }
        if (expected_ecount >= m_Structure.size())
        {
            m_Electrons = std::vector<TElectron>();
            throw EX::TInvalidInput("Parameters lead to expected zero empty states.");
        }

        // Sort state indices by ascending energy
        std::vector<std::size_t> state_sort (m_Structure.size(),0);
        std::iota(state_sort.begin(),state_sort.end(),0);
        std::sort(state_sort.begin(),state_sort.end(),
            [&](std::size_t i, std::size_t j)
            {
                return m_Structure[i].m_Energy < m_Structure[j].m_Energy;
            });
        
        // Place electrons in structure based on sorted states
        for (std::size_t i = 0; i < m_Structure.size(); ++i)
        {
            if (m_Electrons.size() < expected_ecount)
            {
                m_Structure[state_sort[i]].m_ElectronID = static_cast<std::uint32_t>(m_Electrons.size());
                m_Electrons.emplace_back(state_sort[i]);
            }
            else
            {
                m_Structure[state_sort[i]].m_ElectronID = m_ParamSet->m_StateCount;
            }
        }
    }

    // Finalize electron list (including consistency check for inter-linked indices)
    m_Electrons.shrink_to_fit();
    if (m_Electrons.empty())
    {
        throw EX::TInvalidInput("Parameters lead to zero electrons.");
    }
    if (m_Electrons.size() == m_Structure.size())
    {
        m_Electrons = std::vector<TElectron>();
        throw EX::TInvalidInput("Parameters lead to zero empty states.");
    }
    for (std::size_t i = 0; i < m_Structure.size(); ++i)
    {
        if (m_Structure[i].m_ElectronID == m_ParamSet->m_StateCount)
        {
            for (const auto& electron : m_Electrons)
                if (electron.m_CurrentStateID == i)
                    throw EX::TInvalidStatus("Inconsistency in state occupation.",__func__);
        }
        else
        {
            if ((m_Structure[i].m_ElectronID >= m_Electrons.size()) ||
                (m_Electrons[m_Structure[i].m_ElectronID].m_CurrentStateID != i))
                throw EX::TInvalidStatus("Inconsistency in state identification.",__func__);
        }
    }
    for (std::size_t i = 0; i < m_Electrons.size(); ++i)
    {
        if ((m_Electrons[i].m_CurrentStateID >= m_Structure.size()) ||
            (m_Structure[m_Electrons[i].m_CurrentStateID].m_ElectronID != i))
            throw EX::TInvalidStatus("Inconsistency in electron identification.",__func__);
    }
    m_SimResult->m_ElectronCount = static_cast<std::uint32_t>(m_Electrons.size());
    
    if (m_VL >= Verbosity::MEDIUM) 
    {
        // Calculate total energy
        double total_energy = 0.0;
        const double enmax = m_ParamSet->m_MaxStateEnergy;
        const double enfac = m_DOS->GetEnergyFactor();
        for (const auto &electron : m_Electrons)
        {
            total_energy += m_Structure[electron.m_CurrentStateID].m_Energy * enfac + enmax;
        }

        // Print electrons-related overview
        std::cout << "done" << std::endl;
        if (m_ParamSet->m_InitialFDDistrib)
            std::cout << "  Initial electron distribution: Fermi-Dirac function" << std::endl;
        else
            std::cout << "  Initial electron distribution: Step function" << std::endl;
        if ((m_ParamSet->m_EnforceECount) && (m_ParamSet->m_InitialFDDistrib))
            std::cout << "  Enforced adjustment of electron count: " << enforced_electrons << std::endl;
        std::cout << "  Number of electrons: " << m_Electrons.size() << " (expected from states: " 
            << expected_ecount << ", expected from DOS: " << m_DOS->GetElectronCount(m_ThermalEnergy) << ")" << std::endl;
        std::cout << "  Occupation ratio: " 
            << static_cast<double>(m_Electrons.size())/static_cast<double>(m_Structure.size()) << std::endl;
        std::cout << "  Initial total energy (sum of occupied state energies): " << total_energy << " eV" << std::endl;
    }

    m_SimulationReady = false;
    m_ResultsReady = false;
}

// Initialize the simulation data (e.g. randomized path times)
void MC::TEngine::InitializeSimulation()
{
    if (m_VL >= Verbosity::MEDIUM) std::cout << "Initializing simulation data: " << std::flush;
    if (!m_DOS)
        throw EX::TInvalidStatus("DOS is not available.",__func__);
    if (!m_ParamSet)
        throw EX::TInvalidStatus("Parameters are not available.",__func__);
    if ((m_ParametersReady == false) || (m_DOS->IsReady() == false))
        throw EX::TInvalidStatus("Parameters are not properly defined.",__func__);
    if (m_Structure.empty())
        throw EX::TInvalidStatus("Structure is not defined.",__func__);
    for (auto& state : m_Structure)
    {
        if (state.m_Paths.empty())
            throw EX::TInvalidStatus("Paths are not properly defined.",__func__);
        for (auto& path : state.m_Paths)
        {
            if (path.m_Time == 0.0)
                throw EX::TInvalidStatus("Paths times are not defined.",__func__);
        }
    }
    if (m_Electrons.empty())
        throw EX::TInvalidStatus("Electrons are not defined.",__func__);
    if (!m_SimResult)
        throw EX::TInvalidStatus("Result object not available.",__func__);

    // Assign randomized hopping times to all paths (reserve memory for max paths to prevent allocations during simulation)
    bool all_electrons_blocked = true;
    std::uniform_real_distribution<double> probability (0.0,1.0);
    for (auto &electron : m_Electrons)
    {
        electron.m_Disp_x = 0.0;
        electron.m_Disp_y = 0.0;
        electron.m_Disp_z = 0.0;
        electron.m_HopCount = 0U;
        electron.m_AccTimeDecay = 0.0;
        electron.m_LastStateID = electron.m_CurrentStateID;
        electron.m_OscHopCount = 0U;
        electron.m_FirstHopTime = 0.0;

        try
        {
            electron.m_RandomTimes.reserve(m_SimResult->m_MaxPathCount);
            electron.m_RandomTimes.resize(m_Structure[electron.m_CurrentStateID].m_Paths.size(),0.0);
        }
        catch(const std::bad_alloc& e)
        {
            throw EX::TOutOfMemory("Cannot create paths for electrons.",__func__,e.what());
        }

        for (std::uint32_t i = 0; i < electron.m_RandomTimes.size(); ++i)
        {
            if (m_Structure[m_Structure[electron.m_CurrentStateID].m_Paths[i].m_StateID].m_ElectronID 
                < m_ParamSet->m_StateCount)
            {
                electron.m_RandomTimes[i] = Constant::occupied_path;
            }
            else
            {
                electron.m_RandomTimes[i] = 
                    -log(probability(m_RndGen))*m_Structure[electron.m_CurrentStateID].m_Paths[i].m_Time;
                all_electrons_blocked = false;
            }
        }

        auto min_it = std::min_element(electron.m_RandomTimes.cbegin(),electron.m_RandomTimes.cend());
        electron.m_MinIndex = static_cast<std::uint32_t>(std::distance(electron.m_RandomTimes.cbegin(),min_it));
	    electron.m_MinTime = *min_it;
    }
    if (all_electrons_blocked)
    {
        throw EX::TInvalidInput("All paths are blocked for all electrons (try different seed).");
    }

    // Reset hop times and path hopping counters
    for (auto& state : m_Structure)
    {
        state.m_LastHopTime = 0.0;
        state.m_OccTime = 0.0;
        for (auto& path : state.m_Paths)
        {
            path.m_HopCount = 0U;
            path.m_OscHopCount = 0U;
        }
    }

    // Reset simulated timespan
    m_TotalTime = 0.0;

    // Reset progress
    m_SimResult->m_EqProgress.clear();
    m_SimResult->m_Progress.clear();

    // Calculate effective carriers and partial entropy from temperature and chemical potential
    // (already final values if TeffFit parameter is disabled)
    CalculateEffectiveCarrierDensity(false);

    m_SimulationReady = false;
    m_ResultsReady = false;
    if (m_VL >= Verbosity::MEDIUM) 
    {
        std::cout << "done" << std::endl;
        double dos_eff_density = m_DOS->GetEffCarrierDensity(m_ParamSet->m_ChemPot,m_ParamSet->m_Temperature);
        std::cout << "  Effective charge carriers" << ((m_ParamSet->m_TeffFit) ? " (preliminary): " : ": ")
            << m_SimResult->m_EffCarriers << " (expected from DOS: " 
            << dos_eff_density*1.0E-21*pow(m_DOS->GetSpatialFactor(),3.0) << ")" << std::endl;
        std::cout << "  Effective charge carrier density" << ((m_ParamSet->m_TeffFit) ? " (preliminary): " : ": ")
            << m_SimResult->m_EffCarrierDensity << " 1/cm3 (expected from DOS: " 
            << dos_eff_density << " 1/cm3)" << std::endl;
        std::cout << "  Partial entropy" << ((m_ParamSet->m_TeffFit) ? " (preliminary): " : ": ")
            << m_SimResult->m_PartialEntropy << " eV/K (expected from DOS: " 
            << m_DOS->GetPartialEntropy(m_ParamSet->m_ChemPot,m_ParamSet->m_Temperature) << " eV/K)" << std::endl;
    }
}

// Execute the simulation
void MC::TEngine::RunSimulation ()
{
    // Initialize data (also validates that required data is available)
    InitializeSimulation();

    // Hop counter and limits
    std::uint64_t hop_counter {0U};
    const std::uint64_t pre_hop_limit {m_ParamSet->m_PreHopLimit};
    const std::uint64_t eq_hop_limit {m_ParamSet->m_EqHopLimit};
    const std::uint64_t sim_hop_limit {m_ParamSet->m_HopLimit};
    const std::uint64_t total_hop_limit {pre_hop_limit + eq_hop_limit + sim_hop_limit};
    std::uint64_t next_hop_limit {pre_hop_limit};

    // Result of equilibration (nullptr until equilibration completed)
    std::unique_ptr<const TResult> eq_result;

    // Wall-clock time for start of current KMC task
    auto runtime_start = std::chrono::steady_clock::now();

    if ((pre_hop_limit != 0) && (m_VL >= Verbosity::MEDIUM))
    {
        std::cout << "Pre-equilibration: " << std::flush;
    }

    // Outer KMC-Loop
    while (hop_counter < total_hop_limit)
    {
        // Calculate intermediate results (incl. at end of equilibration for any verbosity)
        if ((hop_counter > pre_hop_limit) &&
            ((m_VL >= Verbosity::MEDIUM) || (hop_counter == pre_hop_limit + eq_hop_limit)))
        {
            if (m_ParamSet->m_TeffFit) CalculateEffectiveCarrierDensity(true);
            CalculateResultValues(eq_result.get());
        }

        // Output of intermediate results
        if (m_VL >= Verbosity::MEDIUM)
        {
            // Output: End of pre-equilibration
            if ((hop_counter == pre_hop_limit) && (pre_hop_limit != 0))
            {
                std::uint64_t mob_e_counter {0U};
                std::uint64_t osc_e_counter {0U};
                for (const TElectron& electron : m_Electrons)
                {
                    if (electron.m_HopCount != 0)
                    {
                        if (electron.m_HopCount != electron.m_OscHopCount) 
                            ++mob_e_counter;
                        else
                            ++osc_e_counter;
                    }
                }
                std::cout << "done (" << pre_hop_limit << " hops, " 
                    << mob_e_counter << " mobile electrons + " 
                    << osc_e_counter << " oscillating electrons)" << std::endl;
            }

            // Output: Progress of equilibration or simulation
            if (hop_counter > pre_hop_limit)
            {
                if (hop_counter <= pre_hop_limit + eq_hop_limit)
                {
                    m_SimResult->SaveProgress(
                        0.1*round(1000.0*static_cast<double>(hop_counter - pre_hop_limit)
                        /static_cast<double>(eq_hop_limit)),true);
                    m_SimResult->WriteProgressLine(std::cout,eq_hop_limit, 
                        (m_ParamSet->m_PhiGradient != 0.0),true);
                }
                else
                {
                    m_SimResult->SaveProgress(
                        0.1*round(1000.0*static_cast<double>(hop_counter - pre_hop_limit - eq_hop_limit)
                        /static_cast<double>(sim_hop_limit)),false);
                    m_SimResult->WriteProgressLine(std::cout,sim_hop_limit, 
                        (m_ParamSet->m_PhiGradient != 0.0),false);
                }
            }

            // Output: End of equilibration
            if ((hop_counter == pre_hop_limit + eq_hop_limit) && (eq_hop_limit != 0))
            {
                std::cout << " --- EQUILIBRATION END ---" << std::endl;
                std::cout << "(equilibration runtime: ";
                GF::WriteDuration(std::cout, std::chrono::steady_clock::now() - runtime_start);
                std::cout << ")" << std::endl;
            }
        }

        // Handle finished pre-equilibration
        // - electron distribution remains
        // - movement statistics is reset
        if ((hop_counter == pre_hop_limit) && (pre_hop_limit != 0))
        {
            m_TotalTime = 0.0;
            for (TLocalState& state : m_Structure)
            {
                state.m_LastHopTime = 0.0;
                state.m_OccTime = 0.0;
                for (TPath& path : state.m_Paths)
                {
                    path.m_HopCount = 0U;
                    path.m_OscHopCount = 0U;
                }
            }
            for (TElectron& electron : m_Electrons)
            {
                electron.m_HopCount = 0U;
                electron.m_Disp_x = 0.0;
                electron.m_Disp_y = 0.0;
                electron.m_Disp_z = 0.0;
                electron.m_LastStateID = electron.m_CurrentStateID;
                electron.m_OscHopCount = 0U;
                electron.m_FirstHopTime = 0.0;
            }
        }

        // Handle finished equilibration
        // - movement statistics remains
        // - final time and displacements are saved as starting point for simulation
        // - output of statistics (incl. optional auto-adjust of cut-offs)
        if ((hop_counter == pre_hop_limit + eq_hop_limit) && (eq_hop_limit != 0))
        {
            eq_result = TResult::ValueCopy(m_SimResult.get());
            CalculatePostEquilibrationStatistics();
        }

        // Output of headers for progress tables (start of equilibration or simulation)
        if ((m_VL >= Verbosity::MEDIUM) &&
            ((hop_counter == pre_hop_limit) || (hop_counter == pre_hop_limit + eq_hop_limit)))
        {
            runtime_start = std::chrono::steady_clock::now();
            
            if ((hop_counter == pre_hop_limit) && (eq_hop_limit != 0))
            {
                std::cout << " --- EQUILIBRATION START ---" << std::endl;
                m_SimResult->WriteProgressHeader(std::cout,eq_hop_limit,(m_ParamSet->m_PhiGradient != 0.0));
            }
            else
            {
                std::cout << " --- SIMULATION START ---" << std::endl;
                if (eq_hop_limit != 0) 
                    std::cout << " (time, hops and displacements relative to equilibration)" << std::endl;
                m_SimResult->WriteProgressHeader(std::cout,sim_hop_limit,(m_ParamSet->m_PhiGradient != 0.0));
            }
        }

        // Determine next hop limit (until next intermediate report)
        // -> Start of equilibration
        if ((hop_counter == pre_hop_limit) && (eq_hop_limit != 0))
        {
            next_hop_limit = pre_hop_limit + eq_hop_limit / 40U;
            if ((eq_hop_limit < 100) || (m_VL == Verbosity::MINIMUM))
            {
                next_hop_limit = pre_hop_limit + eq_hop_limit;
            }
        }
        // -> Progress of equilibration
        if ((hop_counter > pre_hop_limit) && (hop_counter < pre_hop_limit + eq_hop_limit))
        {
            next_hop_limit += eq_hop_limit / 40U;
            if (next_hop_limit > pre_hop_limit + eq_hop_limit - eq_hop_limit / 40U)
            {
                next_hop_limit = pre_hop_limit + eq_hop_limit;
            }
        }
        // -> Start of simulation
        if (hop_counter == pre_hop_limit + eq_hop_limit)
        {
            next_hop_limit = pre_hop_limit + eq_hop_limit + sim_hop_limit / 40U;
            if (m_VL == Verbosity::MINIMUM)
            {
                next_hop_limit = total_hop_limit;
            }
        }
        // -> Progress of simulation
        if (hop_counter > pre_hop_limit + eq_hop_limit)
        {
            next_hop_limit += sim_hop_limit / 40U;
            if (next_hop_limit > total_hop_limit - sim_hop_limit / 40U)
            {
                next_hop_limit = total_hop_limit;
            }
        }

        // INNER KMC LOOP
        // Execute hops until next limit is reached
        KMCLoop(hop_counter, next_hop_limit);
    } 

    // Add occupation timespan of states that are occupied at the end of simulation
    for (const auto& state: m_Structure)
    {
        if (state.m_ElectronID < m_ParamSet->m_StateCount)
            state.m_OccTime += m_TotalTime - state.m_LastHopTime;
    }
    
    // Calculate final results
    if (m_ParamSet->m_TeffFit) CalculateEffectiveCarrierDensity(true);
    std::unique_ptr<const TResult> incl_eq_result;
    if ((eq_hop_limit != 0) && (m_VL == Verbosity::MAXIMUM))
    {
        CalculateResultValues(nullptr);
        incl_eq_result = TResult::ValueCopy(m_SimResult.get());
    }
    CalculateResultValues(eq_result.get());
    m_SimulationReady = true;
    m_ResultsReady = false;

    // Output of last line of simulation progress
    if (m_VL >= Verbosity::MEDIUM) 
    {
        m_SimResult->SaveProgress(100.0,false);
        m_SimResult->WriteProgressLine(std::cout,sim_hop_limit,(m_ParamSet->m_PhiGradient != 0.0),false);
        std::cout << " --- SIMULATION END ---" << std::endl;
        std::cout << "(simulation runtime: ";
        GF::WriteDuration(std::cout, std::chrono::steady_clock::now() - runtime_start);
        std::cout << ")" << std::endl;
    }

    // Output of final results (in different variants depending on verbosity)
    if (m_VL == Verbosity::MINIMUM) 
    {
        std::cout << hop_counter - pre_hop_limit - eq_hop_limit << " hops, ";
        std::cout << m_SimResult->m_TotalTime << " s -> ";
        if (m_ParamSet->m_PhiGradient == 0.0)
        {
            std::cout << m_SimResult->m_DiffusionCoefficient << " cm2/s" << std::endl;
        }
        else
        {
            std::cout << m_SimResult->m_DriftConductivity << " S/m, " 
            << m_SimResult->m_DiffusionCoefficientParallel << " cm2/s, "
            << m_SimResult->m_DiffusionCoefficientTransverse << " cm2/s" << std::endl;
        }
    }
    if ((m_VL == Verbosity::MAXIMUM) && (eq_hop_limit != 0))
    {
        std::cout << "Final results:" << std::endl;
        std::cout << "  Simulated hops: " << m_SimResult->m_TotalHops 
            << " [incl. eq.: " << incl_eq_result->m_TotalHops << "]" << std::endl;
        std::cout << "  Oscillating hops: " 
            << m_SimResult->m_TotalHops - m_SimResult->m_NonOscHops << " (" 
            << 100.0 * (1.0 - m_SimResult->m_NonOscHopRatio) << " %" << " of hops) [incl. eq.: " 
            << incl_eq_result->m_TotalHops - incl_eq_result->m_NonOscHops << ", " 
            << 100.0 * (1.0 - incl_eq_result->m_NonOscHopRatio) << " %" << " of hops]" << std::endl;
        std::cout << "  Non-oscillating hops: " << m_SimResult->m_NonOscHops << " (" 
            << 100.0 * m_SimResult->m_NonOscHopRatio << " %" << " of hops) [incl. eq.: " 
            << incl_eq_result->m_NonOscHops << ", " 
            << 100.0 * incl_eq_result->m_NonOscHopRatio << " %" << " of hops]" << std::endl;
        std::cout << "  Simulated timespan: " << m_SimResult->m_TotalTime << " s [incl. eq.: " 
            << incl_eq_result->m_TotalTime << " s]" << std::endl;
        if (m_ParamSet->m_PhiGradient == 0.0)
        {
            std::cout << "  Self-diffusion coefficient D: " << m_SimResult->m_DiffusionCoefficient 
                << " cm2/s [incl. eq.: " << incl_eq_result->m_DiffusionCoefficient << " cm2/s]" << std::endl;
        }
        else
        {
            std::cout << "  Drift conductivity: " << m_SimResult->m_DriftConductivity 
                << " S/m [incl. eq.: " << incl_eq_result->m_DriftConductivity << " S/m]" << std::endl;
            std::cout << "  Drift mobility: " << m_SimResult->m_DriftMobility 
                << " cm2/Vs [incl. eq.: " << incl_eq_result->m_DriftMobility << " cm2/Vs]" << std::endl;
            std::cout << "  Self-diffusion coefficient D: " << m_SimResult->m_DiffusionCoefficient 
                << " cm2/s [incl. eq.: " << incl_eq_result->m_DiffusionCoefficient 
                << " cm2/s] (only valid if Dx approx. equal to Dyz)" << std::endl;
            std::cout << "  Self-diffusion coefficient Dx (parallel): " 
                << m_SimResult->m_DiffusionCoefficientParallel << " cm2/s [incl. eq.: " 
                << incl_eq_result->m_DiffusionCoefficientParallel << " cm2/s]" << std::endl;
            std::cout << "  Self-diffusion coefficient Dyz (transverse): " 
                << m_SimResult->m_DiffusionCoefficientTransverse << " cm2/s [incl. eq.: " 
                << incl_eq_result->m_DiffusionCoefficientTransverse << " cm2/s]" << std::endl;
            std::cout << "  Haven ratio: " << m_SimResult->m_HavenRatio 
                << " [incl. eq.: " << incl_eq_result->m_HavenRatio 
                << "] (only valid if Dx approx. equal to Dyz)" << std::endl;
            std::cout << "  Haven ratio (parallel = Dx / Dchemx): " << m_SimResult->m_HavenRatioParallel 
                << " [incl. eq.: " << incl_eq_result->m_HavenRatioParallel << "]" << std::endl;
            std::cout << "  Haven ratio (transverse = Dyz / Dchemx): " << m_SimResult->m_HavenRatioTransverse 
                << " [incl. eq.: " << incl_eq_result->m_HavenRatioTransverse 
                << "] (only valid for weak electric field)" << std::endl;
        }
        std::cout << "  Partial entropy: " << m_SimResult->m_PartialEntropy 
            << " eV/K (expected from DOS: " 
            << m_DOS->GetPartialEntropy(m_SimResult->m_EffChemPot,m_SimResult->m_EffTemp) << " eV/K)" << std::endl;
        if (m_ParamSet->m_TeffFit)
        {
            std::cout << "  Effective chemical potential: " << m_SimResult->m_EffChemPot 
            << " eV [after eq.: " << eq_result->m_EffChemPot << " eV]" << std::endl;
            std::cout << "  Effective temperature: " << m_SimResult->m_EffTemp 
            << " K [after eq.: " << eq_result->m_EffTemp << " K]" << std::endl;
        }
        std::cout << "  Effective charge carriers: " << m_SimResult->m_EffCarriers 
            << " [after eq.: " << eq_result->m_EffCarriers << "]" << std::endl;
        std::cout << "  Effective charge carrier density: " << m_SimResult->m_EffCarrierDensity 
            << " 1/cm3 (expected from DOS: " 
            << m_DOS->GetEffCarrierDensity(m_SimResult->m_EffChemPot,m_SimResult->m_EffTemp)
            << " 1/cm3) [after eq.: " << eq_result->m_EffCarrierDensity << " 1/cm3]" << std::endl;
        std::cout << "  Mean displacement of eff. charge carriers: (" 
            << m_SimResult->m_MeanDisp_x << ", " << m_SimResult->m_MeanDisp_y << ", " 
            << m_SimResult->m_MeanDisp_z << ") nm [incl. eq.: ("
            << incl_eq_result->m_MeanDisp_x << ", " << incl_eq_result->m_MeanDisp_y << ", " 
            << incl_eq_result->m_MeanDisp_z << ") nm]" << std::endl;
        std::cout << "  Average hops per eff. charge carrier: " 
            << static_cast<double>(m_SimResult->m_TotalHops)/m_SimResult->m_EffCarriers << " [incl. eq.: " 
            << static_cast<double>(incl_eq_result->m_TotalHops)/m_SimResult->m_EffCarriers << "]" << std::endl;
        std::cout << "  Average non-osc. hops per eff. charge carrier: " 
            << static_cast<double>(m_SimResult->m_NonOscHops)/m_SimResult->m_EffCarriers << " [incl. eq.: " 
            << static_cast<double>(incl_eq_result->m_NonOscHops)/m_SimResult->m_EffCarriers << "]" << std::endl;
        if ((m_SimResult->m_NonOscHops != 0) && (incl_eq_result->m_NonOscHops != 0))
        {
            std::cout << "  Percentage of non-osc. hops in +x direction: " << 100.0 * m_SimResult->m_InXDirNonOscRatio 
                << " % [incl. eq.: " << 100.0 * incl_eq_result->m_InXDirNonOscRatio << " %]" << std::endl;
        }
        if (m_ParamSet->m_PhiGradient != 0.0)
        {
            std::cout << "  Average electric field contribution to hop energy difference: " 
                << m_SimResult->m_MeanFieldContribution << " eV [incl. eq.: " 
                << incl_eq_result->m_MeanFieldContribution << " eV]" << std::endl;
        }
        std::cout << "-- all following values include the equilibration --" << std::endl;
        std::cout << "  Mobile electrons (> 0 non-osc. hops): " << m_SimResult->m_MobileElectrons << " ("
            << 100.0 * static_cast<double>(m_SimResult->m_MobileElectrons)/static_cast<double>(m_SimResult->m_ElectronCount) 
            << " %" << " of electrons) [after eq.: "<< eq_result->m_MobileElectrons << ", " 
            << 100.0 * static_cast<double>(eq_result->m_MobileElectrons)/static_cast<double>(m_SimResult->m_MobileElectrons) << " %" << " of final value]" << std::endl;
        const double mob_eff_ratio = static_cast<double>(m_SimResult->m_MobileElectrons) / m_SimResult->m_EffCarriers;
        std::cout << "    Ratio of mobile electrons to eff. charge carriers: " << mob_eff_ratio << std::endl; 
        std::cout << "    Mean displacement of mobile electrons: (" 
            << m_SimResult->m_MeanDisp_x / mob_eff_ratio << ", " << m_SimResult->m_MeanDisp_y / mob_eff_ratio << ", " 
            << m_SimResult->m_MeanDisp_z / mob_eff_ratio << ") nm [incl. eq.: ("
            << incl_eq_result->m_MeanDisp_x / mob_eff_ratio << ", " << incl_eq_result->m_MeanDisp_y / mob_eff_ratio << ", " 
            << incl_eq_result->m_MeanDisp_z / mob_eff_ratio << ") nm]" << std::endl;
        std::cout << "    Average hops per mobile electron: " 
            << static_cast<double>(m_SimResult->m_TotalHops)/static_cast<double>(m_SimResult->m_MobileElectrons) << " [incl. eq.: " 
            << static_cast<double>(incl_eq_result->m_TotalHops)/static_cast<double>(m_SimResult->m_MobileElectrons) << "]" << std::endl;
        std::cout << "    Average non-osc. hops per mobile electron: " 
            << static_cast<double>(m_SimResult->m_NonOscHops)/static_cast<double>(m_SimResult->m_MobileElectrons) << " [incl. eq.: " 
            << static_cast<double>(incl_eq_result->m_NonOscHops)/static_cast<double>(m_SimResult->m_MobileElectrons) << "]" << std::endl;
        std::cout << "  Stationary electrons (0 hops): " << m_SimResult->m_ZeroHopElectrons << " (" 
            << 100.0*static_cast<double>(m_SimResult->m_ZeroHopElectrons)/static_cast<double>(m_SimResult->m_ElectronCount) 
            << " %" << " of electrons) [after eq.: " << eq_result->m_ZeroHopElectrons << ", " 
            << 100.0*static_cast<double>(eq_result->m_ZeroHopElectrons)/static_cast<double>(m_SimResult->m_ZeroHopElectrons) << " %" << " of final value]" << std::endl;
        std::cout << "  Oscillating electrons (only osc. hops): " << m_SimResult->m_OscElectrons << " (" 
            << 100.0*static_cast<double>(m_SimResult->m_OscElectrons)/static_cast<double>(m_SimResult->m_ElectronCount) 
            << " %" << " of electrons) [after eq.: " << eq_result->m_OscElectrons << ", " 
            << 100.0*static_cast<double>(eq_result->m_OscElectrons)/static_cast<double>(m_SimResult->m_OscElectrons) << " %" << " of final value]" << std::endl;
    }
    if ((m_VL != Verbosity::MINIMUM) &&
        ((m_VL != Verbosity::MAXIMUM) || (eq_hop_limit == 0)))
    {
        std::cout << "Final results:" << std::endl;
        std::cout << "  Simulated hops: " << m_SimResult->m_TotalHops;
        std::cout << ", Oscillating hops: " 
            << m_SimResult->m_TotalHops - m_SimResult->m_NonOscHops << " (" 
            << 100.0 * (1.0 - m_SimResult->m_NonOscHopRatio) << " %" << " of hops)";
        std::cout << ", Non-oscillating hops: " << m_SimResult->m_NonOscHops << " (" 
            << 100.0 * m_SimResult->m_NonOscHopRatio << " %" << " of hops)" << std::endl;
        std::cout << "  Simulated timespan: " << m_SimResult->m_TotalTime << " s" << std::endl;
        if (m_ParamSet->m_PhiGradient == 0.0)
        {
            std::cout << "  Self-diffusion coefficient D: " << m_SimResult->m_DiffusionCoefficient << " cm2/s" << std::endl;
        }
        else
        {
            std::cout << "  Drift conductivity: " << m_SimResult->m_DriftConductivity << " S/m" << std::endl;
            std::cout << "  Drift mobility: " << m_SimResult->m_DriftMobility << " cm2/Vs" << std::endl;
            std::cout << "  Self-diffusion coefficient D: " << m_SimResult->m_DiffusionCoefficient 
                << " cm2/s (only valid if Dx approx. equal to Dyz)" << std::endl;
            std::cout << "  Self-diffusion coefficient Dx (parallel): " 
                << m_SimResult->m_DiffusionCoefficientParallel << " cm2/s" << std::endl;
            std::cout << "  Self-diffusion coefficient Dyz (transverse): " 
                << m_SimResult->m_DiffusionCoefficientTransverse << " cm2/s" << std::endl;
            std::cout << "  Haven ratio: " << m_SimResult->m_HavenRatio 
                << " (only valid if Dx approx. equal to Dyz)" << std::endl;
            std::cout << "  Haven ratio (parallel = Dx / Dchemx): " 
                << m_SimResult->m_HavenRatioParallel << std::endl;
            std::cout << "  Haven ratio (transverse = Dyz / Dchemx): " 
                << m_SimResult->m_HavenRatioParallel << " (only valid for weak electric field)" << std::endl;
        }
        std::cout << "  Partial entropy: " << m_SimResult->m_PartialEntropy 
            << " eV/K (expected from DOS: " 
            << m_DOS->GetPartialEntropy(m_SimResult->m_EffChemPot,m_SimResult->m_EffTemp) << " eV/K)" << std::endl;
        if (m_ParamSet->m_TeffFit)
        {
            std::cout << "  Effective chemical potential: " << m_SimResult->m_EffChemPot << " eV" << std::endl;
            std::cout << "  Effective temperature: " << m_SimResult->m_EffTemp << " K" << std::endl;
        }
        std::cout << "  Effective charge carriers: " << m_SimResult->m_EffCarriers << std::endl;
        std::cout << "  Effective charge carrier density: " << m_SimResult->m_EffCarrierDensity 
            << " 1/cm3 (expected from DOS: " 
            << m_DOS->GetEffCarrierDensity(m_SimResult->m_EffChemPot,m_SimResult->m_EffTemp)
            << " 1/cm3)" << std::endl;
        std::cout << "  Mean displacement of eff. charge carriers: (" 
            << m_SimResult->m_MeanDisp_x << ", " << m_SimResult->m_MeanDisp_y << ", " 
            << m_SimResult->m_MeanDisp_z << ") nm" << std::endl;
        std::cout << "  Average hops per eff. charge carrier: " 
            << static_cast<double>(m_SimResult->m_TotalHops)/m_SimResult->m_EffCarriers << std::endl;
        std::cout << "  Average non-osc. hops per eff. charge carrier: " 
            << static_cast<double>(m_SimResult->m_NonOscHops)/m_SimResult->m_EffCarriers << std::endl;
        if (m_SimResult->m_NonOscHops != 0)
        {
            std::cout << "  Percentage of non-osc. hops in +x direction: " 
                << 100.0 * m_SimResult->m_InXDirNonOscRatio << " %" << std::endl;
        }
        if (m_ParamSet->m_PhiGradient != 0.0)
        {
            std::cout << "  Average electric field contribution to hop energy difference: " 
                << m_SimResult->m_MeanFieldContribution << " eV" << std::endl;
        }
        if (eq_hop_limit != 0)
            std::cout << "-- all following values include the equilibration --" << std::endl;
        std::cout << "  Mobile electrons (> 0 non-osc. hops): " << m_SimResult->m_MobileElectrons << " ("
            << 100.0 * static_cast<double>(m_SimResult->m_MobileElectrons)/static_cast<double>(m_SimResult->m_ElectronCount) << " %" << " of electrons)" << std::endl;
        const double mob_eff_ratio = static_cast<double>(m_SimResult->m_MobileElectrons) / m_SimResult->m_EffCarriers;
        std::cout << "    Ratio of mobile electrons to eff. charge carriers: " << mob_eff_ratio << std::endl; 
        std::cout << "    Mean displacement of mobile electrons: (" 
            << m_SimResult->m_MeanDisp_x / mob_eff_ratio << ", " << m_SimResult->m_MeanDisp_y / mob_eff_ratio << ", " 
            << m_SimResult->m_MeanDisp_z / mob_eff_ratio << ") nm" << std::endl;
        std::cout << "    Average hops per mobile electron: " 
            << static_cast<double>(m_SimResult->m_TotalHops)/static_cast<double>(m_SimResult->m_MobileElectrons) << std::endl;
        std::cout << "    Average non-osc. hops per mobile electron: " 
            << static_cast<double>(m_SimResult->m_NonOscHops)/static_cast<double>(m_SimResult->m_MobileElectrons) << std::endl;
        std::cout << "  Stationary electrons (0 hops): " << m_SimResult->m_ZeroHopElectrons << " (" 
            << 100.0*static_cast<double>(m_SimResult->m_ZeroHopElectrons)/static_cast<double>(m_SimResult->m_ElectronCount) << " %" << " of electrons)" << std::endl;
        std::cout << "  Oscillating electrons (only osc. hops): " << m_SimResult->m_OscElectrons << " (" 
            << 100.0*static_cast<double>(m_SimResult->m_OscElectrons)/static_cast<double>(m_SimResult->m_ElectronCount) << " %" << " of electrons)" << std::endl;
    }
}

// Execute KMC hops until hop counter equals the limit
void MC::TEngine::KMCLoop(std::uint64_t& hop_counter, const std::uint64_t hop_limit) noexcept
{
    // No checks because only internally called by RunSimulation

    // This function is not expected to throw exceptions (for example there should be no major memory 
    // allocations because already enough memory was reserved for m_RandomTimes-vectors of 
    // electrons), but the main reason that it is specified as noexcept is to exclude it from 
    // surrounding exception handling for improved performance
    // (calls std::terminate immidiately if an exception occurs in this function)

    // Grab values and references (t-prefix) to avoid repeated *this dereferencing
    double& t_total_time {m_TotalTime};
    const std::uint32_t t_state_count {m_ParamSet->m_StateCount};
    const std::vector<TLocalState>& t_structure {m_Structure};
    std::vector<TElectron>& t_electrons {m_Electrons};
    std::mt19937_64& t_rndgen {m_RndGen};
    std::uniform_real_distribution<double> probability (0.0, 1.0);
    
    // Inner KMC-Loop = performance-critical
    while (hop_counter < hop_limit)
    {
        // Electron with minimal time
        TElectron& min_e {*std::min_element(t_electrons.begin(),t_electrons.end(),
            [] (const TElectron& e_left, const TElectron& e_right) -> bool
            {
                return e_left.m_MinTime < e_right.m_MinTime;
            })};

        // Subtract time step for all electrons and advance total time
        // (also iterates over min_e and invalidates its time values which is no problem)
        {
            const double time_step {min_e.m_MinTime};
            for (TElectron& it_e : t_electrons)
            {
                it_e.m_AccTimeDecay += time_step;
                it_e.m_MinTime -= time_step;
            }
            t_total_time += time_step;
        }

        // Update electrons around the start position (which becomes empty)
        // (also iterates over the minimal path which is no problem because end state is still empty)
        const TLocalState& old_state {t_structure[min_e.m_CurrentStateID]};
        for (const TPath& it_p : old_state.m_Paths)
        {
            const TLocalState& state {t_structure[it_p.m_StateID]};
            if (state.m_ElectronID < t_state_count)
            {
                TElectron& env_e {t_electrons[state.m_ElectronID]};

                // Subtract accumulated time decay from path times before new path is added
                const double& decay {env_e.m_AccTimeDecay};
                for (double& it_t : env_e.m_RandomTimes)
                {
                    it_t -= decay;
                }
                env_e.m_AccTimeDecay = 0.0;

                // Calculate time for the newly available path
                double& reverse_time {env_e.m_RandomTimes[it_p.m_ReverseID]};
                reverse_time = -log(probability(t_rndgen)) * state.m_Paths[it_p.m_ReverseID].m_Time;

                // Check if its the new minimal time (cannot be the old minimum)
                if (reverse_time < env_e.m_MinTime)
                {
                    env_e.m_MinTime = reverse_time;
                    env_e.m_MinIndex = it_p.m_ReverseID;
                }
            }
        }

        // Get reference to chosen path and new state
        const TPath& min_path {old_state.m_Paths[min_e.m_MinIndex]};
        const TLocalState& new_state {t_structure[min_path.m_StateID]};

        // Detect if this is an oscillatory hop (back to previous position)
        if (min_e.m_LastStateID == min_path.m_StateID)
        {
            min_e.m_OscHopCount += 2U;
            ++(min_path.m_OscHopCount);
            ++(new_state.m_Paths[min_path.m_ReverseID].m_OscHopCount);
        }
        else
        {
            min_e.m_LastStateID = min_e.m_CurrentStateID;
        }

        // Transfer electron to new state (and set old state empty)
        min_e.m_CurrentStateID = min_path.m_StateID;
        ++(min_path.m_HopCount);
        new_state.m_ElectronID = old_state.m_ElectronID;
        old_state.m_ElectronID = t_state_count;

        // Add timespan that old state was occupied
        old_state.m_OccTime += t_total_time - old_state.m_LastHopTime;
        old_state.m_LastHopTime = t_total_time;
        new_state.m_LastHopTime = t_total_time;

        // Detect if it is the first hop of this electron
        if (min_e.m_HopCount == 0) min_e.m_FirstHopTime = t_total_time;

        // Calculate displacement and hopping counters
        {
            const double& old_x {old_state.m_Pos_x};
            double new_x {new_state.m_Pos_x};
            if (old_x - new_x > 0.5) new_x += 1.0;
            if (new_x - old_x > 0.5) new_x -= 1.0;
            const double& old_y {old_state.m_Pos_y};
            double new_y {new_state.m_Pos_y};
            if (old_y - new_y > 0.5) new_y += 1.0;
            if (new_y - old_y > 0.5) new_y -= 1.0;
            const double& old_z {old_state.m_Pos_z};
            double new_z {new_state.m_Pos_z};
            if (old_z - new_z > 0.5) new_z += 1.0;
            if (new_z - old_z > 0.5) new_z -= 1.0;

            min_e.m_Disp_x += new_x - old_x;
            min_e.m_Disp_y += new_y - old_y;
            min_e.m_Disp_z += new_z - old_z;
            ++(min_e.m_HopCount);
            ++hop_counter;
        }

        // Generate new randomized transition times and find minimum
        {
            double min_time {Constant::occupied_path};
            std::uint32_t min_index {0U};
            const std::size_t path_count {new_state.m_Paths.size()};
            min_e.m_RandomTimes.resize(path_count,0.0);
            for (std::uint32_t i = 0U; i < path_count; ++i)
            {
                if (t_structure[new_state.m_Paths[i].m_StateID].m_ElectronID < t_state_count)
                {
                    min_e.m_RandomTimes[i] = Constant::occupied_path;
                }
                else
                {
                    min_e.m_RandomTimes[i] = -log(probability(t_rndgen))*new_state.m_Paths[i].m_Time;
                    if (min_e.m_RandomTimes[i] < min_time)
                    {
                        min_time = min_e.m_RandomTimes[i];
                        min_index = i;
                    }
                }
            }
            min_e.m_MinTime = min_time;
            min_e.m_MinIndex = min_index;
            min_e.m_AccTimeDecay = 0.0;
        }

        // Update electrons around the end position (which becomes occupied)
        for (const TPath& it_p : new_state.m_Paths)
        {
            const TLocalState& state {t_structure[it_p.m_StateID]};
            if (state.m_ElectronID < t_state_count)
            {
                TElectron& env_e {t_electrons[state.m_ElectronID]};
                env_e.m_RandomTimes[it_p.m_ReverseID] = Constant::occupied_path;

                // Check if it was the minimal time (cannot remain minimum = find new minimum)
                if (env_e.m_MinIndex == it_p.m_ReverseID)
                {
                    std::vector<double>::const_iterator min_it 
                        {std::min_element(env_e.m_RandomTimes.cbegin(),env_e.m_RandomTimes.cend())};
                    env_e.m_MinIndex = static_cast<std::uint32_t>(std::distance(env_e.m_RandomTimes.cbegin(),min_it));
	                env_e.m_MinTime = *min_it - env_e.m_AccTimeDecay;
                }
            }
        }
    } 
    // end of KMC loop
}

// Generate statistics at end of equilibration (incl. optional auto-adjust of cut-offs)
void MC::TEngine::CalculatePostEquilibrationStatistics()
{
    // No checks because only internally called by RunSimulation

    if ((m_VL < Verbosity::MEDIUM) && (!m_ParamSet->m_CutoffAutoAdjust)) return;

    const double spfac = m_DOS->GetSpatialFactor();
    const double enmax = m_ParamSet->m_MaxStateEnergy;
    const double enfac = m_DOS->GetEnergyFactor();

    // Analysis of used paths
    std::uint64_t total_paths = 0;
    std::uint64_t used_paths = 0;
    GF::TMinMax<double> dist_per_used_path;
    GF::TMinMax<double> ediff_per_used_path;
    for (const auto& state : m_Structure)
    {
        for (const auto& path : state.m_Paths)
        {
            ++total_paths;
            if (path.m_HopCount == 0) continue;
            ++used_paths;

            const double& old_x {state.m_Pos_x};
            double new_x {m_Structure[path.m_StateID].m_Pos_x};
            if (old_x - new_x > 0.5) new_x += 1.0;
            if (new_x - old_x > 0.5) new_x -= 1.0;
            const double& old_y {state.m_Pos_y};
            double new_y {m_Structure[path.m_StateID].m_Pos_y};
            if (old_y - new_y > 0.5) new_y += 1.0;
            if (new_y - old_y > 0.5) new_y -= 1.0;
            const double& old_z {state.m_Pos_z};
            double new_z {m_Structure[path.m_StateID].m_Pos_z};
            if (old_z - new_z > 0.5) new_z += 1.0;
            if (new_z - old_z > 0.5) new_z -= 1.0;
            double t_path_dist = sqrt((new_x - old_x)*(new_x - old_x) + (new_y - old_y)*
                (new_y - old_y) + (new_z - old_z)*(new_z - old_z)) * spfac;

            double t_path_ediff = (m_Structure[path.m_StateID].m_Energy - state.m_Energy) * enfac;

            dist_per_used_path.check(t_path_dist);
            ediff_per_used_path.check(t_path_ediff);
        }
    }
    
    if (m_VL >= Verbosity::MEDIUM)
    {
        std::cout << "Equilibration statistics:" << std::endl;
        std::cout << "  Used paths during eq.: " << used_paths << " (" 
            << 100.0 * static_cast<double>(used_paths) / static_cast<double>(total_paths) << " %" << " of paths)" << std::endl;
        std::cout << "  Path distance range used during eq.: " << dist_per_used_path << " nm" << std::endl;
        std::cout << "  Path energy difference range used during eq.: " << ediff_per_used_path 
            << " eV (without electric field contribution)" << std::endl;
        std::cout << "  State energy range used during eq.: [" << m_SimResult->m_MinUsedStateEnergy << ", "
            << m_SimResult->m_MaxUsedStateEnergy << "] eV" << std::endl;
    }

    // Adjust cut-offs and delete exceeding paths
    if (m_ParamSet->m_CutoffAutoAdjust)
    {
        const double dist_perc = (m_ParamSet->m_DistCutoffAdjustPercentage != 0.0) ? 
            m_ParamSet->m_DistCutoffAdjustPercentage : Constant::autocutoff_inc;
        const double ediff_perc = (m_ParamSet->m_EdiffCutoffAdjustPercentage != 0.0) ? 
            m_ParamSet->m_EdiffCutoffAdjustPercentage : Constant::autocutoff_inc;

        if (m_VL >= Verbosity::MEDIUM)
        {
            std::cout << "Auto-adjustment of cut-offs based on used range (distance: + " << dist_perc 
                << " %, energy diff.: + " << ediff_perc << "%): " << std::flush;
        }
        
        // Calculate new cut-offs
        const double new_max_dist = dist_per_used_path.max 
            * (1.0 + dist_perc / 100.0);
        const double new_max_ediff = std::max(fabs(ediff_per_used_path.min),fabs(ediff_per_used_path.max)) 
            * (1.0 + ediff_perc / 100.0);

        // Check if cut-offs remain unchanged
        if ((new_max_dist >= m_SimResult->m_MaxPathDist) && (new_max_ediff >= m_SimResult->m_MaxPathEdiff))
        {
            if (m_VL >= Verbosity::MEDIUM)
            {
                std::cout << "done" << std::endl;
                std::cout << "  No adjustment: new cut-offs would be equal or greater than current cut-offs." << std::endl;
            }
        }
        else
        {
            // Calculate cut-offs in relative units
            const double new_dist_cutoff = new_max_dist / spfac;
            const double new_ediff_cutoff = new_max_ediff / enfac;

            // Check if at least two paths per state remain
            bool too_few_remaining_paths = false;
            for (const auto& state : m_Structure)
            {
                std::uint32_t new_path_count = 0;
                for (const auto& path : state.m_Paths)
                {
                    const double& old_x {state.m_Pos_x};
                    double new_x {m_Structure[path.m_StateID].m_Pos_x};
                    if (old_x - new_x > 0.5) new_x += 1.0;
                    if (new_x - old_x > 0.5) new_x -= 1.0;
                    const double& old_y {state.m_Pos_y};
                    double new_y {m_Structure[path.m_StateID].m_Pos_y};
                    if (old_y - new_y > 0.5) new_y += 1.0;
                    if (new_y - old_y > 0.5) new_y -= 1.0;
                    const double& old_z {state.m_Pos_z};
                    double new_z {m_Structure[path.m_StateID].m_Pos_z};
                    if (old_z - new_z > 0.5) new_z += 1.0;
                    if (new_z - old_z > 0.5) new_z -= 1.0;
                    double t_path_dist = sqrt((new_x - old_x)*(new_x - old_x) + (new_y - old_y)*
                        (new_y - old_y) + (new_z - old_z)*(new_z - old_z));

                    double t_path_ediff = fabs(m_Structure[path.m_StateID].m_Energy - state.m_Energy);

                    if ((t_path_dist <= new_dist_cutoff) && (t_path_ediff <= new_ediff_cutoff))
                        ++new_path_count;
                }
                if (new_path_count < 2) 
                {
                    too_few_remaining_paths = true;
                    break;
                }
            }

            if (too_few_remaining_paths)
            {
                if (m_VL >= Verbosity::MEDIUM)
                {
                    std::cout << "done" << std::endl;
                    std::cout << "  No adjustment: new cut-offs would lead to states with less than two paths." << std::endl;
                }
            }
            else
            {
                // Delete path pairs that exceed the new cut-offs
                // Attention: update all inter-linking IDs and ordered lists
                std::uint64_t new_total_paths = 0;
                std::uint64_t deleted_paths = 0;
                GF::TMinMaxMean<std::uint32_t> paths_per_state;
                for (TLocalState& state : m_Structure)
                {
                    const double& old_x {state.m_Pos_x};
                    const double& old_y {state.m_Pos_y};
                    const double& old_z {state.m_Pos_z};

                    std::uint32_t k = 0;
                    while (k < state.m_Paths.size())
                    {
                        TLocalState& new_state = m_Structure[state.m_Paths[k].m_StateID];
                        const std::uint32_t rev_k = state.m_Paths[k].m_ReverseID;

                        double new_x {new_state.m_Pos_x};
                        if (old_x - new_x > 0.5) new_x += 1.0;
                        if (new_x - old_x > 0.5) new_x -= 1.0;
                        double new_y {new_state.m_Pos_y};
                        if (old_y - new_y > 0.5) new_y += 1.0;
                        if (new_y - old_y > 0.5) new_y -= 1.0;
                        double new_z {new_state.m_Pos_z};
                        if (old_z - new_z > 0.5) new_z += 1.0;
                        if (new_z - old_z > 0.5) new_z -= 1.0;
                        double t_path_dist = sqrt((new_x - old_x)*(new_x - old_x) + (new_y - old_y)*
                            (new_y - old_y) + (new_z - old_z)*(new_z - old_z));

                        double t_path_ediff = fabs(new_state.m_Energy - state.m_Energy);

                        if ((t_path_dist <= new_dist_cutoff) && (t_path_ediff <= new_ediff_cutoff))
                        {
                            ++k;
                            continue;
                        }
                            
                        // Delete state.m_Paths[k] and new_state.m_Paths[rev_k]
                        state.m_Paths.erase(state.m_Paths.begin() + k);
                        new_state.m_Paths.erase(new_state.m_Paths.begin() + rev_k);
                        deleted_paths += 2;

                        // Update electrons that occupy these states
                        if (state.m_ElectronID < m_ParamSet->m_StateCount)
                        {
                            TElectron& electron = m_Electrons[state.m_ElectronID];

                            // Delete electron.m_RandomTimes[k] (list in same order as m_Paths)
                            electron.m_RandomTimes.erase(electron.m_RandomTimes.begin() + k);

                            // Find new minimum time if deleted path was the minimum
                            if (electron.m_MinIndex == k)
                            {
                                std::vector<double>::const_iterator min_it 
                                    {std::min_element(electron.m_RandomTimes.cbegin(),electron.m_RandomTimes.cend())};
                                electron.m_MinIndex = static_cast<std::uint32_t>(std::distance(electron.m_RandomTimes.cbegin(),min_it));
                                electron.m_MinTime = *min_it - electron.m_AccTimeDecay;
                            }
                            // Reduce minimum index if deleted path was below in the list
                            else if (electron.m_MinIndex > k)
                            {
                                --(electron.m_MinIndex);
                            } 
                        }
                        if (new_state.m_ElectronID < m_ParamSet->m_StateCount)
                        {
                            TElectron& electron = m_Electrons[new_state.m_ElectronID];

                            // Delete electron.m_RandomTimes[rev_k] (list in same order as m_Paths)
                            electron.m_RandomTimes.erase(electron.m_RandomTimes.begin() + rev_k);

                            // Find new minimum time if deleted path was the minimum
                            if (electron.m_MinIndex == rev_k)
                            {
                                std::vector<double>::const_iterator min_it 
                                    {std::min_element(electron.m_RandomTimes.cbegin(),electron.m_RandomTimes.cend())};
                                electron.m_MinIndex = static_cast<std::uint32_t>(std::distance(electron.m_RandomTimes.cbegin(),min_it));
                                electron.m_MinTime = *min_it - electron.m_AccTimeDecay;
                            }
                            // Reduce minimum index if deleted path was below in the list
                            else if (electron.m_MinIndex > rev_k)
                            {
                                --(electron.m_MinIndex);
                            } 
                        }
                        
                        // Update the m_ReverseIDs of paths that are reverse to the paths that shifted in the list (>= k, >= rev_k)
                        // (reduces their m_ReverseIDs by one)
                        for (std::uint32_t corr_k = k; corr_k < state.m_Paths.size(); ++corr_k)
                        {
                            m_Structure[state.m_Paths[corr_k].m_StateID].m_Paths[state.m_Paths[corr_k].m_ReverseID].m_ReverseID = corr_k;
                        }
                        for (std::uint32_t corr_k = rev_k; corr_k < new_state.m_Paths.size(); ++corr_k)
                        {
                            m_Structure[new_state.m_Paths[corr_k].m_StateID].m_Paths[new_state.m_Paths[corr_k].m_ReverseID].m_ReverseID = corr_k;
                        }
                    }
                    paths_per_state.check(static_cast<std::uint32_t>(state.m_Paths.size()));
                    new_total_paths += state.m_Paths.size();
                }

                if (paths_per_state.min < 2) 
                    throw EX::TInvalidStatus("States with less than two paths detected.",__func__);

                m_SimResult->m_MaxPathCount = paths_per_state.max;
                m_SimResult->m_MeanPathCount = paths_per_state.mean;

                // Free memory of deleted paths
                for (TLocalState& state : m_Structure)
                {
                    state.m_Paths.shrink_to_fit();
                }

                // Adjust capacity of electrons time list to avoid allocations during simulation
                for (TElectron& electron : m_Electrons)
                {
                    electron.m_RandomTimes.shrink_to_fit();
                    try
                    {
                        electron.m_RandomTimes.reserve(m_SimResult->m_MaxPathCount);
                    }
                    catch(const std::bad_alloc& e)
                    {
                        throw EX::TOutOfMemory("Cannot reserve space for electron paths.",__func__,e.what());
                    }
                }

                if (m_VL >= Verbosity::MEDIUM)
                {
                    std::cout << "done" << std::endl;
                }
                if (new_max_dist < m_SimResult->m_MaxPathDist)
                {
                    m_SimResult->m_MaxPathDist = new_max_dist;
                    if (m_VL >= Verbosity::MEDIUM)
                        std::cout << "  New cut-off for distance: " << m_SimResult->m_MaxPathDist << " nm" << std::endl;   
                }
                else
                {
                    if (m_VL >= Verbosity::MEDIUM)
                        std::cout << "  No adjustment of distance: would be equal or greater than current cut-off." << std::endl;
                }
                if (new_max_ediff < m_SimResult->m_MaxPathEdiff)
                {
                    m_SimResult->m_MaxPathEdiff = new_max_ediff;
                    if (m_VL >= Verbosity::MEDIUM)
                        std::cout << "  New cut-off for state energy difference: " << m_SimResult->m_MaxPathEdiff << " eV" << std::endl;
                }
                else
                {
                    if (m_VL >= Verbosity::MEDIUM)
                        std::cout << "  No adjustment of energy difference: would be equal or greater than current cut-off." << std::endl;
                }
                if (m_VL >= Verbosity::MEDIUM)
                {
                    std::cout << "  Deleted paths: " << deleted_paths << " (" 
                        << 100.0*static_cast<double>(deleted_paths)/static_cast<double>(total_paths) << " %)" << std::endl;
                    std::cout << "  New total number of paths: " << new_total_paths << " (" 
                        << 100.0*static_cast<double>(new_total_paths)/static_cast<double>(total_paths) << " %)" << std::endl;
                    std::cout << "  New average number of paths per state: " << paths_per_state << std::endl;
                }
            }
        }
    }
}

// Fit state occupation with Fermi-Dirac distribution (optional) and calculate effective charge carrier density
void MC::TEngine::CalculateEffectiveCarrierDensity(bool fit_eff)
{
    // No checks because only internally called by InitializeSimulation and RunSimulation

    const double spfac = m_DOS->GetSpatialFactor();
    const double enmax = m_ParamSet->m_MaxStateEnergy;
    const double enfac = m_DOS->GetEnergyFactor();

    // Lambda function: calculate squared deviation of state occupations (input: chempot in eV, T in K)
    auto calc_deviation = [&] (double i_chempot, double i_T) -> double
    {
        const double t_chempot = (i_chempot - enmax)/enfac;
        const double t_kBT = Constant::kboltz*i_T/enfac;
        double sum = 0.0;
        for (const auto& state: m_Structure)
        {
            if (state.m_ElectronID < m_ParamSet->m_StateCount)
            {
                sum += (1.0 - 1.0/(exp((state.m_Energy - t_chempot)/t_kBT) + 1.0)) *
                    (1.0 - 1.0/(exp((state.m_Energy - t_chempot)/t_kBT) + 1.0));
            }
            else
            {
                sum += 1.0/((exp((state.m_Energy - t_chempot)/t_kBT) + 1.0) *
                    (exp((state.m_Energy - t_chempot)/t_kBT) + 1.0));
            }
        }
        return sum;
    };

    double effT = m_ParamSet->m_Temperature;
    double effChemPot = m_ParamSet->m_ChemPot;
    if (fit_eff)
    {
        // Find appropriate boundaries for golden-section search
        const double refDev = calc_deviation(effChemPot,effT);
        double lowT = effT;
        double lowDev = 0.0;
        double highT = effT;
        double highDev = 0.0;
        do
        {
            lowT *= 0.9;
            lowDev = calc_deviation(m_DOS->GetAdjustedFermiLevel(lowT),lowT);
        } while (lowDev <= refDev);
        do
        {
            highT *= 1.5;
            highDev = calc_deviation(m_DOS->GetAdjustedFermiLevel(highT),highT);
        } while (highDev <= refDev);

        // Golden-section search
        const double lowFactor = (3.0-sqrt(5.0))/2.0;   // approx: 0.382
        const double highFactor = (sqrt(5.0)-1.0)/2.0;  // approx: 0.618
        double lcT = lowT + lowFactor * (highT - lowT);
        double lcDev = calc_deviation(m_DOS->GetAdjustedFermiLevel(lcT),lcT);
        double hcT = lowT + highFactor * (highT - lowT);
        double hcDev = calc_deviation(m_DOS->GetAdjustedFermiLevel(hcT),hcT);
        while (highT - lowT > 0.1)
        {
            if (lcDev < hcDev)
            {
                highT = hcT;
                hcT = lcT;
                hcDev = lcDev;
                lcT = lowT + highFactor * (hcT - lowT);
                lcDev = calc_deviation(m_DOS->GetAdjustedFermiLevel(lcT),lcT);
            }
            else
            {
                lowT = lcT;
                lcT = hcT;
                lcDev = hcDev;
                hcT = lcT + lowFactor * (highT - lcT);
                hcDev = calc_deviation(m_DOS->GetAdjustedFermiLevel(hcT),hcT);
            }
        }
        if (lcDev < hcDev)
            effT = 0.5*(lowT + hcT);
        else
            effT = 0.5*(lcT + highT);
        effChemPot = m_DOS->GetAdjustedFermiLevel(effT);
    }

    // Calculate carrier density (and partial entropy) from states
    const double relChemPot = (effChemPot - enmax)/enfac;
    const double relkBT = Constant::kboltz*effT/enfac;
    double effCount = 0.0;
    double thermCount = 0.0;
    for (const auto& state: m_Structure)
    {
        effCount += (1.0/(exp((state.m_Energy - relChemPot)/relkBT) + 1.0)) *
            (1.0 - 1.0/(exp((state.m_Energy - relChemPot)/relkBT) + 1.0));

        thermCount += (1.0/(exp((state.m_Energy - relChemPot)/relkBT) + 1.0)) *
            (1.0 - 1.0/(exp((state.m_Energy - relChemPot)/relkBT) + 1.0)) *
            (state.m_Energy - relChemPot)/relkBT;
    }
    if (effCount <= 0.0)
        throw EX::TInvalidStatus("Zero or negative effective carrier density.",__func__);
    
    m_SimResult->m_EffChemPot = effChemPot;
    m_SimResult->m_EffTemp = effT;
    m_SimResult->m_EffCarriers = effCount;
    m_SimResult->m_EffCarrierDensity = effCount/pow(spfac,3.0)*1.0E21;
    m_SimResult->m_PartialEntropy = Constant::kboltz * thermCount / effCount;
}

// Extract result values from simulation data (optional: relative to reference result)
void MC::TEngine::CalculateResultValues(const MC::TResult* const ref_result)
{
    // No checks because only internally called by RunSimulation

    const double spfac = m_DOS->GetSpatialFactor();
    const double enmax = m_ParamSet->m_MaxStateEnergy;
    const double enfac = m_DOS->GetEnergyFactor();

    // Evaluation of electron statistics 
    // sum_energy = sum of state energies of occupied states (in eV)
    // nzero, nosc = number of electrons with zero hops and with only osc hops
    // sum_r, sum_x, sum_y, sum_z = sum of r,x,y,z-displacements (in relative units)
    // sum_sqx, sum_sqy, sum_sqz = sum of squared x,y,z-displacements (in relative units)
    double sum_energy = 0.0;
    std::uint32_t nzero = 0, nosc = 0;
    double sum_r = 0.0, sum_x = 0.0, sum_y = 0.0, sum_z = 0.0;
    double sum_sqx = 0.0, sum_sqy = 0.0, sum_sqz = 0.0;
    for (const auto &electron : m_Electrons)
    {
        sum_energy += m_Structure[electron.m_CurrentStateID].m_Energy * enfac + enmax;

        if (electron.m_HopCount == 0)
        {
            nzero++;
            continue;
        }
        if (electron.m_HopCount == electron.m_OscHopCount) nosc++;

        sum_r += sqrt(electron.m_Disp_x * electron.m_Disp_x 
            + electron.m_Disp_y * electron.m_Disp_y
            + electron.m_Disp_z * electron.m_Disp_z);
        sum_x += electron.m_Disp_x;
        sum_y += electron.m_Disp_y;
        sum_z += electron.m_Disp_z;
        sum_sqx += electron.m_Disp_x * electron.m_Disp_x;
        sum_sqy += electron.m_Disp_y * electron.m_Disp_y;
        sum_sqz += electron.m_Disp_z * electron.m_Disp_z;
    }
    
    // Evaluation of path statistics
    // sum_hops = sum of hop counts
    // sum_nonosc_hops = sum of non-oscillatory hops
    // sum_in_xdir_nonosc_hops = sum of non-oscillatory hops in +x direction (dx > 0)
    // avg_field_contrib = average field energy contribution to hop energy (in relative units)
    // max_dist = maximum distance of used paths (in relative units)
    // max_ediff = maximum absolute state energy difference of used paths (in relative units)
    // used_energy_range = min/max of used state energies (in relative units)
    std::uint64_t sum_hops = 0, sum_nonosc_hops = 0, sum_in_xdir_nonosc_hops = 0;
    double avg_field_contrib = 0.0, max_dist = 0.0, max_ediff = 0.0;
    GF::TMinMax<double> used_energy_range;
    for (const auto& state : m_Structure)
    {
        for (const auto& path : state.m_Paths)
        {
            if (path.m_HopCount == 0) continue;

            const double& old_x {state.m_Pos_x};
            double new_x {m_Structure[path.m_StateID].m_Pos_x};
            if (old_x - new_x > 0.5) new_x += 1.0;
            if (new_x - old_x > 0.5) new_x -= 1.0;
            const double& old_y {state.m_Pos_y};
            double new_y {m_Structure[path.m_StateID].m_Pos_y};
            if (old_y - new_y > 0.5) new_y += 1.0;
            if (new_y - old_y > 0.5) new_y -= 1.0;
            const double& old_z {state.m_Pos_z};
            double new_z {m_Structure[path.m_StateID].m_Pos_z};
            if (old_z - new_z > 0.5) new_z += 1.0;
            if (new_z - old_z > 0.5) new_z -= 1.0;
            const double t_path_dist = sqrt((new_x - old_x)*(new_x - old_x) + (new_y - old_y)*
                (new_y - old_y) + (new_z - old_z)*(new_z - old_z));

            const double t_path_ediff = fabs(m_Structure[path.m_StateID].m_Energy - state.m_Energy);

            sum_hops += path.m_HopCount;
            if (m_ParamSet->m_PhiGradient != 0.0)
            {
                avg_field_contrib += (- (new_x - old_x) * m_GradPhiEnergy - avg_field_contrib)
                    * (static_cast<double>(path.m_HopCount) / static_cast<double>(sum_hops));
            }

            sum_nonosc_hops += path.m_HopCount - path.m_OscHopCount;
            if (new_x > old_x)
            {
                sum_in_xdir_nonosc_hops += path.m_HopCount - path.m_OscHopCount;
            }

            if (t_path_dist > max_dist) max_dist = t_path_dist;
            if (t_path_ediff > max_ediff) max_ediff = t_path_ediff;

            used_energy_range.check(state.m_Energy);
            used_energy_range.check(m_Structure[path.m_StateID].m_Energy);
        }
    }

    // Electron-related properties
    m_SimResult->m_MobileElectrons = m_SimResult->m_ElectronCount - nzero - nosc;
    if (m_SimResult->m_MobileElectrons == 0)
        throw EX::TInvalidInput("Parameters lead to only oscillating electrons.");
    m_SimResult->m_ZeroHopElectrons = nzero;
    m_SimResult->m_OscElectrons = nosc;
	m_SimResult->m_MeanDisp = sum_r * spfac / m_SimResult->m_EffCarriers;
	m_SimResult->m_MeanDisp_x = sum_x * spfac / m_SimResult->m_EffCarriers;
	m_SimResult->m_MeanDisp_y = sum_y * spfac / m_SimResult->m_EffCarriers;
	m_SimResult->m_MeanDisp_z = sum_z * spfac / m_SimResult->m_EffCarriers;
	m_SimResult->m_MeanSquaredDisp_x = sum_sqx * spfac * spfac / m_SimResult->m_EffCarriers;
    m_SimResult->m_MeanSquaredDisp_y = sum_sqy * spfac * spfac / m_SimResult->m_EffCarriers;
	m_SimResult->m_MeanSquaredDisp_z = sum_sqz * spfac * spfac / m_SimResult->m_EffCarriers;	    
    m_SimResult->m_MeanDispVariance_x = m_SimResult->m_MeanSquaredDisp_x 
        - m_SimResult->m_EffCarriers / static_cast<double>(m_SimResult->m_MobileElectrons)
        * m_SimResult->m_MeanDisp_x * m_SimResult->m_MeanDisp_x;
    m_SimResult->m_MeanDispVariance_y = m_SimResult->m_MeanSquaredDisp_y 
        - m_SimResult->m_EffCarriers / static_cast<double>(m_SimResult->m_MobileElectrons)
        * m_SimResult->m_MeanDisp_y * m_SimResult->m_MeanDisp_y;
    m_SimResult->m_MeanDispVariance_z = m_SimResult->m_MeanSquaredDisp_z 
        - m_SimResult->m_EffCarriers / static_cast<double>(m_SimResult->m_MobileElectrons)
        * m_SimResult->m_MeanDisp_z * m_SimResult->m_MeanDisp_z;
    m_SimResult->m_TotalEnergy = sum_energy;

    // Path-related properties
    m_SimResult->m_TotalTime = m_TotalTime;
    m_SimResult->m_TotalHops = sum_hops;
    m_SimResult->m_NonOscHops = sum_nonosc_hops;
    m_SimResult->m_NonOscHopRatio = static_cast<double>(sum_nonosc_hops)/static_cast<double>(sum_hops);
    if (sum_nonosc_hops != 0)
    {
        m_SimResult->m_InXDirNonOscRatio = static_cast<double>(sum_in_xdir_nonosc_hops)/static_cast<double>(sum_nonosc_hops);
    }
    else m_SimResult->m_InXDirNonOscRatio = 1.0;
    if (m_ParamSet->m_PhiGradient != 0.0)
    {
        m_SimResult->m_MeanFieldContribution = avg_field_contrib * enfac;
    }
    else m_SimResult->m_MeanFieldContribution = 0.0;
    m_SimResult->m_MaxUsedPathDist = max_dist * spfac;
    m_SimResult->m_MaxUsedPathEdiff = max_ediff * enfac;
    m_SimResult->m_MinUsedStateEnergy = used_energy_range.min * enfac + enmax;
    m_SimResult->m_MaxUsedStateEnergy = used_energy_range.max * enfac + enmax;
    
    // Subtraction of reference values
    if (ref_result != nullptr)
    {
        const double density_factor = ref_result->m_EffCarriers / m_SimResult->m_EffCarriers;

        m_SimResult->m_MeanDisp -= ref_result->m_MeanDisp * density_factor;
        m_SimResult->m_MeanDisp_x -= ref_result->m_MeanDisp_x * density_factor;
        m_SimResult->m_MeanDisp_y -= ref_result->m_MeanDisp_y * density_factor;
        m_SimResult->m_MeanDisp_z -= ref_result->m_MeanDisp_z * density_factor;  
        m_SimResult->m_MeanSquaredDisp_x -= ref_result->m_MeanSquaredDisp_x * density_factor;
        m_SimResult->m_MeanSquaredDisp_y -= ref_result->m_MeanSquaredDisp_y * density_factor;
        m_SimResult->m_MeanSquaredDisp_z -= ref_result->m_MeanSquaredDisp_z * density_factor;
        m_SimResult->m_MeanDispVariance_x -= ref_result->m_MeanDispVariance_x * density_factor;
        m_SimResult->m_MeanDispVariance_y -= ref_result->m_MeanDispVariance_y * density_factor;
        m_SimResult->m_MeanDispVariance_z -= ref_result->m_MeanDispVariance_z * density_factor;
        
        m_SimResult->m_TotalTime -= ref_result->m_TotalTime;
        m_SimResult->m_TotalHops -= ref_result->m_TotalHops;
        m_SimResult->m_NonOscHops -= ref_result->m_NonOscHops;
        m_SimResult->m_NonOscHopRatio = static_cast<double>(m_SimResult->m_NonOscHops)/static_cast<double>(m_SimResult->m_TotalHops);
        if (m_SimResult->m_NonOscHops != 0)
        {
            m_SimResult->m_InXDirNonOscRatio = (static_cast<double>(sum_in_xdir_nonosc_hops) 
                / static_cast<double>(m_SimResult->m_NonOscHops))
                - ref_result->m_InXDirNonOscRatio * (static_cast<double>(ref_result->m_NonOscHops)
                / static_cast<double>(m_SimResult->m_NonOscHops));
        }
        else m_SimResult->m_InXDirNonOscRatio = 1.0;
        if (m_ParamSet->m_PhiGradient != 0.0)
        {
            m_SimResult->m_MeanFieldContribution = avg_field_contrib * enfac * (static_cast<double>(sum_hops)
                / static_cast<double>(m_SimResult->m_TotalHops))
                - ref_result->m_MeanFieldContribution * (static_cast<double>(ref_result->m_TotalHops)
                / static_cast<double>(m_SimResult->m_TotalHops));
        }
    }

    // Calculation of drift properties
    // drift_flux_density = - eff_carrier_density * drift_mobility * grad(phi)
    // drift_conductivity = - electron_charge * eff_carrier_density * drift_mobility
    // drift_flux_density = drift_conductivity/electron_charge * grad(phi)
    // units: 1/cm2s = A/Vcm * 1/C * V/cm
    // drift_flux_density = eff_carrier_density * <disp_x> / total_time
    // units: 1/cm2s = cm/(s*cm3)
    // grad(phi) has to be multiplied by energy_factor / spatial_factor (when using relative value)
    // -> drift conductivity (including scaling):
    // = <disp_x> * electron_charge * eff_carrier_density * spatial_factor / (total_time * grad(phi) * energy_factor)
    // units: A/Vcm = cm * C / (cm3 s V/cm)
    // -> drift mobility (including scaling):
    // = - drift_conductivity / (electron_charge * eff_carrier_density)
    // units: cm2/Vs = A/Vcm * cm3 / C
    if (m_ParamSet->m_PhiGradient != 0.0)
    {
        m_SimResult->m_DriftConductivity = m_SimResult->m_MeanDisp_x 
            * (Constant::echarge * m_SimResult->m_EffCarrierDensity) 
            * spfac / (1.0E14 * m_SimResult->m_TotalTime * m_GradPhiEnergy * enfac);

        m_SimResult->m_DriftMobility = - m_SimResult->m_DriftConductivity 
            / (Constant::echarge * m_SimResult->m_EffCarrierDensity);
    }
    else
    {
        m_SimResult->m_DriftConductivity = 0.0;
	    m_SimResult->m_DriftMobility = 0.0; 
    }

    // Convert conductivity unit from S/cm to S/m
    m_SimResult->m_DriftConductivity *= 100.0;
    
    // Calculation of self-diffusion properties
    // -> self-diffusion coefficient in x-direction (parallel to electric field):
    //    with field:
    //      D_x = <(disp_x - <disp_x>mob)^2> / (2 * total_time)
    //      ATTENTION: This Dx value is not exact because Nmob in general still contains 
    //                 electrons that resist the influence of the electric field
    //    without field:
    //      D_x = <disp_x^2> / (2 * total_time)
    // -> self-diffusion coefficient in yz-direction (perpendicular to electric field, in theory: <disp_y> = <disp_z> = 0):
    //    D_yz = (D_y + D_z) / 2
    //    with field and m_UseXYVariance == true:
    //      D_yz = (<(disp_y - <disp_y>mob)^2> + <(disp_z - <disp_z>mob)^2>) / (4 * total_time)
    //    without field or m_UseXYVariance == false:
    //      D_yz = (<disp_y^2> + <disp_z^2>) / (4 * total_time)
    // -> isotropic self-diffusion coefficient (without field or for weak field where Dx approx. equal to Dyz):
    //    D = (D_x + 2 * D_yz) / 3 
    // units: cm2/s
    if (m_ParamSet->m_PhiGradient != 0.0)
    {
        m_SimResult->m_DiffusionCoefficientParallel = m_SimResult->m_MeanDispVariance_x 
            / (2.0E14 * m_SimResult->m_TotalTime);
    }
    else
    {
        m_SimResult->m_DiffusionCoefficientParallel = m_SimResult->m_MeanSquaredDisp_x 
            / (2.0E14 * m_SimResult->m_TotalTime);
    }
    if ((m_ParamSet->m_PhiGradient != 0.0) && (m_ParamSet->m_UseYZVariance))
    {
        m_SimResult->m_DiffusionCoefficientTransverse = 
            (m_SimResult->m_MeanDispVariance_y + m_SimResult->m_MeanDispVariance_z) 
            / (4.0E14 * m_SimResult->m_TotalTime);
    }
    else
    {
        m_SimResult->m_DiffusionCoefficientTransverse = 
            (m_SimResult->m_MeanSquaredDisp_y + m_SimResult->m_MeanSquaredDisp_z) 
            / (4.0E14 * m_SimResult->m_TotalTime);
    }
    m_SimResult->m_DiffusionCoefficient = (m_SimResult->m_DiffusionCoefficientParallel + 
        2.0 * m_SimResult->m_DiffusionCoefficientTransverse) / 3.0;

    // Calculation of Haven ratio (self-diffusion coefficient divided by chemical diffusion coefficient)
    // Nernst-Einstein relation: Dchem_x = - k_boltz * eff_T * drift_mobility / electron_charge
    // units: cm2/s = eV/K * K * cm2/Vs / e
    // -> isotropic Haven ratio (only valid for weak electric field where Dx approx. equal to Dyz)
    //    HR = D / Dchem_x
    // -> parallel Haven ratio (parallel to electric field; for high field)
    //    HR_x = D_x / Dchem_x
    // -> transverse Haven ratio (only for weak field; D_yz = isotropic self-diffusion coefficient 
    //                            at weak field because of inaccuracy in calculation of D_x, see above)
    //    HR_yz = D_yz / Dchem_x
    if ((m_ParamSet->m_PhiGradient != 0.0) && (m_SimResult->m_DriftMobility != 0.0))
    {
        const double Dchemx = - Constant::kboltz * m_SimResult->m_EffTemp * m_SimResult->m_DriftMobility;
        m_SimResult->m_HavenRatio = m_SimResult->m_DiffusionCoefficient / Dchemx;
        m_SimResult->m_HavenRatioParallel = m_SimResult->m_DiffusionCoefficientParallel / Dchemx;
        m_SimResult->m_HavenRatioTransverse = m_SimResult->m_DiffusionCoefficientTransverse / Dchemx;
    }
    else
    {
        m_SimResult->m_HavenRatio = 0.0;
        m_SimResult->m_HavenRatioParallel = 0.0;
        m_SimResult->m_HavenRatioTransverse = 0.0;
    }
}

// Generate simulation results from simulation data
void MC::TEngine::GenerateResults()
{
    if (!m_SimulationReady)
        throw EX::TInvalidStatus("Simulation is not finished.",__func__);
    if (!m_SimResult)
        throw EX::TInvalidStatus("Result object not available.",__func__);

    const double spfac = m_DOS->GetSpatialFactor();
    const double enmax = m_ParamSet->m_MaxStateEnergy;
    const double enfac = m_DOS->GetEnergyFactor();
    
    // Analyze paths (e.g. min/max for histogram boundaries)
    GF::TMinMax<double> time_per_path;
    GF::TMinMax<double> time_per_hop;
    std::uint64_t total_paths = 0;
    std::uint64_t osc_only_paths = 0;
    std::uint64_t nonosc_only_paths = 0;
    GF::TMinMaxMean<std::uint64_t> hops_per_path;
    GF::TMinMaxMean<std::uint64_t> oscs_per_path;
    GF::TMinMaxMean<std::uint64_t> nonoscs_per_path;
    GF::TMinMaxMean<std::uint32_t> used_paths_per_state;
    GF::TMinMaxMean<std::uint64_t> hops_per_state;
    GF::TMinMaxMean<std::uint64_t> oscs_per_state;
    GF::TMinMaxMean<std::uint64_t> nonoscs_per_state;
    std::uint32_t single_osc_states = 0;
    std::uint32_t multi_osc_states = 0;
    GF::TMinMaxMean<std::uint32_t> osc_paths_per_state;
    for (const auto& state : m_Structure)
    {
        std::uint32_t t_upaths = 0;
        std::uint32_t t_opaths = 0;
        std::uint64_t t_shops = 0;
        std::uint64_t t_soscs = 0;
        std::uint64_t t_snonoscs = 0;
        for (const auto& path : state.m_Paths)
        {
            ++total_paths;
            time_per_path.check(path.m_Time);
            if (path.m_HopCount != 0)
            {
                time_per_hop.check(path.m_Time);
                hops_per_path.check(path.m_HopCount);
                ++t_upaths;
                t_shops += path.m_HopCount;
                t_soscs += path.m_OscHopCount;
                t_snonoscs += path.m_HopCount - path.m_OscHopCount;

                if (path.m_OscHopCount != 0)
                {
                    oscs_per_path.check(path.m_OscHopCount);
                    ++t_opaths;
                }

                if (path.m_HopCount != path.m_OscHopCount)
                {
                    nonoscs_per_path.check(path.m_HopCount - path.m_OscHopCount);
                    
                    if (path.m_OscHopCount == 0) ++nonosc_only_paths;
                }
                else ++osc_only_paths;
            }
        }
        if (t_upaths != 0) used_paths_per_state.check(t_upaths);
        if (t_shops != 0) hops_per_state.check(t_shops);
        if (t_soscs != 0) oscs_per_state.check(t_soscs);
        if (t_snonoscs != 0) nonoscs_per_state.check(t_snonoscs);
        if (t_opaths != 0)
        {
            osc_paths_per_state.check(t_opaths);
            if (t_opaths == 1)
                ++single_osc_states;
            else
                ++multi_osc_states;
        }
    }

    // Analyze electrons (e.g. min/max for histogram boundaries)
    GF::TMinMaxIdx<double> disp_per_electron;
    GF::TMinMax<double> dispx_per_electron;
    GF::TMinMax<double> dispy_per_electron;
    GF::TMinMax<double> dispz_per_electron;
    GF::TMinMaxIdx<std::uint64_t> hops_per_electron;
    GF::TMinMaxIdx<std::uint64_t> oscs_per_electron;
    GF::TMinMaxIdx<std::uint64_t> nonoscs_per_electron;
    std::uint32_t single_hop_electrons = 0;
    std::uint32_t single_nonosc_electrons = 0;
    std::uint32_t nonosc_only_electrons = 0;
    GF::TMinMaxMean<std::uint64_t> excess_hops_per_osc_state;
    GF::TMinMaxMean<std::uint64_t> excess_oscs_per_osc_state;
    std::vector<bool> electron_is_blocked (m_Electrons.size(),true);
    std::uint32_t blocked_mobile_electrons = 0;
    std::uint32_t blocked_zerohop_electrons = 0;
    std::uint32_t blocked_osc_electrons = 0;
    GF::TMinMax<double> next_time_per_electron;
    for (std::size_t i = 0; i < m_Electrons.size(); ++i)
    {
        const TElectron& electron = m_Electrons[i];

        for (const auto& path : m_Structure[electron.m_CurrentStateID].m_Paths)
        {
            if (m_Structure[path.m_StateID].m_ElectronID == m_ParamSet->m_StateCount)
            {
                electron_is_blocked[i] = false;
                break;
            }
        }
        if (!electron_is_blocked[i])
            next_time_per_electron.check(electron.m_MinTime);
        else if (electron.m_HopCount == 0)
            ++blocked_zerohop_electrons;
        else if (electron.m_HopCount == electron.m_OscHopCount)
            ++blocked_osc_electrons;
        else
            ++blocked_mobile_electrons;

        if (electron.m_HopCount == 0) continue;
        
        double t_electron_disp = sqrt(electron.m_Disp_x * electron.m_Disp_x + electron.m_Disp_y * 
            electron.m_Disp_y + electron.m_Disp_z * electron.m_Disp_z) * spfac;
        disp_per_electron.check(t_electron_disp,i);
        dispx_per_electron.check(electron.m_Disp_x * spfac);
        dispy_per_electron.check(electron.m_Disp_y * spfac);
        dispz_per_electron.check(electron.m_Disp_z * spfac);
        hops_per_electron.check(electron.m_HopCount,i);
        if (electron.m_OscHopCount != 0)
            oscs_per_electron.check(electron.m_OscHopCount,i);
        if (electron.m_HopCount != electron.m_OscHopCount)
            nonoscs_per_electron.check(electron.m_HopCount - electron.m_OscHopCount,i);

        if (electron.m_HopCount == 1)
            single_hop_electrons++;
        else if (electron.m_HopCount - electron.m_OscHopCount == 1) 
            single_nonosc_electrons++;
        if (electron.m_OscHopCount == 0) nonosc_only_electrons++;

        if (electron.m_HopCount == electron.m_OscHopCount)
        {
            std::uint64_t t_ophops = 0;
            std::uint64_t t_oposcs = 0;
            for (const auto& path : m_Structure[electron.m_CurrentStateID].m_Paths)
            {
                t_ophops += path.m_HopCount;
                t_oposcs += path.m_OscHopCount;
            }
            excess_hops_per_osc_state.check(t_ophops - (electron.m_HopCount/2U));
            excess_oscs_per_osc_state.check(t_oposcs - (electron.m_HopCount/2U));
        }
    }

    if (m_VL >= Verbosity::MEDIUM)
    {
        std::cout << "  Electrons with a single hop: " << single_hop_electrons 
            << " (" << 100.0*static_cast<double>(single_hop_electrons)/static_cast<double>(m_SimResult->m_MobileElectrons) << " %" << " of mobile electrons)" << std::endl;
        std::cout << "  Electrons with a single hop except osc.: " << single_nonosc_electrons 
            << " (" << 100.0*static_cast<double>(single_nonosc_electrons)/static_cast<double>(m_SimResult->m_MobileElectrons) << " %" << " of mobile electrons)" << std::endl;
        std::cout << "  Electrons with only non-osc. hops: " << nonosc_only_electrons 
            << " (" << 100.0*static_cast<double>(nonosc_only_electrons)/static_cast<double>(m_SimResult->m_MobileElectrons) << " %" << " of mobile electrons)" << std::endl;
        std::cout << "  Currently blocked mobile electrons (> 0 non-osc. hops, all paths occupied): " << blocked_mobile_electrons 
            << " (" << 100.0*static_cast<double>(blocked_mobile_electrons)/static_cast<double>(m_SimResult->m_MobileElectrons) << " %" << " of mobile electrons)" << std::endl;
        if (m_SimResult->m_ZeroHopElectrons != 0)
        {
            std::cout << "  Currently blocked stationary electrons (0 hops, all paths occupied): " << blocked_zerohop_electrons 
                << " (" << 100.0*static_cast<double>(blocked_zerohop_electrons)/static_cast<double>(m_SimResult->m_ZeroHopElectrons) << " %" << " of stationary electrons)" << std::endl;
        }
        if (m_SimResult->m_OscElectrons != 0)
        {
            std::cout << "  Currently blocked oscillating electrons (only osc. hops, all paths occupied): " << blocked_osc_electrons 
                << " (" << 100.0*static_cast<double>(blocked_osc_electrons)/static_cast<double>(m_SimResult->m_OscElectrons) << " %" << " of oscillating electrons)" << std::endl;
        }
        std::cout << "  Maximum displacement per electron (Electron-ID: " 
            << disp_per_electron.max_idx + 1 << "): " << disp_per_electron.max << " nm" << std::endl;
        std::cout << "  Maximum hops per electron (Electron-ID: " 
            << hops_per_electron.max_idx + 1 << "): " << hops_per_electron.max << std::endl;
        std::cout << "  Maximum osc. hops per electron";
        if (oscs_per_electron.has_data)
            std::cout << " (Electron-ID: " << oscs_per_electron.max_idx + 1 << "): " << oscs_per_electron.max << std::endl;
        else
            std::cout << ": 0" << std::endl;
        std::cout << "  Maximum non-osc. hops per electron";
        if (nonoscs_per_electron.has_data)
            std::cout << " (Electron-ID: " << nonoscs_per_electron.max_idx + 1 << "): " << nonoscs_per_electron.max << std::endl;
        else
            std::cout << ": 0" << std::endl;
        std::cout << "  Average additional outgoing hops from states with osc.-only electrons: ";
        if (excess_hops_per_osc_state.count != 0)
            std::cout << excess_hops_per_osc_state << std::endl;
        else
            std::cout << "0 (no osc.-only electrons)" << std::endl;
        std::cout << "  Average additional outgoing osc. hops from states with osc.-only electrons: ";
        if (excess_oscs_per_osc_state.count != 0)
            std::cout << excess_oscs_per_osc_state << std::endl;
        else
            std::cout << "0 (no osc.-only electrons)" << std::endl;

        std::cout << "  Used states: " << used_paths_per_state.count << " (" 
            << 100.0 * static_cast<double>(used_paths_per_state.count) / static_cast<double>(m_ParamSet->m_StateCount) << " %" << " of states)" << std::endl;
        std::cout << "    Average used paths per used state: " << used_paths_per_state << std::endl;
        std::cout << "    Average hops per used state: " << hops_per_state << std::endl;
        std::cout << "    Average osc. hops per used state: " << oscs_per_state << std::endl;
        std::cout << "    Average non-osc. hops per used state: " << nonoscs_per_state << std::endl;
        std::cout << "    Used states without osc.-paths: " << used_paths_per_state.count - osc_paths_per_state.count << " ("
            << 100.0 * static_cast<double>(used_paths_per_state.count - osc_paths_per_state.count)/static_cast<double>(used_paths_per_state.count) << " %" << " of used states)" << std::endl;
        std::cout << "    Used states with osc.-paths: " << osc_paths_per_state.count << " ("
            << 100.0 * static_cast<double>(osc_paths_per_state.count)/static_cast<double>(used_paths_per_state.count) << " %" << " of used states)" << std::endl;
        if (osc_paths_per_state.count != 0)
        {
            std::cout << "      States with one osc.-path: " << single_osc_states << std::endl;
            std::cout << "      States with more than one osc.-path: " << multi_osc_states << std::endl;
            std::cout << "      Average osc.-paths of states with osc.-paths: " << osc_paths_per_state << std::endl;
        }

        std::cout << "  Used paths: " << hops_per_path.count << " (" 
            << 100.0 * static_cast<double>(hops_per_path.count) / static_cast<double>(total_paths) << " %" << " of paths)" << std::endl;
        std::cout << "    Paths with osc. hops: " << oscs_per_path.count << " (" 
            << 100.0 * static_cast<double>(oscs_per_path.count) / static_cast<double>(hops_per_path.count) << " %" << " of used paths)" << std::endl;
        std::cout << "      Paths with only osc. hops: " << osc_only_paths << " (" 
            << 100.0 * static_cast<double>(osc_only_paths) / static_cast<double>(hops_per_path.count) << " %" << " of used paths)" << std::endl;
        std::cout << "    Paths with non-osc. hops: " << nonoscs_per_path.count << " (" 
            << 100.0 * static_cast<double>(nonoscs_per_path.count) / static_cast<double>(hops_per_path.count) << " %" << " of used paths)" << std::endl;
        std::cout << "      Paths with only non-osc. hops: " << nonosc_only_paths << " (" 
            << 100.0 * static_cast<double>(nonosc_only_paths) / static_cast<double>(hops_per_path.count) << " %" << " of used paths)" << std::endl;
        std::cout << "    Paths with both osc. and non-osc. hops: " << hops_per_path.count - osc_only_paths - nonosc_only_paths << " (" 
            << 100.0 * static_cast<double>(hops_per_path.count - osc_only_paths - nonosc_only_paths) / static_cast<double>(hops_per_path.count) << " %" << " of used paths)" << std::endl;
        std::cout << "    Average hops per used path: " << hops_per_path << std::endl;
        std::cout << "    Average osc. hops of paths with osc. hops: ";
        if (oscs_per_path.count != 0)
            std::cout << oscs_per_path << std::endl;
        else
            std::cout << "0 (no osc. hops)" << std::endl;
        std::cout << "    Average non-osc. hops of paths with non-osc. hops: ";
        if (nonoscs_per_path.count != 0)
            std::cout << nonoscs_per_path << std::endl;
        else
            std::cout << "0 (no non-osc. hops)" << std::endl;
        
        std::cout << "Collecting histogram data: " << std::flush;
    }

    // -- CONFIGURE HISTOGRAMS --
	m_SimResult->m_HStateEnergy.ConfigureLinear("Energy","eV",m_ParamSet->m_MinStateEnergy,
        m_ParamSet->m_MaxStateEnergy,50);
    m_SimResult->m_HStateEnergy.AddProperty(HP::ST, HT::COUNT, "States");
    m_SimResult->m_HStateEnergy.AddProperty(HP::USED_ST, HT::COUNT, "UsedStates");
    m_SimResult->m_HStateEnergy.AddProperty(HP::EL, HT::COUNT, "Electrons", "", 
        "energetic distribution of electrons refers to the end of the simulation");
    m_SimResult->m_HStateEnergy.AddProperty(HP::OCC_R, HT::AVG, "OccRatio");
    m_SimResult->m_HStateEnergy.AddProperty(HP::OCC_FIT, HT::VALUE, "OccFit");
    m_SimResult->m_HStateEnergy.AddProperty(HP::MOB_EL, HT::COUNT, "MobileElectrons");
    m_SimResult->m_HStateEnergy.AddProperty(HP::MOB_R, HT::AVG, "MobileRatio");
    m_SimResult->m_HStateEnergy.AddProperty(HP::ZEROHOP_EL, HT::COUNT, "ZeroHopElectrons");
    m_SimResult->m_HStateEnergy.AddProperty(HP::OSC_EL, HT::COUNT, "OscElectrons");
    m_SimResult->m_HStateEnergy.AddProperty(HP::ZERODISP_EL, HT::COUNT, "ZeroDispElectrons", "",
        "zero-hop electrons and electrons with only osc. are excluded from zero-disp. electrons");
    m_SimResult->m_HStateEnergy.AddProperty(HP::HOPS_PER_EL, HT::AVG_STDDEV, "Hops/Electron");
    m_SimResult->m_HStateEnergy.AddProperty(HP::OSC_PER_EL, HT::AVG_STDDEV, "Osc/Electron");
    m_SimResult->m_HStateEnergy.AddProperty(HP::NONOSC_PER_EL, HT::AVG_STDDEV, "NonOsc/Electron", "",
        "averages per electron refer to mobile electrons");
    if (m_ParamSet->m_EqHopLimit != 0)
    {
        m_SimResult->m_HStateEnergy.AddProperty(HP::NEWMOB_EL, HT::COUNT, "MobileElectrons(**)", "",
            "** = electrons that had their first hop after the equilibration");
        m_SimResult->m_HStateEnergy.AddProperty(HP::NEWOSC_EL, HT::COUNT, "OscElectrons(**)");
    }
    m_SimResult->m_HStateEnergy.AddProperty(HP::OUT_HOPS, HT::COUNT, "OutgoingHops");
    m_SimResult->m_HStateEnergy.AddProperty(HP::OUT_OSC, HT::COUNT, "OutgoingOscHops");
    m_SimResult->m_HStateEnergy.AddProperty(HP::OUT_NONOSC, HT::COUNT, "OutgoingNonOscHops");
    m_SimResult->m_HStateEnergy.AddProperty(HP::IN_HOPS, HT::COUNT, "IncomingHops");
    m_SimResult->m_HStateEnergy.AddProperty(HP::IN_OSC, HT::COUNT, "IncomingOscHops");
    m_SimResult->m_HStateEnergy.AddProperty(HP::IN_NONOSC, HT::COUNT, "IncomingNonOscHops");
    m_SimResult->m_HStateEnergy.AddProperty(HP::OCC_RTIME, HT::AVG_STDDEV, "RelOccTime/State");
    m_SimResult->m_HStateEnergy.AddProperty(HP::EFF_EL_CALC, HT::SUM, "EffCarriers(calc)");
    m_SimResult->m_HStateEnergy.AddProperty(HP::EFF_EL_A, HT::COUNT, "EffCarriers(a)","",
        "a = effective carriers identified from occupied states with maximum outgoing hops");
    m_SimResult->m_HStateEnergy.AddProperty(HP::EFF_EL_B, HT::COUNT, "EffCarriers(b)","",
        "b = effective carriers identified from occupied states with maximum outgoing non-oscillating hops");
    m_SimResult->m_HStateEnergy.AddProperty(HP::EFF_EL_C, HT::COUNT, "EffCarriers(c)","",
        "c = effective carriers identified from occupied states with minimal difference between occupied and unoccupied timespan");


    m_SimResult->m_HPathTime.ConfigureLogarithmic("Time","s",time_per_path.min,time_per_path.max,50);
    m_SimResult->m_HPathTime.AddProperty(HP::PT, HT::COUNT, "Paths");
    m_SimResult->m_HPathTime.AddProperty(HP::USED_PT, HT::COUNT, "UsedPaths");
    m_SimResult->m_HPathTime.AddProperty(HP::USED_R, HT::AVG, "UsedRatio");
    m_SimResult->m_HPathTime.AddProperty(HP::HOPS, HT::COUNT, "Hops");
    m_SimResult->m_HPathTime.AddProperty(HP::OSC, HT::COUNT, "OscHops");
    m_SimResult->m_HPathTime.AddProperty(HP::NONOSC, HT::COUNT, "NonOscHops");
    m_SimResult->m_HPathTime.AddProperty(HP::OSC_R, HT::AVG, "OscRatio");
    m_SimResult->m_HPathTime.AddProperty(HP::HOPS_PER_PT, HT::AVG_STDDEV, "Hops/Path", "",
        "unused paths are excluded from average hops and osc. hops");
    m_SimResult->m_HPathTime.AddProperty(HP::OSC_PER_PT, HT::AVG_STDDEV, "OscHops/Path");
    m_SimResult->m_HPathTime.AddProperty(HP::DIST_PER_PT, HT::AVG_STDDEV, "Dist/Path", "nm");
    m_SimResult->m_HPathTime.AddProperty(HP::EDIFF_PER_PT, HT::AVG_STDDEV, "Ediff/Path", "eV",
        "energy difference per path includes the contribution from electric field");


    m_SimResult->m_HHopTime.ConfigureLogarithmic("Time","s",time_per_hop.min,time_per_hop.max,50);
    m_SimResult->m_HHopTime.AddProperty(HP::HOPS, HT::COUNT, "Hops");
    m_SimResult->m_HHopTime.AddProperty(HP::OSC, HT::COUNT, "OscHops");
    m_SimResult->m_HHopTime.AddProperty(HP::NONOSC, HT::COUNT, "NonOscHops");
    m_SimResult->m_HHopTime.AddProperty(HP::OSC_R, HT::AVG, "OscRatio");
    m_SimResult->m_HHopTime.AddProperty(HP::HOPS_PER_PT, HT::AVG_STDDEV, "Hops/Path", "",
        "unused paths are excluded from averages per path");
    m_SimResult->m_HHopTime.AddProperty(HP::OSC_PER_PT, HT::AVG_STDDEV, "OscHops/Path");
    m_SimResult->m_HHopTime.AddProperty(HP::PT, HT::COUNT, "Paths");
    m_SimResult->m_HHopTime.AddProperty(HP::USED_PT, HT::COUNT, "UsedPaths");
    m_SimResult->m_HHopTime.AddProperty(HP::USED_R, HT::AVG, "UsedRatio");
    m_SimResult->m_HHopTime.AddProperty(HP::DIST_PER_PT, HT::AVG_STDDEV, "Dist/Path", "nm");
    m_SimResult->m_HHopTime.AddProperty(HP::EDIFF_PER_PT, HT::AVG_STDDEV, "Ediff/Path", "eV",
        "energy difference per path includes the contribution from electric field");


	m_SimResult->m_HDisp.ConfigureLinear("Disp","nm",disp_per_electron.min,disp_per_electron.max,50);
    m_SimResult->m_HDisp.AddProperty(HP::EL, HT::COUNT, "Electrons","",
        "\"electron\" refers to mobile electrons (> 0 non-osc. hops)");
    m_SimResult->m_HDisp.AddProperty(HP::INCOSC_EL, HT::COUNT, "Electrons(*)", "",
        "* = including oscillating electrons");
    if (m_ParamSet->m_EqHopLimit != 0)
    {
        m_SimResult->m_HDisp.AddProperty(HP::NEW_EL, HT::COUNT, "Electrons(**)", "",
            "** = electrons that had their first hop after the equilibration");
        m_SimResult->m_HDisp.AddProperty(HP::NEW_INCOSC_EL, HT::COUNT, "Electrons(*,**)");
    }
    m_SimResult->m_HDisp.AddProperty(HP::HOPS_PER_EL, HT::AVG_STDDEV, "Hops/Electron");
    m_SimResult->m_HDisp.AddProperty(HP::OSC_PER_EL, HT::AVG_STDDEV, "Osc/Electron");
    m_SimResult->m_HDisp.AddProperty(HP::NONOSC_PER_EL, HT::AVG_STDDEV, "NonOsc/Electron");


    m_SimResult->m_HDisp_x.ConfigureLinear("xDisp","nm",dispx_per_electron.min,dispx_per_electron.max,50);
    m_SimResult->m_HDisp_x.AddProperty(HP::EL, HT::COUNT, "Electrons","",
        "\"electron\" refers to mobile electrons (> 0 non-osc. hops)");
    m_SimResult->m_HDisp_x.AddProperty(HP::INCOSC_EL, HT::COUNT, "Electrons(*)", "",
        "* = including oscillating electrons");
    if (m_ParamSet->m_EqHopLimit != 0)
    {
        m_SimResult->m_HDisp_x.AddProperty(HP::NEW_EL, HT::COUNT, "Electrons(**)", "",
            "** = electrons that had their first hop after the equilibration");
        m_SimResult->m_HDisp_x.AddProperty(HP::NEW_INCOSC_EL, HT::COUNT, "Electrons(*,**)");
    }
    m_SimResult->m_HDisp_x.AddProperty(HP::HOPS_PER_EL, HT::AVG_STDDEV, "Hops/Electron");
    m_SimResult->m_HDisp_x.AddProperty(HP::OSC_PER_EL, HT::AVG_STDDEV, "Osc/Electron");
    m_SimResult->m_HDisp_x.AddProperty(HP::NONOSC_PER_EL, HT::AVG_STDDEV, "NonOsc/Electron");


    m_SimResult->m_HDisp_y.ConfigureLinear("yDisp","nm",dispy_per_electron.min,dispy_per_electron.max,50);
    m_SimResult->m_HDisp_y.AddProperty(HP::EL, HT::COUNT, "Electrons","",
        "\"electron\" refers to mobile electrons (> 0 non-osc. hops)");
    m_SimResult->m_HDisp_y.AddProperty(HP::INCOSC_EL, HT::COUNT, "Electrons(*)", "",
        "* = including oscillating electrons");
    if (m_ParamSet->m_EqHopLimit != 0)
    {
        m_SimResult->m_HDisp_y.AddProperty(HP::NEW_EL, HT::COUNT, "Electrons(**)", "",
            "** = electrons that had their first hop after the equilibration");
        m_SimResult->m_HDisp_y.AddProperty(HP::NEW_INCOSC_EL, HT::COUNT, "Electrons(*,**)");
    }
    m_SimResult->m_HDisp_y.AddProperty(HP::HOPS_PER_EL, HT::AVG_STDDEV, "Hops/Electron");
    m_SimResult->m_HDisp_y.AddProperty(HP::OSC_PER_EL, HT::AVG_STDDEV, "Osc/Electron");
    m_SimResult->m_HDisp_y.AddProperty(HP::NONOSC_PER_EL, HT::AVG_STDDEV, "NonOsc/Electron");


    m_SimResult->m_HDisp_z.ConfigureLinear("zDisp","nm",dispz_per_electron.min,dispz_per_electron.max,50);
    m_SimResult->m_HDisp_z.AddProperty(HP::EL, HT::COUNT, "Electrons","",
        "\"electron\" refers to mobile electrons (> 0 non-osc. hops)");
    m_SimResult->m_HDisp_z.AddProperty(HP::INCOSC_EL, HT::COUNT, "Electrons(*)", "",
        "* = including oscillating electrons");
    if (m_ParamSet->m_EqHopLimit != 0)
    {
        m_SimResult->m_HDisp_z.AddProperty(HP::NEW_EL, HT::COUNT, "Electrons(**)", "",
            "** = electrons that had their first hop after the equilibration");
        m_SimResult->m_HDisp_z.AddProperty(HP::NEW_INCOSC_EL, HT::COUNT, "Electrons(*,**)");
    }
    m_SimResult->m_HDisp_z.AddProperty(HP::HOPS_PER_EL, HT::AVG_STDDEV, "Hops/Electron");
    m_SimResult->m_HDisp_z.AddProperty(HP::OSC_PER_EL, HT::AVG_STDDEV, "Osc/Electron");
    m_SimResult->m_HDisp_z.AddProperty(HP::NONOSC_PER_EL, HT::AVG_STDDEV, "NonOsc/Electron");
	

	m_SimResult->m_HHopCount.ConfigureLogarithmic("Hops","",
        static_cast<double>(std::min(hops_per_electron.min,hops_per_state.min)),
        static_cast<double>(std::max(hops_per_electron.max,hops_per_state.max)),50);
    m_SimResult->m_HHopCount.AddProperty(HP::EL, HT::COUNT, "Electrons","",
        "\"electron\" refers to mobile electrons (> 0 non-osc. hops)");
    m_SimResult->m_HHopCount.AddProperty(HP::INCOSC_EL, HT::COUNT, "Electrons(*)", "",
        "* = including oscillating electrons");
    m_SimResult->m_HHopCount.AddProperty(HP::OSC_EL, HT::COUNT, "OscElectrons");
    if (m_ParamSet->m_EqHopLimit != 0)
    {
        m_SimResult->m_HHopCount.AddProperty(HP::NEW_EL, HT::COUNT, "Electrons(**)", "",
            "** = electrons that had their first hop after the equilibration");
        m_SimResult->m_HHopCount.AddProperty(HP::NEW_INCOSC_EL, HT::COUNT, "Electrons(*,**)");
        m_SimResult->m_HHopCount.AddProperty(HP::NEWOSC_EL, HT::COUNT, "OscElectrons(**)");
    }
    m_SimResult->m_HHopCount.AddProperty(HP::ZERODISP_EL, HT::COUNT, "ZeroDispElectrons", "",
        "zero-hop electrons and electrons with only osc. are excluded from zero-disp. electrons");
    m_SimResult->m_HHopCount.AddProperty(HP::OSC_R, HT::AVG_STDDEV, "OscRatio");
    m_SimResult->m_HHopCount.AddProperty(HP::DISP_PER_EL, HT::AVG_STDDEV, "Disp/Electron", "nm");
    m_SimResult->m_HHopCount.AddProperty(HP::DISP_PER_INCOSC_EL, HT::AVG_STDDEV, "Disp/Electron(*)", "nm");
    m_SimResult->m_HHopCount.AddProperty(HP::ST, HT::COUNT, "States", "", 
        "for states the outgoing hops are meant");
    m_SimResult->m_HHopCount.AddProperty(HP::OCC_ST, HT::COUNT, "OccupiedStates");
    m_SimResult->m_HHopCount.AddProperty(HP::OSC_ST, HT::COUNT, "OscOccStates");

    if ((oscs_per_electron.has_data) && (oscs_per_state.count != 0))
    {
        m_SimResult->m_HOscHopCount.ConfigureLogarithmic("OscHops","",
            static_cast<double>(std::min(oscs_per_electron.min,oscs_per_state.min)),
            static_cast<double>(std::max(oscs_per_electron.max,oscs_per_state.max)),50);
        m_SimResult->m_HOscHopCount.AddProperty(HP::EL, HT::COUNT, "Electrons","",
            "\"electron\" refers to mobile electrons (> 0 non-osc. hops)");
        m_SimResult->m_HOscHopCount.AddProperty(HP::INCOSC_EL, HT::COUNT, "Electrons(*)", "",
            "* = including oscillating electrons");
        m_SimResult->m_HOscHopCount.AddProperty(HP::OSC_EL, HT::COUNT, "OscElectrons");
        if (m_ParamSet->m_EqHopLimit != 0)
        {
            m_SimResult->m_HOscHopCount.AddProperty(HP::NEW_EL, HT::COUNT, "Electrons(**)", "",
                "** = electrons that had their first hop after the equilibration");
            m_SimResult->m_HOscHopCount.AddProperty(HP::NEW_INCOSC_EL, HT::COUNT, "Electrons(*,**)");
            m_SimResult->m_HOscHopCount.AddProperty(HP::NEWOSC_EL, HT::COUNT, "OscElectrons(**)");
        }
        m_SimResult->m_HOscHopCount.AddProperty(HP::ZERODISP_EL, HT::COUNT, "ZeroDispElectrons", "",
            "zero-hop electrons and electrons with only osc. are excluded from zero-disp. electrons");
        m_SimResult->m_HOscHopCount.AddProperty(HP::OSC_R, HT::AVG_STDDEV, "OscRatio");
        m_SimResult->m_HOscHopCount.AddProperty(HP::DISP_PER_EL, HT::AVG_STDDEV, "Disp/Electron", "nm");
        m_SimResult->m_HOscHopCount.AddProperty(HP::DISP_PER_INCOSC_EL, HT::AVG_STDDEV, "Disp/Electron(*)", "nm");
        m_SimResult->m_HOscHopCount.AddProperty(HP::ST, HT::COUNT, "States", "", 
            "for states the outgoing hops are meant");
        m_SimResult->m_HOscHopCount.AddProperty(HP::OCC_ST, HT::COUNT, "OccupiedStates");
        m_SimResult->m_HOscHopCount.AddProperty(HP::OSC_ST, HT::COUNT, "OscOccStates");
    }
    else m_SimResult->m_HOscHopCount = THistogram();


    if ((nonoscs_per_electron.has_data) && (nonoscs_per_state.count != 0))
    {
        m_SimResult->m_HNonOscHopCount.ConfigureLogarithmic("NonOscHops","",
            static_cast<double>(std::min(nonoscs_per_electron.min,nonoscs_per_state.min)),
            static_cast<double>(std::max(nonoscs_per_electron.max,nonoscs_per_state.max)),50);
        m_SimResult->m_HNonOscHopCount.AddProperty(HP::EL, HT::COUNT, "Electrons","",
            "\"electron\" refers to mobile electrons (> 0 non-osc. hops)");
        if (m_ParamSet->m_EqHopLimit != 0)
        {
            m_SimResult->m_HNonOscHopCount.AddProperty(HP::NEW_EL, HT::COUNT, "Electrons(**)", "",
                "** = electrons that had their first hop after the equilibration");
        }
        m_SimResult->m_HNonOscHopCount.AddProperty(HP::ZERODISP_EL, HT::COUNT, "ZeroDispElectrons");
        m_SimResult->m_HNonOscHopCount.AddProperty(HP::NONOSC_R, HT::AVG_STDDEV, "NonOscRatio");
        m_SimResult->m_HNonOscHopCount.AddProperty(HP::DISP_PER_EL, HT::AVG_STDDEV, "Disp/Electron", "nm");
        m_SimResult->m_HNonOscHopCount.AddProperty(HP::ST, HT::COUNT, "States", "", 
            "for states the outgoing hops are meant");
        m_SimResult->m_HNonOscHopCount.AddProperty(HP::OCC_ST, HT::COUNT, "OccupiedStates");
        m_SimResult->m_HNonOscHopCount.AddProperty(HP::OSC_ST, HT::COUNT, "OscOccStates");
    }
    else m_SimResult->m_HNonOscHopCount = THistogram();


    m_SimResult->m_HStateEnergyDifference.ConfigureLinear("Ediff","eV",
        -m_SimResult->m_MaxPathEdiff,m_SimResult->m_MaxPathEdiff,50);
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::HOPS, HT::COUNT, "Hops");
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::OSC, HT::COUNT, "OscHops");
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::NONOSC, HT::COUNT, "NonOscHops");
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::OSC_R, HT::AVG, "OscRatio");
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::HOPS_PER_PT, HT::AVG_STDDEV, "Hops/Path", "",
        "unused paths are excluded from averages per path");
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::OSC_PER_PT, HT::AVG_STDDEV, "OscHops/Path");
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::PT, HT::COUNT, "Paths");
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::USED_PT, HT::COUNT, "UsedPaths");
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::USED_R, HT::AVG, "UsedRatio");
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::DIST_PER_HOP, HT::AVG_STDDEV, "Dist/Hop", "nm");
    m_SimResult->m_HStateEnergyDifference.AddProperty(HP::DIST_PER_NONOSC, HT::AVG_STDDEV, "Dist/NonOscHop", "nm");


    if (m_ParamSet->m_PhiGradient != 0.0)
    {
        double max_field_contrib = fabs(m_ParamSet->m_PhiGradient * m_SimResult->m_MaxPathDist * 1.0E-7);
        m_SimResult->m_HFieldEnergyContribution.ConfigureLinear("Efield","eV",-max_field_contrib,
            max_field_contrib,50);
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::HOPS, HT::COUNT, "Hops");
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::OSC, HT::COUNT, "OscHops");
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::NONOSC, HT::COUNT, "NonOscHops");
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::OSC_R, HT::AVG, "OscRatio");
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::HOPS_PER_PT, HT::AVG_STDDEV, "Hops/Path", "",
            "unused paths are excluded from averages per path");
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::OSC_PER_PT, HT::AVG_STDDEV, "OscHops/Path");
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::PT, HT::COUNT, "Paths");
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::USED_PT, HT::COUNT, "UsedPaths");
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::USED_R, HT::AVG, "UsedRatio");
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::DIST_PER_HOP, HT::AVG_STDDEV, "Dist/Hop", "nm");
        m_SimResult->m_HFieldEnergyContribution.AddProperty(HP::DIST_PER_NONOSC, HT::AVG_STDDEV, "Dist/NonOscHop", "nm");
    }
    else m_SimResult->m_HFieldEnergyContribution = THistogram();


    m_SimResult->m_HDistance.ConfigureLinear("Distance","nm",0.0,m_SimResult->m_MaxPathDist,50);
    m_SimResult->m_HDistance.AddProperty(HP::HOPS, HT::COUNT, "Hops");
    m_SimResult->m_HDistance.AddProperty(HP::OSC, HT::COUNT, "OscHops");
    m_SimResult->m_HDistance.AddProperty(HP::NONOSC, HT::COUNT, "NonOscHops");
    m_SimResult->m_HDistance.AddProperty(HP::OSC_R, HT::AVG, "OscRatio");
    m_SimResult->m_HDistance.AddProperty(HP::HOPS_PER_PT, HT::AVG_STDDEV, "Hops/Path", "",
        "unused paths are excluded from averages per path");
    m_SimResult->m_HDistance.AddProperty(HP::OSC_PER_PT, HT::AVG_STDDEV, "OscHops/Path");
    m_SimResult->m_HDistance.AddProperty(HP::PT, HT::COUNT, "Paths");
    m_SimResult->m_HDistance.AddProperty(HP::USED_PT, HT::COUNT, "UsedPaths");
    m_SimResult->m_HDistance.AddProperty(HP::USED_R, HT::AVG, "UsedRatio");
    m_SimResult->m_HDistance.AddProperty(HP::EDIFF_PER_HOP, HT::AVG_STDDEV, "Ediff/Hop", "eV",
        "energy difference per hop includes the contribution from electric field");
    m_SimResult->m_HDistance.AddProperty(HP::EDIFF_PER_NONOSC, HT::AVG_STDDEV, "Ediff/NonOscHop", "eV");


    m_SimResult->m_HNextTime.ConfigureLogarithmic("Time","s",next_time_per_electron.min,
        next_time_per_electron.max,50);
    m_SimResult->m_HNextTime.AddProperty(HP::EL, HT::COUNT, "Electrons");
    m_SimResult->m_HNextTime.AddProperty(HP::MOB_EL, HT::COUNT, "MobileElectrons");
    m_SimResult->m_HNextTime.AddProperty(HP::ZEROHOP_EL, HT::COUNT, "ZeroHopElectrons");
    m_SimResult->m_HNextTime.AddProperty(HP::OSC_EL, HT::COUNT, "OscElectrons");


    m_SimResult->m_HRelOccTime.ConfigureLinear("RelOcc","",0.0,1.0,50);
    m_SimResult->m_HRelOccTime.AddProperty(HP::ST, HT::COUNT, "States");
    m_SimResult->m_HRelOccTime.AddProperty(HP::OCC_ST, HT::COUNT, "OccupiedStates");
    m_SimResult->m_HRelOccTime.AddProperty(HP::OSC_ST, HT::COUNT, "OscOccStates");
    m_SimResult->m_HRelOccTime.AddProperty(HP::ST_EGY, HT::AVG_STDDEV, "Energy/State", "eV");
    m_SimResult->m_HRelOccTime.AddProperty(HP::OUT_HOPS, HT::AVG_STDDEV, "OutHops/State");
    m_SimResult->m_HRelOccTime.AddProperty(HP::OUT_OSC, HT::AVG_STDDEV, "OutOscHops/State");
    m_SimResult->m_HRelOccTime.AddProperty(HP::OUT_NONOSC, HT::AVG_STDDEV, "OutNonOscHops/State");
    
    
    // -- FILL HISTOGRAMS --
    
    // Additional statistics on energies and distances
    std::uint64_t total_hops = 0;
    std::vector<std::uint64_t> state_hops (m_Structure.size(), 0);
    std::vector<std::uint64_t> state_oscs (m_Structure.size(), 0);
    std::vector<std::uint64_t> state_nonoscs (m_Structure.size(), 0);
    GF::TOuterCounter<std::uint64_t> outer_distance_hops (0.0, 0.9*m_SimResult->m_MaxPathDist);
    GF::TOuterCounter<std::uint64_t> outer_ediff_hops (-0.9*m_SimResult->m_MaxPathEdiff,0.9*m_SimResult->m_MaxPathEdiff);
    GF::TOuterCounter<std::uint64_t> outer_energy_hops (-0.95, -0.05);
    GF::TOuterCounter<std::uint32_t> outer_electrons (-0.95, -0.05);
    GF::TOuterCounter<std::uint32_t> outer_mobile_electrons (-0.95, -0.05);
    std::vector<bool> state_is_used (m_Structure.size(), false);
    GF::TOuterCounter<std::uint32_t> outer_states (-0.95, -0.05);
    GF::TOuterCounter<std::uint32_t> outer_used_states (-0.95, -0.05);
    GF::TMinMax<double> dist_per_used_path;
    GF::TMinMax<double> ediff_per_used_path;

    // Paths
    for (std::size_t i = 0; i < m_Structure.size(); ++i)
    {
        const TLocalState& state = m_Structure[i];
        for (const auto& path : state.m_Paths)
        {
            const double& old_x {state.m_Pos_x};
            double new_x {m_Structure[path.m_StateID].m_Pos_x};
            if (old_x - new_x > 0.5) new_x += 1.0;
            if (new_x - old_x > 0.5) new_x -= 1.0;
            const double& old_y {state.m_Pos_y};
            double new_y {m_Structure[path.m_StateID].m_Pos_y};
            if (old_y - new_y > 0.5) new_y += 1.0;
            if (new_y - old_y > 0.5) new_y -= 1.0;
            const double& old_z {state.m_Pos_z};
            double new_z {m_Structure[path.m_StateID].m_Pos_z};
            if (old_z - new_z > 0.5) new_z += 1.0;
            if (new_z - old_z > 0.5) new_z -= 1.0;
            double t_path_dist = sqrt((new_x - old_x)*(new_x - old_x) + (new_y - old_y)*
                (new_y - old_y) + (new_z - old_z)*(new_z - old_z)) * spfac;

            double t_path_ediff = (m_Structure[path.m_StateID].m_Energy - state.m_Energy) * enfac;

            double t_field_contrib = 0.0;
            if (m_ParamSet->m_PhiGradient != 0.0)
            {
                t_field_contrib = - (new_x - old_x) * m_GradPhiEnergy * enfac;
            }

            // Any path
            m_SimResult->m_HPathTime.IncCount(HP::PT, path.m_Time);
            m_SimResult->m_HHopTime.IncCount(HP::PT, path.m_Time);
            m_SimResult->m_HPathTime.AddValue(HP::DIST_PER_PT, path.m_Time, t_path_dist);
            m_SimResult->m_HPathTime.AddValue(HP::EDIFF_PER_PT, path.m_Time, t_path_ediff + t_field_contrib);
            //
            m_SimResult->m_HStateEnergyDifference.IncCount(HP::PT, t_path_ediff);
            //
            if (m_SimResult->m_HFieldEnergyContribution.m_IsConfigured)
            {
                m_SimResult->m_HFieldEnergyContribution.IncCount(HP::PT, t_field_contrib);
            }
            //
            m_SimResult->m_HDistance.IncCount(HP::PT, t_path_dist);

            // Unused path (0 hops)
            if (path.m_HopCount == 0) 
            {
                m_SimResult->m_HPathTime.AddValue(HP::USED_R, path.m_Time, 0.0);
                m_SimResult->m_HHopTime.AddValue(HP::USED_R, path.m_Time, 0.0);
                //
                m_SimResult->m_HStateEnergyDifference.AddValue(HP::USED_R, t_path_ediff, 0.0);
                //
                if (m_SimResult->m_HFieldEnergyContribution.m_IsConfigured)
                {
                    m_SimResult->m_HFieldEnergyContribution.AddValue(HP::USED_R, t_field_contrib, 0.0);
                }
                //
                m_SimResult->m_HDistance.AddValue(HP::USED_R, t_path_dist, 0.0);

                continue;
            }
            
            // Used path (> 0 hops)
            state_is_used[i] = true;
            state_is_used[path.m_StateID] = true;
            total_hops += path.m_HopCount;
            state_hops[i] += path.m_HopCount;
            state_oscs[i] += path.m_OscHopCount;
            state_nonoscs[i] += path.m_HopCount - path.m_OscHopCount;
            outer_energy_hops.check(state.m_Energy, path.m_HopCount);
            outer_distance_hops.check(t_path_dist, path.m_HopCount);
            outer_ediff_hops.check(t_path_ediff, path.m_HopCount);
            dist_per_used_path.check(t_path_dist);
            ediff_per_used_path.check(t_path_ediff);
            //
            m_SimResult->m_HStateEnergy.IncCount(HP::OUT_HOPS, state.m_Energy * enfac + enmax, path.m_HopCount);
            m_SimResult->m_HStateEnergy.IncCount(HP::OUT_OSC, state.m_Energy * enfac + enmax, path.m_OscHopCount);
            m_SimResult->m_HStateEnergy.IncCount(HP::OUT_NONOSC, state.m_Energy * enfac + enmax, path.m_HopCount - path.m_OscHopCount);
            m_SimResult->m_HStateEnergy.IncCount(HP::IN_HOPS, m_Structure[path.m_StateID].m_Energy * enfac + enmax, path.m_HopCount);
            m_SimResult->m_HStateEnergy.IncCount(HP::IN_OSC, m_Structure[path.m_StateID].m_Energy * enfac + enmax, path.m_OscHopCount);
            m_SimResult->m_HStateEnergy.IncCount(HP::IN_NONOSC, m_Structure[path.m_StateID].m_Energy * enfac + enmax, path.m_HopCount - path.m_OscHopCount);
            //
            m_SimResult->m_HPathTime.IncCount(HP::USED_PT, path.m_Time);
            m_SimResult->m_HHopTime.IncCount(HP::USED_PT, path.m_Time);
            m_SimResult->m_HPathTime.AddValue(HP::USED_R, path.m_Time, 1.0);
            m_SimResult->m_HHopTime.AddValue(HP::USED_R, path.m_Time, 1.0);
            m_SimResult->m_HPathTime.IncCount(HP::HOPS, path.m_Time, path.m_HopCount);
            m_SimResult->m_HHopTime.IncCount(HP::HOPS, path.m_Time, path.m_HopCount);
            m_SimResult->m_HPathTime.IncCount(HP::OSC, path.m_Time, path.m_OscHopCount);
            m_SimResult->m_HHopTime.IncCount(HP::OSC, path.m_Time, path.m_OscHopCount);
            m_SimResult->m_HPathTime.IncCount(HP::NONOSC, path.m_Time, path.m_HopCount - path.m_OscHopCount);
            m_SimResult->m_HHopTime.IncCount(HP::NONOSC, path.m_Time, path.m_HopCount - path.m_OscHopCount);
            m_SimResult->m_HPathTime.AddValue(HP::OSC_R, path.m_Time, 1.0, path.m_OscHopCount);
            m_SimResult->m_HPathTime.AddValue(HP::OSC_R, path.m_Time, 0.0, path.m_HopCount - path.m_OscHopCount);
            m_SimResult->m_HHopTime.AddValue(HP::OSC_R, path.m_Time, 1.0, path.m_OscHopCount);
            m_SimResult->m_HHopTime.AddValue(HP::OSC_R, path.m_Time, 0.0, path.m_HopCount - path.m_OscHopCount);
            m_SimResult->m_HPathTime.AddValue(HP::HOPS_PER_PT, path.m_Time, static_cast<double>(path.m_HopCount));
            m_SimResult->m_HHopTime.AddValue(HP::HOPS_PER_PT, path.m_Time, static_cast<double>(path.m_HopCount));
            m_SimResult->m_HPathTime.AddValue(HP::OSC_PER_PT, path.m_Time, static_cast<double>(path.m_OscHopCount));
            m_SimResult->m_HHopTime.AddValue(HP::OSC_PER_PT, path.m_Time, static_cast<double>(path.m_OscHopCount));       
            m_SimResult->m_HHopTime.AddValue(HP::DIST_PER_PT, path.m_Time, t_path_dist);
            m_SimResult->m_HHopTime.AddValue(HP::EDIFF_PER_PT, path.m_Time, t_path_ediff + t_field_contrib);
            //
            m_SimResult->m_HStateEnergyDifference.IncCount(HP::USED_PT, t_path_ediff);
            m_SimResult->m_HStateEnergyDifference.AddValue(HP::USED_R, t_path_ediff, 1.0);
            m_SimResult->m_HStateEnergyDifference.IncCount(HP::HOPS, t_path_ediff, path.m_HopCount);
            m_SimResult->m_HStateEnergyDifference.IncCount(HP::OSC, t_path_ediff, path.m_OscHopCount);
            m_SimResult->m_HStateEnergyDifference.IncCount(HP::NONOSC, t_path_ediff, path.m_HopCount - path.m_OscHopCount);
            m_SimResult->m_HStateEnergyDifference.AddValue(HP::OSC_R, t_path_ediff, 1.0, path.m_OscHopCount);
            m_SimResult->m_HStateEnergyDifference.AddValue(HP::OSC_R, t_path_ediff, 0.0, path.m_HopCount - path.m_OscHopCount);
            m_SimResult->m_HStateEnergyDifference.AddValue(HP::HOPS_PER_PT, t_path_ediff, static_cast<double>(path.m_HopCount));
            m_SimResult->m_HStateEnergyDifference.AddValue(HP::OSC_PER_PT, t_path_ediff, static_cast<double>(path.m_OscHopCount));
            m_SimResult->m_HStateEnergyDifference.AddValue(HP::DIST_PER_HOP, t_path_ediff, t_path_dist, path.m_HopCount);
            m_SimResult->m_HStateEnergyDifference.AddValue(HP::DIST_PER_NONOSC, t_path_ediff, t_path_dist, path.m_HopCount - path.m_OscHopCount);
            //
            if (m_SimResult->m_HFieldEnergyContribution.m_IsConfigured)
            {
                m_SimResult->m_HFieldEnergyContribution.IncCount(HP::USED_PT, t_field_contrib);
                m_SimResult->m_HFieldEnergyContribution.AddValue(HP::USED_R, t_field_contrib, 1.0);
                m_SimResult->m_HFieldEnergyContribution.IncCount(HP::HOPS, t_field_contrib, path.m_HopCount);
                m_SimResult->m_HFieldEnergyContribution.IncCount(HP::OSC, t_field_contrib, path.m_OscHopCount);
                m_SimResult->m_HFieldEnergyContribution.IncCount(HP::NONOSC, t_field_contrib, path.m_HopCount - path.m_OscHopCount);
                m_SimResult->m_HFieldEnergyContribution.AddValue(HP::OSC_R, t_field_contrib, 1.0, path.m_OscHopCount);
                m_SimResult->m_HFieldEnergyContribution.AddValue(HP::OSC_R, t_field_contrib, 0.0, path.m_HopCount - path.m_OscHopCount);
                m_SimResult->m_HFieldEnergyContribution.AddValue(HP::HOPS_PER_PT, t_field_contrib, static_cast<double>(path.m_HopCount));
                m_SimResult->m_HFieldEnergyContribution.AddValue(HP::OSC_PER_PT, t_field_contrib, static_cast<double>(path.m_OscHopCount));
                m_SimResult->m_HFieldEnergyContribution.AddValue(HP::DIST_PER_HOP, t_field_contrib, t_path_dist, path.m_HopCount);
                m_SimResult->m_HFieldEnergyContribution.AddValue(HP::DIST_PER_NONOSC, t_field_contrib, t_path_dist, path.m_HopCount - path.m_OscHopCount);
            }
            //
            m_SimResult->m_HDistance.IncCount(HP::USED_PT, t_path_dist);
            m_SimResult->m_HDistance.AddValue(HP::USED_R, t_path_dist, 1.0);
            m_SimResult->m_HDistance.IncCount(HP::HOPS, t_path_dist, path.m_HopCount);
            m_SimResult->m_HDistance.IncCount(HP::OSC, t_path_dist, path.m_OscHopCount);
            m_SimResult->m_HDistance.IncCount(HP::NONOSC, t_path_dist, path.m_HopCount - path.m_OscHopCount);
            m_SimResult->m_HDistance.AddValue(HP::OSC_R, t_path_dist, 1.0, path.m_OscHopCount);
            m_SimResult->m_HDistance.AddValue(HP::OSC_R, t_path_dist, 0.0, path.m_HopCount - path.m_OscHopCount);
            m_SimResult->m_HDistance.AddValue(HP::HOPS_PER_PT, t_path_dist, static_cast<double>(path.m_HopCount));
            m_SimResult->m_HDistance.AddValue(HP::OSC_PER_PT, t_path_dist, static_cast<double>(path.m_OscHopCount));
            m_SimResult->m_HDistance.AddValue(HP::EDIFF_PER_HOP, t_path_dist, t_path_ediff + t_field_contrib, path.m_HopCount);
            m_SimResult->m_HDistance.AddValue(HP::EDIFF_PER_NONOSC, t_path_dist, t_path_ediff + t_field_contrib, path.m_HopCount - path.m_OscHopCount);
        }
    }

    // States
    for (std::size_t i = 0; i < m_Structure.size(); ++i)
    {
        const TLocalState& state = m_Structure[i];
        double t_state_energy = state.m_Energy * enfac + enmax;

        // Any state
        outer_states.check(state.m_Energy);
        //
        m_SimResult->m_HStateEnergy.IncCount(HP::ST, t_state_energy);
        m_SimResult->m_HStateEnergy.AddValue(HP::OCC_RTIME, t_state_energy, state.m_OccTime/m_TotalTime);
        //
        m_SimResult->m_HRelOccTime.IncCount(HP::ST, state.m_OccTime/m_TotalTime);
        m_SimResult->m_HRelOccTime.AddValue(HP::ST_EGY, state.m_OccTime/m_TotalTime, t_state_energy);
        m_SimResult->m_HRelOccTime.AddValue(HP::OUT_HOPS, state.m_OccTime/m_TotalTime, state_hops[i]);
        m_SimResult->m_HRelOccTime.AddValue(HP::OUT_OSC, state.m_OccTime/m_TotalTime, state_oscs[i]);
        m_SimResult->m_HRelOccTime.AddValue(HP::OUT_NONOSC, state.m_OccTime/m_TotalTime, state_nonoscs[i]);
        
        // Used state
        if (state_is_used[i])
        {
            outer_used_states.check(state.m_Energy);
            //
            m_SimResult->m_HStateEnergy.IncCount(HP::USED_ST, t_state_energy);
        }

        // Occupied state
        if (state.m_ElectronID < m_ParamSet->m_StateCount)
        {
            m_SimResult->m_HStateEnergy.AddValue(HP::OCC_R, t_state_energy, 1.0);
            //
            m_SimResult->m_HRelOccTime.IncCount(HP::OCC_ST, state.m_OccTime/m_TotalTime);
            if ((m_Electrons[state.m_ElectronID].m_HopCount != 0) &&
                (m_Electrons[state.m_ElectronID].m_HopCount == m_Electrons[state.m_ElectronID].m_OscHopCount))
                m_SimResult->m_HRelOccTime.IncCount(HP::OSC_ST, state.m_OccTime/m_TotalTime);
        }
        // Unoccupied state
        else
        {
            m_SimResult->m_HStateEnergy.AddValue(HP::OCC_R, t_state_energy, 0.0);
        }

        // State has outgoing hops
        if (state_hops[i] != 0) 
        {
            m_SimResult->m_HHopCount.IncCount(HP::ST, state_hops[i]);
            if (state.m_ElectronID < m_ParamSet->m_StateCount)
            {
                m_SimResult->m_HHopCount.IncCount(HP::OCC_ST, state_hops[i]);
                if ((m_Electrons[state.m_ElectronID].m_HopCount != 0) &&
                    (m_Electrons[state.m_ElectronID].m_HopCount == m_Electrons[state.m_ElectronID].m_OscHopCount))
                    m_SimResult->m_HHopCount.IncCount(HP::OSC_ST, state_hops[i]);
            }
        }

        // State has outgoing osc. hops
        if (state_oscs[i] != 0) 
        {
            m_SimResult->m_HOscHopCount.IncCount(HP::ST, state_oscs[i]);
            if (state.m_ElectronID < m_ParamSet->m_StateCount)
            {
                m_SimResult->m_HOscHopCount.IncCount(HP::OCC_ST, state_oscs[i]);
                if ((m_Electrons[state.m_ElectronID].m_HopCount != 0) &&
                    (m_Electrons[state.m_ElectronID].m_HopCount == m_Electrons[state.m_ElectronID].m_OscHopCount))
                    m_SimResult->m_HOscHopCount.IncCount(HP::OSC_ST, state_oscs[i]);
            }
        }

        // State has outgoing non-osc. hops
        if (state_nonoscs[i] != 0) 
        {
            m_SimResult->m_HNonOscHopCount.IncCount(HP::ST, state_nonoscs[i]);
            if (state.m_ElectronID < m_ParamSet->m_StateCount)
            {
                m_SimResult->m_HNonOscHopCount.IncCount(HP::OCC_ST, state_nonoscs[i]);
                if ((m_Electrons[state.m_ElectronID].m_HopCount != 0) &&
                    (m_Electrons[state.m_ElectronID].m_HopCount == m_Electrons[state.m_ElectronID].m_OscHopCount))
                    m_SimResult->m_HNonOscHopCount.IncCount(HP::OSC_ST, state_nonoscs[i]);
            }
        }

        // Calculated effective carriers
        m_SimResult->m_HStateEnergy.AddValue(HP::EFF_EL_CALC, t_state_energy, 
            (1.0/(exp((t_state_energy - m_SimResult->m_EffChemPot)/(Constant::kboltz*m_SimResult->m_EffTemp)) + 1.0)) *
            (1.0 - 1.0/(exp((t_state_energy - m_SimResult->m_EffChemPot)/(Constant::kboltz*m_SimResult->m_EffTemp)) + 1.0)));
    }

    // Occupation fit curve
    {
        auto state_energy_histo_axis = m_SimResult->m_HStateEnergy.GetBinCenters();
        for (const double& st_en : state_energy_histo_axis)
        {
            m_SimResult->m_HStateEnergy.AddValue(HP::OCC_FIT, st_en, 
                1.0/(exp((st_en - m_SimResult->m_EffChemPot)/(Constant::kboltz*m_SimResult->m_EffTemp)) + 1.0));
        }
    }

    // Sort state indices according to:
    // 0: max hops
    // 1: max non-osc hops (if equal max hops)
    // 2: max product of occupied and unoccupied time
    std::vector<std::vector<std::size_t>> state_sort (3,std::vector<std::size_t>(m_Structure.size(),0));
    std::iota(state_sort[0].begin(),state_sort[0].end(),0);
    std::iota(state_sort[1].begin(),state_sort[1].end(),0);
    std::iota(state_sort[2].begin(),state_sort[2].end(),0);
    std::sort(state_sort[0].begin(),state_sort[0].end(),
        [&](std::size_t i, std::size_t j)
        {
            return state_hops[i] > state_hops[j];
        });
    std::sort(state_sort[1].begin(),state_sort[1].end(),
        [&](std::size_t i, std::size_t j)
        {
            if (state_nonoscs[i] == state_nonoscs[j]) return state_hops[i] > state_hops[j];
            return state_nonoscs[i] > state_nonoscs[j];
        });
    std::sort(state_sort[2].begin(),state_sort[2].end(),
        [&](std::size_t i, std::size_t j)
        {
            return (m_Structure[i].m_OccTime/m_TotalTime)*(1.0 - m_Structure[i].m_OccTime/m_TotalTime) 
                > (m_Structure[j].m_OccTime/m_TotalTime)*(1.0 - m_Structure[j].m_OccTime/m_TotalTime);
        });

    // States (sorted by max sum of hops)
    std::uint32_t assigned_neff = 0;
    for (std::size_t i = 0; (i < m_Structure.size()) && (assigned_neff < static_cast<std::uint32_t>(m_SimResult->m_EffCarriers + 0.5)); ++i)
    {
        if (m_Structure[state_sort[0][i]].m_ElectronID < m_ParamSet->m_StateCount)
        {
            m_SimResult->m_HStateEnergy.IncCount(HP::EFF_EL_A, m_Structure[state_sort[0][i]].m_Energy * enfac + enmax);
            ++assigned_neff;
        }
    }

    // States (sorted by max sum of non-osc. hops)
    assigned_neff = 0;
    for (std::size_t i = 0; (i < m_Structure.size()) && (assigned_neff < static_cast<std::uint32_t>(m_SimResult->m_EffCarriers + 0.5)); ++i)
    {
        if (m_Structure[state_sort[1][i]].m_ElectronID < m_ParamSet->m_StateCount)
        {
            m_SimResult->m_HStateEnergy.IncCount(HP::EFF_EL_B, m_Structure[state_sort[1][i]].m_Energy * enfac + enmax);
            ++assigned_neff;
        }
    }

    // States (sorted by max product of occupied and unoccupied time)
    assigned_neff = 0;
    for (std::size_t i = 0; (i < m_Structure.size()) && (assigned_neff < static_cast<std::uint32_t>(m_SimResult->m_EffCarriers + 0.5)); ++i)
    {
        if (m_Structure[state_sort[2][i]].m_ElectronID < m_ParamSet->m_StateCount)
        {
            m_SimResult->m_HStateEnergy.IncCount(HP::EFF_EL_C, m_Structure[state_sort[2][i]].m_Energy * enfac + enmax);
            ++assigned_neff;
        }
    }

    // Electrons
    for (std::size_t i = 0; i < m_Electrons.size(); ++i)
    {
        const TElectron& electron = m_Electrons[i];
        double t_electron_energy = m_Structure[electron.m_CurrentStateID].m_Energy * enfac + enmax;
        double t_electron_disp = sqrt(electron.m_Disp_x * electron.m_Disp_x + electron.m_Disp_y * 
            electron.m_Disp_y + electron.m_Disp_z * electron.m_Disp_z) * spfac;

        // Any electron
        outer_electrons.check(m_Structure[electron.m_CurrentStateID].m_Energy);
        //
        m_SimResult->m_HStateEnergy.IncCount(HP::EL, t_electron_energy);
        //
        if (!electron_is_blocked[i])
        {
            m_SimResult->m_HNextTime.IncCount(HP::EL, electron.m_MinTime);
        }

        // Stationary electron (0 hops)
        if (electron.m_HopCount == 0)
        {
            m_SimResult->m_HStateEnergy.IncCount(HP::ZEROHOP_EL, t_electron_energy);
            m_SimResult->m_HStateEnergy.AddValue(HP::MOB_R, t_electron_energy, 0.0);
            //
            if (!electron_is_blocked[i])
            {
                m_SimResult->m_HNextTime.IncCount(HP::ZEROHOP_EL, electron.m_MinTime);
            }

            continue;
        }

        // Mobile or oscillating electron (> 0 hops)
        m_SimResult->m_HDisp.IncCount(HP::INCOSC_EL, t_electron_disp);
        m_SimResult->m_HDisp_x.IncCount(HP::INCOSC_EL, electron.m_Disp_x * spfac);
        m_SimResult->m_HDisp_y.IncCount(HP::INCOSC_EL, electron.m_Disp_y * spfac);
        m_SimResult->m_HDisp_z.IncCount(HP::INCOSC_EL, electron.m_Disp_z * spfac);
        //
        m_SimResult->m_HHopCount.IncCount(HP::INCOSC_EL, static_cast<double>(electron.m_HopCount));
        m_SimResult->m_HHopCount.AddValue(HP::DISP_PER_INCOSC_EL, static_cast<double>(electron.m_HopCount), t_electron_disp);
        //
        if (m_SimResult->m_HOscHopCount.m_IsConfigured)
        {
            m_SimResult->m_HOscHopCount.IncCount(HP::INCOSC_EL, static_cast<double>(electron.m_OscHopCount));
            m_SimResult->m_HOscHopCount.AddValue(HP::DISP_PER_INCOSC_EL, static_cast<double>(electron.m_OscHopCount), t_electron_disp);
        }

        // Mobile or oscillating electron (> 0 hops) that became mobile after equilibration
        if ((m_ParamSet->m_EqHopLimit != 0) && (electron.m_FirstHopTime > m_TotalTime - m_SimResult->m_TotalTime))
        {
            m_SimResult->m_HDisp.IncCount(HP::NEW_INCOSC_EL, t_electron_disp);
            m_SimResult->m_HDisp_x.IncCount(HP::NEW_INCOSC_EL, electron.m_Disp_x * spfac);
            m_SimResult->m_HDisp_y.IncCount(HP::NEW_INCOSC_EL, electron.m_Disp_y * spfac);
            m_SimResult->m_HDisp_z.IncCount(HP::NEW_INCOSC_EL, electron.m_Disp_z * spfac);
            //
            m_SimResult->m_HHopCount.IncCount(HP::NEW_INCOSC_EL, static_cast<double>(electron.m_HopCount));
            //
            if (m_SimResult->m_HOscHopCount.m_IsConfigured)
            {
                m_SimResult->m_HOscHopCount.IncCount(HP::NEW_INCOSC_EL, static_cast<double>(electron.m_OscHopCount));
            }
        }

        // Oscillating electron (only osc. hops)
        if (electron.m_HopCount == electron.m_OscHopCount)
        {
            m_SimResult->m_HStateEnergy.IncCount(HP::OSC_EL, t_electron_energy);
            //
            m_SimResult->m_HHopCount.IncCount(HP::OSC_EL, static_cast<double>(electron.m_HopCount));
            //
            if (m_SimResult->m_HOscHopCount.m_IsConfigured)
            {
                m_SimResult->m_HOscHopCount.IncCount(HP::OSC_EL, static_cast<double>(electron.m_OscHopCount));
            }
            //
            if (!electron_is_blocked[i])
            {
                m_SimResult->m_HNextTime.IncCount(HP::OSC_EL, electron.m_MinTime);
            }

            // Oscillating electron (only osc. hops) that became mobile after equilibration
            if ((m_ParamSet->m_EqHopLimit != 0) && (electron.m_FirstHopTime > m_TotalTime - m_SimResult->m_TotalTime))
            {
                m_SimResult->m_HStateEnergy.IncCount(HP::NEWOSC_EL, t_electron_energy);
                //
                m_SimResult->m_HHopCount.IncCount(HP::NEWOSC_EL, static_cast<double>(electron.m_HopCount));
                //
                if (m_SimResult->m_HOscHopCount.m_IsConfigured)
                {
                    m_SimResult->m_HOscHopCount.IncCount(HP::NEWOSC_EL, static_cast<double>(electron.m_OscHopCount));
                }
            }

            continue;
        }

        // Zero-displacement electron (but not zero-hop and not osc.-only)
        if (t_electron_disp < 0.5 * dist_per_used_path.min)
        {
            m_SimResult->m_HStateEnergy.IncCount(HP::ZERODISP_EL, t_electron_energy);
            //
            m_SimResult->m_HHopCount.IncCount(HP::ZERODISP_EL, static_cast<double>(electron.m_HopCount));
            //
            if (m_SimResult->m_HOscHopCount.m_IsConfigured)
            {
                m_SimResult->m_HOscHopCount.IncCount(HP::ZERODISP_EL, static_cast<double>(electron.m_OscHopCount));
            }
            //
            if (m_SimResult->m_HNonOscHopCount.m_IsConfigured)
            {
                m_SimResult->m_HNonOscHopCount.IncCount(HP::ZERODISP_EL, 
                    static_cast<double>(electron.m_HopCount - electron.m_OscHopCount));
            }
        }
        
        // Mobile electron (> 0 non-osc. hops)
        outer_mobile_electrons.check(m_Structure[electron.m_CurrentStateID].m_Energy);
        //
        m_SimResult->m_HStateEnergy.IncCount(HP::MOB_EL, t_electron_energy);
        m_SimResult->m_HStateEnergy.AddValue(HP::MOB_R, t_electron_energy, 1.0);
        m_SimResult->m_HStateEnergy.AddValue(HP::HOPS_PER_EL, t_electron_energy, static_cast<double>(electron.m_HopCount));
        m_SimResult->m_HStateEnergy.AddValue(HP::OSC_PER_EL, t_electron_energy, static_cast<double>(electron.m_OscHopCount));
        m_SimResult->m_HStateEnergy.AddValue(HP::NONOSC_PER_EL, t_electron_energy, static_cast<double>(electron.m_HopCount - electron.m_OscHopCount));
        //
        m_SimResult->m_HDisp.IncCount(HP::EL, t_electron_disp);
        m_SimResult->m_HDisp_x.IncCount(HP::EL, electron.m_Disp_x * spfac);
        m_SimResult->m_HDisp_y.IncCount(HP::EL, electron.m_Disp_y * spfac);
        m_SimResult->m_HDisp_z.IncCount(HP::EL, electron.m_Disp_z * spfac);
        m_SimResult->m_HDisp.AddValue(HP::HOPS_PER_EL, t_electron_disp, static_cast<double>(electron.m_HopCount));
        m_SimResult->m_HDisp_x.AddValue(HP::HOPS_PER_EL, electron.m_Disp_x * spfac, static_cast<double>(electron.m_HopCount));
        m_SimResult->m_HDisp_y.AddValue(HP::HOPS_PER_EL, electron.m_Disp_y * spfac, static_cast<double>(electron.m_HopCount));
        m_SimResult->m_HDisp_z.AddValue(HP::HOPS_PER_EL, electron.m_Disp_z * spfac, static_cast<double>(electron.m_HopCount));
        m_SimResult->m_HDisp.AddValue(HP::OSC_PER_EL, t_electron_disp, static_cast<double>(electron.m_OscHopCount));
        m_SimResult->m_HDisp_x.AddValue(HP::OSC_PER_EL, electron.m_Disp_x * spfac, static_cast<double>(electron.m_OscHopCount));
        m_SimResult->m_HDisp_y.AddValue(HP::OSC_PER_EL, electron.m_Disp_y * spfac, static_cast<double>(electron.m_OscHopCount));
        m_SimResult->m_HDisp_z.AddValue(HP::OSC_PER_EL, electron.m_Disp_z * spfac, static_cast<double>(electron.m_OscHopCount));
        m_SimResult->m_HDisp.AddValue(HP::NONOSC_PER_EL, t_electron_disp, 
            static_cast<double>(electron.m_HopCount - electron.m_OscHopCount));
        m_SimResult->m_HDisp_x.AddValue(HP::NONOSC_PER_EL, electron.m_Disp_x * spfac, 
            static_cast<double>(electron.m_HopCount - electron.m_OscHopCount));
        m_SimResult->m_HDisp_y.AddValue(HP::NONOSC_PER_EL, electron.m_Disp_y * spfac, 
            static_cast<double>(electron.m_HopCount - electron.m_OscHopCount));
        m_SimResult->m_HDisp_z.AddValue(HP::NONOSC_PER_EL, electron.m_Disp_z * spfac, 
            static_cast<double>(electron.m_HopCount - electron.m_OscHopCount));
        //
        m_SimResult->m_HHopCount.IncCount(HP::EL, static_cast<double>(electron.m_HopCount));
        m_SimResult->m_HHopCount.AddValue(HP::OSC_R, static_cast<double>(electron.m_HopCount), 1.0, electron.m_OscHopCount);
        m_SimResult->m_HHopCount.AddValue(HP::OSC_R, static_cast<double>(electron.m_HopCount), 0.0, electron.m_HopCount - electron.m_OscHopCount);
        m_SimResult->m_HHopCount.AddValue(HP::DISP_PER_EL, static_cast<double>(electron.m_HopCount), t_electron_disp);
        //
        if (m_SimResult->m_HOscHopCount.m_IsConfigured)
        {
            m_SimResult->m_HOscHopCount.IncCount(HP::EL, static_cast<double>(electron.m_OscHopCount));
            m_SimResult->m_HOscHopCount.AddValue(HP::OSC_R, static_cast<double>(electron.m_OscHopCount), 1.0, electron.m_OscHopCount);
            m_SimResult->m_HOscHopCount.AddValue(HP::OSC_R, static_cast<double>(electron.m_OscHopCount), 0.0, electron.m_HopCount - electron.m_OscHopCount);
            m_SimResult->m_HOscHopCount.AddValue(HP::DISP_PER_EL, static_cast<double>(electron.m_OscHopCount), t_electron_disp);
        }
        //
        if (m_SimResult->m_HNonOscHopCount.m_IsConfigured)
        {
            m_SimResult->m_HNonOscHopCount.IncCount(HP::EL, 
                static_cast<double>(electron.m_HopCount - electron.m_OscHopCount));
            m_SimResult->m_HNonOscHopCount.AddValue(HP::NONOSC_R, 
                static_cast<double>(electron.m_HopCount - electron.m_OscHopCount),
                0.0, electron.m_OscHopCount);
            m_SimResult->m_HNonOscHopCount.AddValue(HP::NONOSC_R, 
                static_cast<double>(electron.m_HopCount - electron.m_OscHopCount),
                1.0, electron.m_HopCount - electron.m_OscHopCount);
            m_SimResult->m_HNonOscHopCount.AddValue(HP::DISP_PER_EL, 
                static_cast<double>(electron.m_HopCount - electron.m_OscHopCount), t_electron_disp);
        }
        //
        if (!electron_is_blocked[i])
        {
            m_SimResult->m_HNextTime.IncCount(HP::MOB_EL, electron.m_MinTime);
        }

        // Mobile electron (> 0 non-osc. hops) that became mobile after equilibration
        if ((m_ParamSet->m_EqHopLimit != 0) && (electron.m_FirstHopTime > m_TotalTime - m_SimResult->m_TotalTime))
        {
            m_SimResult->m_HStateEnergy.IncCount(HP::NEWMOB_EL, t_electron_energy);
            //
            m_SimResult->m_HDisp.IncCount(HP::NEW_EL, t_electron_disp);
            m_SimResult->m_HDisp_x.IncCount(HP::NEW_EL, electron.m_Disp_x * spfac);
            m_SimResult->m_HDisp_y.IncCount(HP::NEW_EL, electron.m_Disp_y * spfac);
            m_SimResult->m_HDisp_z.IncCount(HP::NEW_EL, electron.m_Disp_z * spfac);
            //
            m_SimResult->m_HHopCount.IncCount(HP::NEW_EL, static_cast<double>(electron.m_HopCount));
            //
            if (m_SimResult->m_HOscHopCount.m_IsConfigured)
            {
                m_SimResult->m_HOscHopCount.IncCount(HP::NEW_EL, static_cast<double>(electron.m_OscHopCount));
            }
            //
            if (m_SimResult->m_HNonOscHopCount.m_IsConfigured)
            {
                m_SimResult->m_HNonOscHopCount.IncCount(HP::NEW_EL, 
                    static_cast<double>(electron.m_HopCount - electron.m_OscHopCount));
            }
        }
    }

    // Evaluate which minimum path count would have produced a sufficient distance cutoff
    std::uint32_t sufficient_path_count = 0;
    for (const auto& state : m_Structure)
    {
        std::uint32_t shorter_paths = 0;
        for (const auto& path : state.m_Paths)
        {
            const double& old_x {state.m_Pos_x};
            double new_x {m_Structure[path.m_StateID].m_Pos_x};
            if (old_x - new_x > 0.5) new_x += 1.0;
            if (new_x - old_x > 0.5) new_x -= 1.0;
            const double& old_y {state.m_Pos_y};
            double new_y {m_Structure[path.m_StateID].m_Pos_y};
            if (old_y - new_y > 0.5) new_y += 1.0;
            if (new_y - old_y > 0.5) new_y -= 1.0;
            const double& old_z {state.m_Pos_z};
            double new_z {m_Structure[path.m_StateID].m_Pos_z};
            if (old_z - new_z > 0.5) new_z += 1.0;
            if (new_z - old_z > 0.5) new_z -= 1.0;
            double path_dist = sqrt((new_x - old_x)*(new_x - old_x) + (new_y - old_y)*
                (new_y - old_y) + (new_z - old_z)*(new_z - old_z)) * spfac;

            if (path_dist <= dist_per_used_path.max) shorter_paths++;
        }
        if ((sufficient_path_count == 0) || (shorter_paths < sufficient_path_count))
        {
            sufficient_path_count = shorter_paths + 1;
        }
    }

    if (m_VL >= Verbosity::MEDIUM) 
    {
        std::cout << "done" << std::endl;
        
        std::cout << "  Used path distance range: " << dist_per_used_path
            << " nm (sufficient minimum path count: " << sufficient_path_count << ")" << std::endl;
        std::cout << "  Used path energy difference range: " << ediff_per_used_path 
            << " eV (without electric field contribution)" << std::endl;
        std::cout << "  Used path time range: " << time_per_hop 
            << " s (without random factor; available path time range: " << time_per_path << " s)" << std::endl;
        std::cout << "  Used state energy range: [" << m_SimResult->m_MinUsedStateEnergy << ", " 
            << m_SimResult->m_MaxUsedStateEnergy << "] eV" << std::endl;

        std::cout << "  Hops with distance > 90 %" << " of cut-off: " << outer_distance_hops.upper << " (" 
            << 100.0 * static_cast<double>(outer_distance_hops.upper) / static_cast<double>(total_hops) << " %" << " of hops)" << std::endl;
        std::cout << "  Hops with state energy difference in outer 5 %" << "s of cut-off range: " << outer_ediff_hops.lower << " (" 
            << 100.0 * static_cast<double>(outer_ediff_hops.lower) / static_cast<double>(total_hops) << " %)"
            << " <--> " << outer_ediff_hops.upper << " (" 
            << 100.0 * static_cast<double>(outer_ediff_hops.upper) / static_cast<double>(total_hops) << " %) [% = of hops]" << std::endl;
        std::cout << "  Hops starting in outer 5 %" << "s of state energy range: " << outer_energy_hops.lower << " (" 
            << 100.0 * static_cast<double>(outer_energy_hops.lower) / static_cast<double>(total_hops) << " %)"
            << " <--> " << outer_energy_hops.upper << " (" 
            << 100.0 * static_cast<double>(outer_energy_hops.upper) / static_cast<double>(total_hops) << " %) [% = of hops]" << std::endl;
        
        std::cout << "  Mobile electrons in outer 5 %" << "s of energy range: " << outer_mobile_electrons.lower;
        if (outer_electrons.lower != 0) std::cout << " (" << 100.0 * static_cast<double>(outer_mobile_electrons.lower) / static_cast<double>(outer_electrons.lower) << " %)";
        std::cout << " <--> " << outer_mobile_electrons.upper;
        if (outer_electrons.upper != 0) std::cout << " (" << 100.0 * static_cast<double>(outer_mobile_electrons.upper) / static_cast<double>(outer_electrons.upper) << " %)";
        if ((outer_electrons.lower != 0) || (outer_electrons.upper != 0)) std::cout << " [% = of outer electrons]";
        std::cout << std::endl;
        std::cout << "  Used states in outer 5 %" << "s of energy range: " << outer_used_states.lower;
        if (outer_states.lower != 0) std::cout << " (" << 100.0 * static_cast<double>(outer_used_states.lower) / static_cast<double>(outer_states.lower) << " %)";
        std::cout << " <--> " << outer_used_states.upper;
        if (outer_states.upper != 0) std::cout << " (" << 100.0 * static_cast<double>(outer_used_states.upper) / static_cast<double>(outer_states.upper) << " %)";
        if ((outer_states.lower != 0) || (outer_states.upper != 0)) std::cout << " [% = of outer states]";
        std::cout << std::endl;
    }

    // Additional analysis of different electron types
    if (m_VL >= Verbosity::MAXIMUM) 
    {
        std::cout << "Analysis of electrons with zero displacement (but > 0 non-osc. hops):" << std::endl;
        std::uint32_t zd_electrons = 0;
        GF::TMinMaxMean<std::uint64_t> zd_hops_per_electron;
        GF::TMinMaxMean<std::uint64_t> zd_oscs_per_electron;
        GF::TMinMaxMean<std::uint64_t> zd_nonoscs_per_electron;
        for (const auto& electron : m_Electrons)
        {
            if ((electron.m_HopCount == 0) || (electron.m_HopCount == electron.m_OscHopCount)) continue;
            
            double disp = sqrt(electron.m_Disp_x * electron.m_Disp_x + electron.m_Disp_y * 
                electron.m_Disp_y + electron.m_Disp_z * electron.m_Disp_z) * spfac;

            if (disp < 0.5 * dist_per_used_path.min)
            {
                zd_electrons++;
                zd_hops_per_electron.check(electron.m_HopCount);
                zd_oscs_per_electron.check(electron.m_OscHopCount);
                zd_nonoscs_per_electron.check(electron.m_HopCount - electron.m_OscHopCount);
            }
        }
        std::cout << "  Electrons: " << zd_electrons << std::endl;
        if (zd_electrons != 0)
        {
            std::cout << "  Average hops per electron: " << zd_hops_per_electron << std::endl;
            std::cout << "  Average osc. hops per electron: " << zd_oscs_per_electron << std::endl;
            std::cout << "  Average non-osc. hops per electron: " << zd_nonoscs_per_electron << std::endl;
        }

        if (m_ParamSet->m_EqHopLimit != 0)
        {
            std::cout << "Analysis of electrons with first hop after equilibration:" << std::endl;
            std::uint32_t nw_electrons = 0, nw_osc = 0, nw_single = 0, nw_single_osc = 0;
            GF::TMinMaxMean<std::uint64_t> nw_hops_per_electron;
            GF::TMinMaxMean<std::uint64_t> nw_oscs_per_electron;
            GF::TMinMaxMean<std::uint64_t> nw_nonoscs_per_electron;
            GF::TMinMaxMean<double> nw_disp_per_electron;
            for (const auto& electron : m_Electrons)
            {
                if (electron.m_HopCount == 0) continue;

                if (electron.m_FirstHopTime > m_TotalTime - m_SimResult->m_TotalTime)
                {
                    nw_electrons++;
                    if (electron.m_HopCount == electron.m_OscHopCount) nw_osc++;
                    if (electron.m_HopCount == 1) nw_single++;
                    if ((electron.m_HopCount - electron.m_OscHopCount == 1) && (electron.m_OscHopCount != 0)) nw_single_osc++;
                    nw_hops_per_electron.check(electron.m_HopCount);
                    nw_oscs_per_electron.check(electron.m_OscHopCount);
                    nw_nonoscs_per_electron.check(electron.m_HopCount - electron.m_OscHopCount);
                    nw_disp_per_electron.check(sqrt(electron.m_Disp_x * electron.m_Disp_x + electron.m_Disp_y * 
                        electron.m_Disp_y + electron.m_Disp_z * electron.m_Disp_z) * spfac);
                }
            }
            std::cout << "  Electrons: " << nw_electrons << std::endl;
            if (nw_electrons != 0)
            {
                std::cout << "  - mobile: " << nw_electrons - nw_osc << std::endl;
                std::cout << "    - one hop: " << nw_single << std::endl;
                std::cout << "    - one hop except osc. hops: " << nw_single_osc << std::endl;
                std::cout << "  - oscillating: " << nw_osc << std::endl;
                std::cout << "  Average hops per electron: " << nw_hops_per_electron << std::endl;
                std::cout << "  Average osc. hops per electron: " << nw_oscs_per_electron << std::endl;
                std::cout << "  Average non-osc. hops per electron: " << nw_nonoscs_per_electron << std::endl;
                std::cout << "  Average displacement per electron: " << nw_disp_per_electron << " nm" << std::endl;
            }
        }
    }

    m_ResultsReady = true;
}