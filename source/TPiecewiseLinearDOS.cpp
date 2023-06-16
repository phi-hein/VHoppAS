#include "TPiecewiseLinearDOS.hpp"

#include <utility>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <functional>
#include <cmath>
#include <regex>

#include "Constants.hpp"
#include "GlobalFunctions.hpp"
#include "CustomExceptions.hpp"

// Type of this DOS sub-class
const std::string MC::TPiecewiseLinearDOS::m_Type = "PiecewiseLinear";

// Parameter descriptors
const std::array<std::string,2> MC::TPiecewiseLinearDOS::s_Type = {"Type",""};
const std::array<std::string,2> MC::TPiecewiseLinearDOS::s_RefTemp = {"RefTemp","K"};

// Default constructor
MC::TPiecewiseLinearDOS::TPiecewiseLinearDOS()
    : m_HasDOS(false), m_Ready(false), m_RefTemp(0.0), m_MaxEnergy(0.0), m_MinEnergy(0.0), 
    m_SpatialFactor(0.0), m_FermiLevel(0.0)
{
    
}

// Get copy of the internal object
std::unique_ptr<MC::TDOS> MC::TPiecewiseLinearDOS::Copy()
{
    try
    {
        return std::unique_ptr<TDOS>(new TPiecewiseLinearDOS(*this));
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot create copy of DOS.",__func__,e.what());
    }
}

// Set DOS values (from file content)
void MC::TPiecewiseLinearDOS::SpecifyDOS(const std::string& dos_content)
{
    double ref_temp;
    std::vector<double> energy_values;
    std::vector<double> dos_values;
    std::smatch match;

    // Read reference temperature
    if (std::regex_search(dos_content, match,
        std::regex(GF::DescRegex(s_RefTemp) + "\\s*=\\s*(" + Constant::dblex + ")")))
    {
        if (!GF::StringToDouble(match[1],ref_temp))
        {
            throw EX::TInvalidInput("Cannot read reference temperature from DOS file.");
        }
    }
    else throw EX::TInvalidInput("DOS file does not contain reference temperature.");

    // Read DOS values
    try
    {
        std::string::const_iterator start(dos_content.cbegin());
        while (std::regex_search(start,dos_content.cend(), match,
            std::regex("(" + Constant::dblex + ")\\s+(" + Constant::dblex + ")")))
        {
            double energy;
            double density;
            if ((!GF::StringToDouble(match[1],energy)) || (!GF::StringToDouble(match[2],density)))
            {
                throw EX::TInvalidInput("Cannot read DOS values.");
            }
            energy_values.push_back(energy);
            dos_values.push_back(density);
            start = match.suffix().first;
        }
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot read all DOS values.",__func__,e.what());
    }

    // Validate DOS data
    if (energy_values.size() == 0)
    {
        throw EX::TInvalidInput("No DOS data found.");
    }
    if (energy_values.size() == 1)
    {
        throw EX::TInvalidInput("DOS has only one point.");
    }
    if (std::adjacent_find(energy_values.cbegin(),energy_values.cend(), 
        std::greater_equal<double>()) != energy_values.cend())
    {
        throw EX::TInvalidInput("Non-ascending DOS energies.");
    }
    if (std::any_of(dos_values.cbegin(), dos_values.cend(), [](double d){ return d < 0.0; }))
    {
        throw EX::TInvalidInput("Negative DOS values.");
    }
    if (std::all_of(dos_values.cbegin(), dos_values.cend(), [](double d){ return d == 0.0; }))
    {
        throw EX::TInvalidInput("DOS is zero everywhere.");
    }
    if (ref_temp < 0.0)
    {
        throw EX::TInvalidInput("Negative reference temperature of DOS.");
    }
    if (ref_temp == 0.0)
    {
        if ((energy_values.front() >= 0.0) || (energy_values.back() < 0.0))
        {
            throw EX::TInvalidInput("Fermi level is outside of DOS energy axis (Ref-T = 0 K).");
        }
        auto it = std::lower_bound(energy_values.cbegin(),energy_values.cend(),0.0);
        ++it;
        if (std::all_of(dos_values.cbegin(), dos_values.cbegin() + std::distance(energy_values.cbegin(),it), 
            [](double d){ return d == 0.0; }))
        {
            throw EX::TInvalidInput("DOS is zero everywhere below Fermi level (Ref-T = 0 K).");
        }
    }

    m_EnergyValues = std::move(energy_values);
    m_DOSValues = std::move(dos_values);
    m_RefTemp = ref_temp;
    m_HasDOS = true;
    DeleteParameters();

    std::cout << "  Type: piecewise-linear" << std::endl;
    std::cout << "  Energy range: [" << m_EnergyValues.front() << ", " << m_EnergyValues.back() << "] eV" << std::endl;
    std::cout << "  Reference temperature: " << m_RefTemp << " K" << std::endl;
}

// Validate additional parameters for DOS
void MC::TPiecewiseLinearDOS::ValidateParameters(const TParamSet& pset) const
{
    if (m_HasDOS == false)
    {
        throw EX::TInvalidStatus("No DOS data present.",__func__);
    }
    if (pset.m_MinStateEnergy >= pset.m_MaxStateEnergy)
    {
        throw EX::TInvalidInput("Minimum energy equal or higher than maximum energy.");
    }
    if (pset.m_MinStateEnergy < m_EnergyValues.front())
    {
        throw EX::TInvalidInput("Minimum state energy exceeds DOS.");
    }
    if (pset.m_MaxStateEnergy > m_EnergyValues.back())
    {
        throw EX::TInvalidInput("Maximum state energy exceeds DOS.");
    }
    auto range_start = m_DOSValues.cbegin() + std::distance(m_EnergyValues.cbegin(),
        std::lower_bound(m_EnergyValues.cbegin(),m_EnergyValues.cend(),pset.m_MinStateEnergy));
    if (range_start != m_DOSValues.cbegin()) --range_start;
    auto range_end = m_DOSValues.cbegin() + std::distance(m_EnergyValues.cbegin(),
        std::upper_bound(m_EnergyValues.cbegin(),m_EnergyValues.cend(),pset.m_MaxStateEnergy));
    if (range_end != m_DOSValues.cend()) ++range_end;
    if (std::all_of(range_start, range_end, [](double d){ return d == 0.0; }))
    {
        throw EX::TInvalidInput("DOS is zero in selected range.");
    }
}

// Set additional parameters for DOS
void MC::TPiecewiseLinearDOS::ApplyParameters(const TParamSet& pset, std::string& output)
{
    // Adjust energy axis to temperature
    // (set reference temperature in any case for later consistency with GetAdjustedFermiLevel)
    std::stringstream sstr;
    if (pset.m_EFTAdjust)
    {
        m_Ready = false;
        double new_EF = GetAdjustedFermiLevel(pset.m_Temperature);
        if (new_EF != 0.0)
        {
            for (double& energy: m_EnergyValues) energy -= new_EF;
            
            sstr << "Adjustment of charge neutral Fermi level (= zero on DOS energy axis) from " 
                << m_RefTemp << " K to " << pset.m_Temperature << " K: " << new_EF << " eV";
        }
    }
    else
    {
        sstr << "! Warning: temperature-adjustment of charge neutral Fermi level (= zero on DOS energy axis) is disabled -> results might belong to a different chemical potential !";
    }
    m_RefTemp = pset.m_Temperature;
    output = sstr.str();

    // Validate DOS-related parameters
    ValidateParameters(pset);

    // Set energy extrema
    m_MinEnergy = pset.m_MinStateEnergy;
    m_MaxEnergy = pset.m_MaxStateEnergy;
    m_FermiLevel = (pset.m_ChemPot - m_MaxEnergy)/(m_MaxEnergy - m_MinEnergy);

    try
    {    
        // Generate cropped and scaled energy axis (-1 to 0) and DOS values
        m_RelEnergyValues = std::vector<double>();
        m_RelDOSValues = std::vector<double>();
        for (std::size_t i = 1; i < m_EnergyValues.size(); ++i)
        {
            if (m_EnergyValues[i] > m_MinEnergy)
            {
                if (m_RelEnergyValues.empty())
                {
                    m_RelEnergyValues.push_back(-1.0);
                    m_RelDOSValues.push_back(m_DOSValues[i-1] + (m_DOSValues[i] - m_DOSValues[i-1]) * 
                        (m_MinEnergy - m_EnergyValues[i-1]) / (m_EnergyValues[i] - m_EnergyValues[i-1]));
                }
                if (m_EnergyValues[i] < m_MaxEnergy)
                {
                    m_RelEnergyValues.push_back((m_EnergyValues[i] - m_MinEnergy)/(m_MaxEnergy-m_MinEnergy)-1.0);
                    m_RelDOSValues.push_back(m_DOSValues[i]);
                }
                else
                {
                    m_RelEnergyValues.push_back(0.0);
                    m_RelDOSValues.push_back(m_DOSValues[i-1] + (m_DOSValues[i] - m_DOSValues[i-1]) * 
                        (m_MaxEnergy - m_EnergyValues[i-1]) / (m_EnergyValues[i] - m_EnergyValues[i-1]));
                    break;
                }
            }
        }

        // Set piecewise linear probability distribution
        m_DOSDistribution = std::piecewise_linear_distribution<double>(
            m_RelEnergyValues.cbegin(),m_RelEnergyValues.cend(),m_RelDOSValues.cbegin());
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot create scaled DOS distribution.",__func__,e.what());
    }

    // Calculate spatial scaling factor: simulation cell is a unit cube in relative units (1.0 x 1.0 x 1.0)
    auto energy_it = m_RelEnergyValues.cbegin();
    auto dos_it = m_RelDOSValues.cbegin();
    auto energy_next = std::next(energy_it);
    auto dos_next = std::next(dos_it);
    double integral = 0.0;
    while (energy_next != m_RelEnergyValues.cend())
    {
        integral += 0.5*(*dos_it*1.0E-21 + *dos_next*1.0E-21)*(*energy_next - *energy_it);
        ++energy_it;
        ++energy_next;
        ++dos_it;
        ++dos_next;
    }
    m_SpatialFactor = std::cbrt(static_cast<double>(pset.mO_StateCount())/(integral*(m_MaxEnergy - m_MinEnergy)));

    m_Ready = true;
}

// Remove applied parameters (except the raw DOS data)
void MC::TPiecewiseLinearDOS::DeleteParameters()
{
    m_Ready = false; 
    m_MaxEnergy = 0.0;
    m_MinEnergy = 0.0;
    m_FermiLevel = 0.0;
    m_SpatialFactor = 0.0;
    m_RelEnergyValues = std::vector<double>();
    m_RelDOSValues = std::vector<double>();
    m_DOSDistribution = std::piecewise_linear_distribution<double>();
}

// Returns whether DOS is properly defined
bool MC::TPiecewiseLinearDOS::HasDOS() const
{
    return m_HasDOS;
}

// Returns whether DOS and its internal parameters (e.g. scaling) are properly defined
bool MC::TPiecewiseLinearDOS::IsReady() const
{
    return (m_HasDOS && m_Ready);
}

// Scaling factor for spatial differences (in nm): absolute = factor*relative
double MC::TPiecewiseLinearDOS::GetSpatialFactor() const
{
    return m_SpatialFactor;
}

// Scaling factor for energy differences (in eV): absolute = factor*relative
double MC::TPiecewiseLinearDOS::GetEnergyFactor() const
{
    return m_MaxEnergy - m_MinEnergy;
}

// Calculate F-D occupation probability based on state energy and kBT (both in relative energy)
double MC::TPiecewiseLinearDOS::GetOccupationProbability(double energy, double kBT) const
{
    return 1.0/(exp((energy - m_FermiLevel)/kBT) + 1.0);
}

// Generate random state energy (relative: -1 to 0) according to DOS
double MC::TPiecewiseLinearDOS::GetRandomEnergy(std::mt19937_64 &rng)
{
    return m_DOSDistribution(rng);
}

// Calculate expected number of electrons based on selected DOS range and Fermi level (kBT in relative energy)
std::uint32_t MC::TPiecewiseLinearDOS::GetElectronCount(double kBT) const
{
    if (IsReady() == false)
    {
        throw EX::TInvalidStatus("DOS object is not ready.",__func__);
    }

    auto energy_it = m_RelEnergyValues.cbegin();
    auto dos_it = m_RelDOSValues.cbegin();
    auto energy_next = std::next(energy_it);
    auto dos_next = std::next(dos_it);
    double integral = 0.0;

    // DOS times F-D distribution is integrated in at least deltaE kT intervals (in relative energy)
    double en_val = *energy_it;
    double dos_val = *dos_it*1.0E-21;
    while (energy_next != m_RelEnergyValues.cend())
    {
        bool loop_switch = true;
        while (loop_switch)
        {
            double next_en_val = en_val + Constant::deltaE * kBT;
            double next_dos_val = 0.0;

            if (next_en_val < *energy_next)
            {
                next_dos_val = *dos_it*1.0E-21 + (*dos_next*1.0E-21 - *dos_it*1.0E-21)*
                    (next_en_val - *energy_it)/(*energy_next - *energy_it);
            
            }
            else
            {
                next_en_val = *energy_next;
                next_dos_val = *dos_next*1.0E-21;
                loop_switch = false;
            }

            integral += 0.5*(next_en_val - en_val)*(
                dos_val/(1.0 + exp((en_val - m_FermiLevel)/kBT)) +
                next_dos_val/(1.0 + exp((next_en_val - m_FermiLevel)/kBT)));

            en_val = next_en_val;
            dos_val = next_dos_val;
        }
        ++energy_it;
        ++energy_next;
        ++dos_it;
        ++dos_next;
    }
    
    return static_cast<std::uint32_t>(integral*(m_MaxEnergy-m_MinEnergy)*pow(m_SpatialFactor,3.0) + 0.5);
}

// Calculate expected effective carrier density (in 1/cm3) based on chemical potential (in eV) and temperature (in K)
double MC::TPiecewiseLinearDOS::GetEffCarrierDensity(double chem_pot, double temp) const
{
    if (m_HasDOS == false)
    {
        throw EX::TInvalidStatus("DOS data not present.",__func__);
    }

    auto energy_it = m_EnergyValues.cbegin();
    auto dos_it = m_DOSValues.cbegin();
    auto energy_next = std::next(energy_it);
    auto dos_next = std::next(dos_it);
    double integral = 0.0;
    const double kBT = Constant::kboltz * temp;

    // DOS times f*(1-f) is integrated in at least deltaE kT intervals
    double en_val = *energy_it;
    double dos_val = *dos_it*1.0E-21;
    while (energy_next != m_EnergyValues.cend())
    {
        bool loop_switch = true;
        while (loop_switch)
        {
            double next_en_val = en_val + Constant::deltaE * kBT;
            double next_dos_val = 0.0;

            if (next_en_val < *energy_next)
            {
                next_dos_val = *dos_it*1.0E-21 + (*dos_next*1.0E-21 - *dos_it*1.0E-21)*
                    (next_en_val - *energy_it)/(*energy_next - *energy_it);
            
            }
            else
            {
                next_en_val = *energy_next;
                next_dos_val = *dos_next*1.0E-21;
                loop_switch = false;
            }

            integral += 0.5 * (next_en_val - en_val) * (
                dos_val / (1.0 + exp((en_val - chem_pot)/kBT)) 
                * (1.0 - 1.0 / (1.0 + exp((en_val - chem_pot)/kBT))) +
                next_dos_val / (1.0 + exp((next_en_val - chem_pot)/kBT))
                * (1.0 - 1.0 / (1.0 + exp((next_en_val - chem_pot)/kBT))));

            en_val = next_en_val;
            dos_val = next_dos_val;
        }
        ++energy_it;
        ++energy_next;
        ++dos_it;
        ++dos_next;
    }
    
    return integral*1.0E21;
}

// Calculate expected partial entropy (in eV/K) based on chemical potential (in eV) and temperature (in K)
double MC::TPiecewiseLinearDOS::GetPartialEntropy(double chem_pot, double temp) const
{
    if (m_HasDOS == false)
    {
        throw EX::TInvalidStatus("DOS data not present.",__func__);
    }

    auto energy_it = m_EnergyValues.cbegin();
    auto dos_it = m_DOSValues.cbegin();
    auto energy_next = std::next(energy_it);
    auto dos_next = std::next(dos_it);
    double neff_integral = 0.0;
    double ntherm_integral = 0.0;
    const double kBT = Constant::kboltz * temp;

    // DOS times f*(1-f) and times f*(1-f)*(E-EF)/kT is integrated in at least deltaE kT intervals
    double en_val = *energy_it;
    double dos_val = *dos_it*1.0E-21;
    while (energy_next != m_EnergyValues.cend())
    {
        bool loop_switch = true;
        while (loop_switch)
        {
            double next_en_val = en_val + Constant::deltaE * kBT;
            double next_dos_val = 0.0;

            if (next_en_val < *energy_next)
            {
                next_dos_val = *dos_it*1.0E-21 + (*dos_next*1.0E-21 - *dos_it*1.0E-21)*
                    (next_en_val - *energy_it)/(*energy_next - *energy_it);
            
            }
            else
            {
                next_en_val = *energy_next;
                next_dos_val = *dos_next*1.0E-21;
                loop_switch = false;
            }

            neff_integral += 0.5 * (next_en_val - en_val) * (
                dos_val / (1.0 + exp((en_val - chem_pot)/kBT)) 
                * (1.0 - 1.0 / (1.0 + exp((en_val - chem_pot)/kBT))) +
                next_dos_val / (1.0 + exp((next_en_val - chem_pot)/kBT))
                * (1.0 - 1.0 / (1.0 + exp((next_en_val - chem_pot)/kBT))));

            ntherm_integral += 0.5 * (next_en_val - en_val) * (
                dos_val / (1.0 + exp((en_val - chem_pot)/kBT)) 
                * (1.0 - 1.0 / (1.0 + exp((en_val - chem_pot)/kBT)))
                * (en_val - chem_pot)/kBT +
                next_dos_val / (1.0 + exp((next_en_val - chem_pot)/kBT))
                * (1.0 - 1.0 / (1.0 + exp((next_en_val - chem_pot)/kBT)))
                * (next_en_val - chem_pot)/kBT);

            en_val = next_en_val;
            dos_val = next_dos_val;
        }
        ++energy_it;
        ++energy_next;
        ++dos_it;
        ++dos_next;
    }
    
    return Constant::kboltz * ntherm_integral / neff_integral;
}

// Calculate Fermi level (in eV) adjusted to another temperature (in K)
double MC::TPiecewiseLinearDOS::GetAdjustedFermiLevel(double newT) const
{
    if (m_HasDOS == false)
    {
        throw EX::TInvalidStatus("DOS data not present.",__func__);
    }
    if (newT <= 0.0)
    {
        throw EX::TInvalidStatus("Non-positive temperature.",__func__);
    }

    // Lambda function: integrate DOS (in 1/nm3) up to Fermi level (for electron density at 0 K)
    auto calc_dos_integral = [&] (double i_EF) -> double
    {
        auto energy_it = m_EnergyValues.cbegin();
        auto dos_it = m_DOSValues.cbegin();
        auto energy_next = std::next(energy_it);
        auto dos_next = std::next(dos_it);
        double integral = 0.0;

        if ((i_EF <= m_EnergyValues.front()) || (i_EF >= m_EnergyValues.back())) return 0.0;

        while (*energy_next < i_EF)
        {
            integral += 0.5*(*dos_it*1.0E-21 + *dos_next*1.0E-21)*(*energy_next - *energy_it);
            ++energy_it;
            ++energy_next;
            ++dos_it;
            ++dos_next;
        }
        
        integral += 0.5*(*dos_it*1.0E-21 + *dos_it*1.0E-21 + 
            (*dos_next*1.0E-21 - *dos_it*1.0E-21)*((i_EF - *energy_it)/(*energy_next - *energy_it)))*(i_EF - *energy_it);

        return integral;
    };

    // Lambda function: integrate DOS times F-D distribution (in 1/nm3), interpolated on at least deltaE kT
    auto calc_electron_density = [&] (double i_EF, double i_T) -> double
    {
        auto energy_it = m_EnergyValues.cbegin();
        auto dos_it = m_DOSValues.cbegin();
        auto energy_next = std::next(energy_it);
        auto dos_next = std::next(dos_it);
        double integral = 0.0;

        double en_val = *energy_it;
        double dos_val = *dos_it*1.0E-21;
        while (energy_next != m_EnergyValues.cend())
        {
            bool loop_switch = true;
            while (loop_switch)
            {
                double next_en_val = en_val + Constant::deltaE * Constant::kboltz * i_T;
                double next_dos_val = 0.0;

                if (next_en_val < *energy_next)
                {
                    next_dos_val = *dos_it*1.0E-21 + (*dos_next*1.0E-21 - *dos_it*1.0E-21)*
                        (next_en_val - *energy_it)/(*energy_next - *energy_it);
                }
                else
                {
                    next_en_val = *energy_next;
                    next_dos_val = *dos_next*1.0E-21;
                    loop_switch = false;
                }

                integral += 0.5*(next_en_val - en_val)*(
                    dos_val/(1.0 + exp((en_val - i_EF)/(Constant::kboltz * i_T))) +
                    next_dos_val/(1.0 + exp((next_en_val - i_EF)/(Constant::kboltz * i_T))));

                en_val = next_en_val;
                dos_val = next_dos_val;
            }
            ++energy_it;
            ++energy_next;
            ++dos_it;
            ++dos_next;
        }
        
        return integral;
    };   
    
    double refEF = 0.0;
    if (m_Ready) refEF = m_FermiLevel*(m_MaxEnergy - m_MinEnergy) + m_MaxEnergy;
    if (newT == m_RefTemp) return refEF;

    // Calculate reference density (in 1/nm3)
    double refDensity = 0.0;
    if (m_RefTemp == 0.0)
        refDensity = calc_dos_integral(refEF);
    else
        refDensity = calc_electron_density(refEF, m_RefTemp);

    if (refDensity <= 0.0)
    {
        throw EX::TInvalidStatus("Invalid reference electron density.",__func__);
    }

    // Find appropriate boundaries for bisection search
    double lowerEF = refEF;
    double lowerDensity = calc_electron_density(refEF, newT);
    double upperEF = refEF;
    double upperDensity = 0.0;
    if (lowerDensity == refDensity) return refEF;
    if (lowerDensity < refDensity)
    {
        do
        {
            upperEF += 0.1;
            upperDensity = calc_electron_density(upperEF, newT);
        } while (upperDensity < refDensity);
        if (upperDensity == refDensity) return upperEF;
    }
    else
    {
        upperDensity = lowerDensity;
        do
        {
            lowerEF -= 0.1;
            lowerDensity = calc_electron_density(lowerEF, newT);
        } while (lowerDensity > refDensity);
        if (lowerDensity == refDensity) return lowerEF;
    }
    
    // Bisection search
    while (upperEF - lowerEF > 0.00001)
    {
        const double centerEF = 0.5*(upperEF + lowerEF);
        const double centerDensity = calc_electron_density(centerEF, newT);
        if (centerDensity > refDensity)
        {
            upperEF = centerEF;
            upperDensity = centerDensity;
        }
        else
        {
            lowerEF = centerEF;
            lowerDensity = centerDensity;
        }
    }
    
    return 0.5*(upperEF + lowerEF);
}