#ifndef TPiecewiseLinearDOS_H_
#define TPiecewiseLinearDOS_H_

#include <random>
#include <vector>
#include <memory>
#include <array>

#include "TDOS.hpp"

namespace MC
{

class TPiecewiseLinearDOS : public TDOS
{
public:
    // Type of this DOS sub-class
    static const std::string m_Type;

    // Parameter descriptors and units (without white-spaces)
    static const std::array<std::string,2> s_Type;
    static const std::array<std::string,2> s_RefTemp;
    static const std::string s_Data;

    // Default constructor
    TPiecewiseLinearDOS();

    // Get copy of the internal object
    std::unique_ptr<TDOS> Copy();

    // Set DOS values (from file content)
    void SpecifyDOS(const std::string& dos_content);

    // Validate additional parameters for DOS
    void ValidateParameters(const TParamSet& pset) const;

    // Set additional parameters for DOS
    void ApplyParameters(const TParamSet& pset, std::string& output);

    // Remove applied parameters (except the raw DOS data)
    void DeleteParameters();

    // Returns whether DOS is properly defined
    bool HasDOS() const;

    // Returns whether DOS and its internal parameters (e.g. scaling) are properly defined
    bool IsReady() const;

	// Scaling factor for spatial differences (in nm): absolute = factor*relative
    double GetSpatialFactor() const;

    // Scaling factor for energy differences (in eV): absolute = factor*relative
    double GetEnergyFactor() const;

    // Calculate F-D occupation probability based on state energy and kBT (both in relative energy)
    double GetOccupationProbability(double energy, double kBT) const;
    
    // Generate random state energy (relative: -1 to 0) according to DOS
	double GetRandomEnergy(std::mt19937_64 &rng);

    // Calculate expected number of electrons based on selected DOS range and Fermi level (kBT in relative energy)
    std::uint32_t GetElectronCount(double kBT) const;

    // Calculate expected effective carrier density (in 1/cm3) based on chemical potential (in eV) and temperature (in K)
    double GetEffCarrierDensity(double chem_pot, double temp) const;

    // Calculate expected partial entropy (in eV/K) based on chemical potential (in eV) and temperature (in K)
    double GetPartialEntropy(double chem_pot, double temp) const;

    // Calculate Fermi level (in eV) adjusted to another temperature (in K)
    double GetAdjustedFermiLevel(double newT) const;

private:
    // Ready switch: true = DOS values are properly specified
    bool m_HasDOS;

    // Ready switch: true = DOS and scaling parameters are properly specified
    bool m_Ready;

    // Reference temperature (in K; temperature of the Fermi level which is the zero reference of the absolute energy axis)
    double m_RefTemp;

    // Maximum energy (absolute in eV; corresponds to relative 0)
    double m_MaxEnergy;

    // Minimum energy (absolute in eV; corresponds to relative -1)
    double m_MinEnergy;

    // Fermi level (on relative energy scale; at the reference temperature)
    double m_FermiLevel;

    // Spatial conversion factor (absolute in nm = factor * relative)
    double m_SpatialFactor;

    // Energy axis values (in eV; zero = Fermi level at reference temperature)
    std::vector<double> m_EnergyValues;

    // DOS values (in 1/cm3eV)
    std::vector<double> m_DOSValues;

    // Relevant energy axis (cropped and scaled: -1 to 0)
    std::vector<double> m_RelEnergyValues;

    // Relevant DOS values (in 1/cm3eV)
    std::vector<double> m_RelDOSValues;

    // DOS random number distribution (normalized)
    std::piecewise_linear_distribution<double> m_DOSDistribution;
};

} // MC namespace
#endif  // TPiecewiseLinearDOS_H_
