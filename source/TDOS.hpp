#ifndef TDOS_H_
#define TDOS_H_

#include <random>
#include <memory>

#include "TParamSet.hpp"

namespace MC
{

class TDOS 
{
public:
    // Destructor (necessary for std::unique_ptr<TDOS> to call derived destructors)
    virtual ~TDOS() { }

    // Get copy of the internal object
    virtual std::unique_ptr<TDOS> Copy() = 0;

    // Validate additional parameters for DOS
    virtual void ValidateParameters(const TParamSet& pset) const = 0;

    // Set additional parameters for DOS
    virtual void ApplyParameters(const TParamSet& pset, std::string& output) = 0;

    // Remove applied parameters (except the raw DOS data)
    virtual void DeleteParameters() = 0;

    // Returns whether DOS is properly defined
    virtual bool HasDOS() const = 0;

    // Returns whether DOS and its internal parameters (e.g. scaling) are properly defined
    virtual bool IsReady() const = 0;

	// Scaling factor for spatial differences (in nm): absolute = factor*relative
    virtual double GetSpatialFactor() const = 0;

    // Scaling factor for energy differences (in eV): absolute = factor*relative
    virtual double GetEnergyFactor() const = 0;

    // Calculate F-D occupation probability based on state energy and kBT (both in relative energy)
    virtual double GetOccupationProbability(double energy, double kBT) const = 0;
    
    // Generate random state energy (relative: -1 to 0) according to DOS
	virtual double GetRandomEnergy(std::mt19937_64 &rng) = 0;

    // Calculate expected number of electrons based on selected DOS range and Fermi level (kBT in relative energy)
    virtual std::uint32_t GetElectronCount(double kBT) const = 0;

    // Calculate expected effective carrier density (in 1/cm3) based on chemical potential (in eV) and temperature (in K)
    virtual double GetEffCarrierDensity(double chem_pot, double temp) const = 0;

    // Calculate expected partial entropy (in eV/K) based on chemical potential (in eV) and temperature (in K)
    virtual double GetPartialEntropy(double chem_pot, double temp) const = 0;

    // Calculate Fermi level (in eV) adjusted to another temperature (in K)
    virtual double GetAdjustedFermiLevel(double newT) const = 0;
};

} // MC namespace
#endif  // TDOS_H_
