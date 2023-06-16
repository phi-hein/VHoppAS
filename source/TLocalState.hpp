#ifndef TLocalState_H_
#define TLocalState_H_

#include <vector>
#include <string>
#include <cstdint>

#include "TPath.hpp"

namespace MC
{

class TLocalState 
{
public:
	// Position vector (relative: 0 to 1)
	const double m_Pos_x, m_Pos_y, m_Pos_z;

    // Energy (relative: -1 to 0)
    const double m_Energy;

	// Occupying electron (index in electron list; equals mI_StateCount if empty)
	mutable std::uint32_t m_ElectronID;

	// Time of last occupation change (incoming or outgoing hop; in s)
	mutable double m_LastHopTime;

	// Occupied timespan (in s; no valid value until end of simulation)
	mutable double m_OccTime;

    // Hopping paths
    std::vector<TPath> m_Paths;

	// (disabled) Default constructor
	TLocalState () = delete;

    // Initialization constructor (position, energy and electron index)
	TLocalState (double x, double y, double z, double energy, std::uint32_t electron_id);
	
	// Move constructor
	TLocalState (TLocalState&& other);
	
	// (disabled) Move assignment operator
	TLocalState& operator= (TLocalState&& other) = delete;

    // (disabled) Copy assignment operator
	TLocalState& operator= (const TLocalState& other) = delete;

	// String output (x, y, z, energy in relative units)
	std::string Write(const std::string& separator = ", ") const;

	// String output (x, y, z, energy in absolute units)
	std::string Write(double spatial_scale, double max_energy, double energy_scale, const std::string& separator = ", ") const;
};

} // MC namespace
#endif  // TLocalState_H_