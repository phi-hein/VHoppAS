#ifndef TElectron_H_
#define TElectron_H_

#include <vector>
#include <cstdint>

namespace MC
{

class TElectron 
{
public:
	// Time step for fastest transition (already including acc_time_decay for fast sorting of electrons)
	double m_MinTime;

	// Accumulated time decay (time elapsed since the last update of this electrons' transitions)
	double m_AccTimeDecay;

	// Index of fastest transition
	std::uint32_t m_MinIndex;

	// Randomized hopping times for all paths (in the same order as the paths in the occupied local state)
	std::vector<double> m_RandomTimes;

	// State index of this electron's current location
	std::uint32_t m_CurrentStateID;
	
	// Displacement vector
	double m_Disp_x, m_Disp_y, m_Disp_z;
	
	// Hopping count
	std::uint64_t m_HopCount;

	// State index of this electron's last location
	std::uint32_t m_LastStateID;

	// Oscillatory hop count (back-and-forth = 2 hops)
	std::uint64_t m_OscHopCount;

	// Simulation time when this electron made it's first hop
	double m_FirstHopTime;

	// (disabled) Default constructor
	TElectron () = delete;

	// Initialization constructor
	TElectron (std::uint32_t state_id);
	
	// Move constructor
	TElectron (TElectron&& other);
	
	// Move assignment operator
	TElectron& operator= (TElectron&& other);
};

} // MC namespace
#endif  // TElectron_H_