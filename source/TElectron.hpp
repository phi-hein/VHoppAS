#ifndef TElectron_H_
#define TElectron_H_

#include <vector>
#include <cstdint>

namespace MC
{

class TElectron 
{
public:
	// Time step for fastest transition (in s since m_LastHopTime)
	double m_MinTime;

	// Time step for next fastest transition (in s since m_LastHopTime)
	double m_NextMinTime;

	// Simulation time when this electron made it's last hop
	double m_LastHopTime;

	// Index of fastest transition (lowest of m_RandomTimes)
	std::uint32_t m_MinIndex;

	// Index of next fastest transition (2nd lowest of m_RandomTimes)
	std::uint32_t m_NextMinIndex;

	// Target state index of the fastest transition
	std::uint32_t m_MinStateID;

	// Randomized hopping times for all paths 
	// (in s since m_LastHopTime; in the same order as the paths in the occupied local state)
	std::vector<double> m_RandomTimes;

	// State index of this electron's current location
	std::uint32_t m_CurrentStateID;

	// State index of this electron's last location
	std::uint32_t m_LastStateID;
	
	// Displacement vector
	double m_Disp_x, m_Disp_y, m_Disp_z;
	
	// Hopping count
	std::uint64_t m_HopCount;

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