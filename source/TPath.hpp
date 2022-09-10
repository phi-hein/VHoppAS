#ifndef TPath_H_
#define TPath_H_

#include <cstdint>

namespace MC
{

class TPath 
{
public:
    // Hopping path (index of the destination state)
    std::uint32_t m_StateID;

    // Reverse path (index of the reverse path in the path list of the destination state)
    std::uint32_t m_ReverseID;

    // Hopping time (inverse hopping rate)
    double m_Time;

    // Hopping counter (for statistics)
    mutable std::uint64_t m_HopCount;

    // Oscillation counter (for statistics)
    mutable std::uint64_t m_OscHopCount;

    // (disabled) Default constructor
	TPath () = delete;

	// Initialization constructor
	TPath (std::uint32_t state_id, std::uint32_t reverse_id, double time);
	
	// (default) Move constructor
	TPath (TPath&& other) = default;
	
	// (default) Move assignment operator
	TPath& operator= (TPath&& other) = default;
};

} // MC namespace
#endif  // TPath_H_