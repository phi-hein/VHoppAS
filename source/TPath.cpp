#include "TPath.hpp"

// Initialization constructor
MC::TPath::TPath(std::uint32_t state_id, std::uint32_t reverse_id, double time)
    : m_StateID(state_id), m_ReverseID(reverse_id), m_Time(time), m_HopCount(0U), m_OscHopCount(0U)
{

}