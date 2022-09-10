#include "TLocalState.hpp"

#include <utility>
#include <sstream>
#include <iomanip>

// Initialization constructor (position, energy and electron index)
MC::TLocalState::TLocalState(double x, double y, double z, double energy, std::uint32_t electron_id)
    : m_Pos_x(x), m_Pos_y(y), m_Pos_z(z), m_Energy(energy), m_ElectronID(electron_id),
    m_LastHopTime(0.0), m_OccTime(0.0)
{
	
}

// Move constructor
MC::TLocalState::TLocalState(TLocalState&& other)
	: m_Pos_x(other.m_Pos_x), m_Pos_y(other.m_Pos_y), m_Pos_z(other.m_Pos_z), m_Energy(other.m_Energy),
	m_ElectronID(other.m_ElectronID), m_Paths(std::move(other.m_Paths)), 
    m_LastHopTime(other.m_LastHopTime), m_OccTime(other.m_OccTime)
{
	
}

// String output (x, y, z, energy in relative units)
std::string MC::TLocalState::Write(const std::string& separator) const
{
    std::ostringstream osstr;
    osstr << std::setprecision(8) << m_Pos_x << separator << m_Pos_y << separator << m_Pos_z <<
        separator << m_Energy;
    return osstr.str();
}

// String output (x, y, z, energy in absolute units)
std::string MC::TLocalState::Write(double spatial_scale, double max_energy, double energy_scale, const std::string& separator) const
{
    std::ostringstream osstr;
    osstr << std::setprecision(8) << spatial_scale*m_Pos_x << separator << 
        spatial_scale*m_Pos_y << separator << spatial_scale*m_Pos_z << separator << 
        max_energy + energy_scale*m_Energy;
    return osstr.str();
}