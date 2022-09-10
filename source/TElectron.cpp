#include "TElectron.hpp"

#include <utility>

// Initialization constructor
MC::TElectron::TElectron(std::uint32_t state_id)
	: m_CurrentStateID(state_id), m_Disp_x(0.0), m_Disp_y(0.0), m_Disp_z(0.0), m_HopCount(0), 
	m_AccTimeDecay(0.0), m_MinIndex(0), m_MinTime(0.0), m_LastStateID(state_id), m_OscHopCount(0), 
	m_FirstHopTime(0.0)
{
	
}

// Move constructor
MC::TElectron::TElectron(TElectron&& other)
	: m_CurrentStateID(other.m_CurrentStateID), m_Disp_x(other.m_Disp_x), m_Disp_y(other.m_Disp_y), 
	m_Disp_z(other.m_Disp_z), m_HopCount(other.m_HopCount), m_AccTimeDecay(other.m_AccTimeDecay), 
	m_MinIndex(other.m_MinIndex), m_MinTime(other.m_MinTime), m_LastStateID(other.m_LastStateID),
	m_OscHopCount(other.m_OscHopCount), m_FirstHopTime(other.m_FirstHopTime),
	m_RandomTimes(std::move(other.m_RandomTimes))
{
	
}

// Move assignment operator
MC::TElectron& MC::TElectron::operator=(TElectron&& other)
{
	// Move assign fields of other to this
	m_CurrentStateID = other.m_CurrentStateID;
	m_Disp_x = other.m_Disp_x;
	m_Disp_y = other.m_Disp_y;
	m_Disp_z = other.m_Disp_z;
	m_HopCount = other.m_HopCount;
	m_AccTimeDecay = other.m_AccTimeDecay;
	m_MinIndex = other.m_MinIndex;
	m_MinTime = other.m_MinTime;
	m_LastStateID = other.m_LastStateID;
	m_OscHopCount = other.m_OscHopCount;
	m_FirstHopTime = other.m_FirstHopTime;
	m_RandomTimes = std::move(other.m_RandomTimes);
	
	return *this;
}
