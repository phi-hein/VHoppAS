#include "TResult.hpp"

#include <sstream>
#include <regex>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "Constants.hpp"
#include "GlobalFunctions.hpp"
#include "CustomExceptions.hpp"

// Result descriptors and units (without white-spaces)
const std::array<std::string,2> MC::TResult::s_RepID = {"RepID",""};
const std::array<std::string,2> MC::TResult::s_DriftConductivity = {"DriftCond","S/m"};
const std::array<std::string,2> MC::TResult::s_DriftMobility = {"DriftMob","cm2/Vs"};
const std::array<std::string,2> MC::TResult::s_DiffusionCoefficient = {"DiffCoeff","cm2/s"};
const std::array<std::string,2> MC::TResult::s_DiffusionCoefficientParallel = {"DiffCoeffP","cm2/s"};
const std::array<std::string,2> MC::TResult::s_DiffusionCoefficientTransverse = {"DiffCoeffT","cm2/s"};
const std::array<std::string,2> MC::TResult::s_HavenRatio = {"HavenRatio",""};
const std::array<std::string,2> MC::TResult::s_HavenRatioParallel = {"HavenRatioP",""};
const std::array<std::string,2> MC::TResult::s_HavenRatioTransverse = {"HavenRatioT",""};
const std::array<std::string,2> MC::TResult::s_PartialEntropy = {"PartialEntropy","eV/K"};
const std::array<std::string,2> MC::TResult::s_EffChemPot = {"EffChemPot","eV"};
const std::array<std::string,2> MC::TResult::s_EffTemp = {"EffTemp","K"};
const std::array<std::string,2> MC::TResult::s_EffCarriers = {"EffCarriers",""};
const std::array<std::string,2> MC::TResult::s_EffCarrierDensity = {"EffCarrierDensity","1/cm3"};
const std::array<std::string,2> MC::TResult::s_ElectronCount = {"Electrons",""};
const std::array<std::string,2> MC::TResult::s_MobileElectrons = {"MobileElectrons",""};
const std::array<std::string,2> MC::TResult::s_ZeroHopElectrons = {"ZeroHopElectrons",""};
const std::array<std::string,2> MC::TResult::s_OscElectrons = {"OscElectrons",""};
const std::array<std::string,2> MC::TResult::s_TotalTime = {"Time","s"};
const std::array<std::string,2> MC::TResult::s_MeanDisp = {"MeanDisp","nm"};
const std::array<std::string,2> MC::TResult::s_MeanDisp_x = {"MeanDispx","nm"};
const std::array<std::string,2> MC::TResult::s_MeanDisp_y = {"MeanDispy","nm"};
const std::array<std::string,2> MC::TResult::s_MeanDisp_z = {"MeanDispz","nm"};
const std::array<std::string,2> MC::TResult::s_MeanDispVariance_x = {"MeanDispVarx","nm2"};
const std::array<std::string,2> MC::TResult::s_MeanDispVariance_y = {"MeanDispVary","nm2"};
const std::array<std::string,2> MC::TResult::s_MeanDispVariance_z = {"MeanDispVarz","nm2"};
const std::array<std::string,2> MC::TResult::s_MeanSquaredDisp_x = {"MeanSqDispx","nm2"};
const std::array<std::string,2> MC::TResult::s_MeanSquaredDisp_y = {"MeanSqDispy","nm2"};
const std::array<std::string,2> MC::TResult::s_MeanSquaredDisp_z = {"MeanSqDispz","nm2"};	
const std::array<std::string,2> MC::TResult::s_NonOscHops = {"NonOscHops",""};
const std::array<std::string,2> MC::TResult::s_NonOscHopRatio = {"NonOscRatio",""};
const std::array<std::string,2> MC::TResult::s_InXDirNonOscRatio = {"NonOscXDirRatio",""};
const std::array<std::string,2> MC::TResult::s_MeanFieldContribution = {"MeanFieldContrib","eV"};
const std::array<std::string,2> MC::TResult::s_TotalEnergy = {"TotalEnergy","eV"};
const std::array<std::string,2> MC::TResult::s_CellSize = {"CellSize","nm"};
const std::array<std::string,2> MC::TResult::s_MaxPathDist = {"MaxPathDist","nm"};
const std::array<std::string,2> MC::TResult::s_MaxUsedPathDist = {"MaxUsedPathDist","nm"};
const std::array<std::string,2> MC::TResult::s_MaxPathEdiff = {"MaxPathEdiff","eV"};
const std::array<std::string,2> MC::TResult::s_MaxUsedPathEdiff = {"MaxUsedPathEdiff","eV"};
const std::array<std::string,2> MC::TResult::s_MaxPathCount = {"MaxPaths",""};
const std::array<std::string,2> MC::TResult::s_MeanPathCount = {"MeanPaths",""};
const std::array<std::string,2> MC::TResult::s_MinUsedStateEnergy = {"MinUsedStateE","eV"};
const std::array<std::string,2> MC::TResult::s_MaxUsedStateEnergy = {"MaxUsedStateE","eV"};

// Default constructor
MC::TResult::TResult()
    : m_TotalHops(0), m_RepID(0), 
    m_DriftConductivity(0.0), m_DriftMobility(0.0), m_DiffusionCoefficient(0.0),
    m_DiffusionCoefficientParallel(0.0), m_DiffusionCoefficientTransverse(0.0), 
    m_HavenRatio(0.0), m_HavenRatioParallel(0.0), m_HavenRatioTransverse(0.0), m_PartialEntropy(0.0), 
    m_EffChemPot(0.0), m_EffTemp(0.0), m_EffCarriers(0.0), m_EffCarrierDensity(0.0), 
    m_ElectronCount(0), m_MobileElectrons(0), m_ZeroHopElectrons(0), m_OscElectrons(0), 
    m_TotalTime(0.0), 
    m_MeanDisp(0.0), m_MeanDisp_x(0.0), m_MeanDisp_y(0.0), m_MeanDisp_z(0.0), 
    m_MeanDispVariance_x(0.0), m_MeanDispVariance_y(0.0), m_MeanDispVariance_z(0.0), 
    m_MeanSquaredDisp_x(0.0), m_MeanSquaredDisp_y(0.0), m_MeanSquaredDisp_z(0.0),
    m_NonOscHops(0), m_NonOscHopRatio(0.0), m_InXDirNonOscRatio(0.0), m_MeanFieldContribution(0.0),
    m_TotalEnergy(0.0), m_CellSize(0.0), 
    m_MaxPathDist(0.0), m_MaxUsedPathDist(0.0), m_MaxPathEdiff(0.0), m_MaxUsedPathEdiff(0.0),
    m_MaxPathCount(0), m_MeanPathCount(0.0), m_MinUsedStateEnergy(0.0), m_MaxUsedStateEnergy(0.0)
{

}

// Copy constructor (ignores histograms)
std::unique_ptr<MC::TResult> MC::TResult::ValueCopy(const TResult* const result)
{
    if (result == nullptr) return nullptr;

    std::unique_ptr<TResult> temporary;
    try
    {
        temporary = std::make_unique<TResult>();
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot create copy of results.",__func__,e.what());
    }
    TResult& temp = *temporary;

    temp.m_TotalHops = result->m_TotalHops;
    temp.m_RepID = result->m_RepID;
    temp.m_DriftConductivity = result->m_DriftConductivity;
    temp.m_DriftMobility = result->m_DriftMobility;
    temp.m_DiffusionCoefficient = result->m_DiffusionCoefficient;
	temp.m_DiffusionCoefficientParallel = result->m_DiffusionCoefficientParallel;
	temp.m_DiffusionCoefficientTransverse = result->m_DiffusionCoefficientTransverse;
    temp.m_HavenRatio = result->m_HavenRatio;
    temp.m_HavenRatioParallel = result->m_HavenRatioParallel;
    temp.m_HavenRatioTransverse = result->m_HavenRatioTransverse;
    temp.m_PartialEntropy = result->m_PartialEntropy;
    temp.m_EffChemPot = result->m_EffChemPot;
    temp.m_EffTemp = result->m_EffTemp;
    temp.m_EffCarriers = result->m_EffCarriers;
	temp.m_EffCarrierDensity = result->m_EffCarrierDensity;
    temp.m_ElectronCount = result->m_ElectronCount;
    temp.m_MobileElectrons = result->m_MobileElectrons;
	temp.m_ZeroHopElectrons = result->m_ZeroHopElectrons;
    temp.m_OscElectrons = result->m_OscElectrons;
	temp.m_TotalTime = result->m_TotalTime;
    temp.m_MeanDisp = result->m_MeanDisp;
	temp.m_MeanDisp_x = result->m_MeanDisp_x;
	temp.m_MeanDisp_y = result->m_MeanDisp_y;
	temp.m_MeanDisp_z = result->m_MeanDisp_z;
	temp.m_MeanDispVariance_x = result->m_MeanDispVariance_x;
    temp.m_MeanDispVariance_y = result->m_MeanDispVariance_y;
    temp.m_MeanDispVariance_z = result->m_MeanDispVariance_z;
	temp.m_MeanSquaredDisp_x = result->m_MeanSquaredDisp_x;
	temp.m_MeanSquaredDisp_y = result->m_MeanSquaredDisp_y;
	temp.m_MeanSquaredDisp_z = result->m_MeanSquaredDisp_z;
    temp.m_NonOscHops = result->m_NonOscHops;
    temp.m_NonOscHopRatio = result->m_NonOscHopRatio;
    temp.m_InXDirNonOscRatio = result->m_InXDirNonOscRatio;
    temp.m_MeanFieldContribution = result->m_MeanFieldContribution;
    temp.m_TotalEnergy = result->m_TotalEnergy;
    temp.m_CellSize = result->m_CellSize;
	temp.m_MaxPathDist = result->m_MaxPathDist;
    temp.m_MaxUsedPathDist = result->m_MaxUsedPathDist;
    temp.m_MaxPathEdiff = result->m_MaxPathEdiff;
    temp.m_MaxUsedPathEdiff = result->m_MaxUsedPathEdiff;
	temp.m_MaxPathCount = result->m_MaxPathCount;
	temp.m_MeanPathCount = result->m_MeanPathCount;
    temp.m_MinUsedStateEnergy = result->m_MinUsedStateEnergy;
    temp.m_MaxUsedStateEnergy = result->m_MaxUsedStateEnergy;
    
    return temporary;
}

// Write results block
void MC::TResult::Write(std::ostream& o_str) const
{
    o_str << "<Results>" << std::endl;
    o_str << GF::CombineDescUnit(s_RepID) << " = " << m_RepID << std::endl;
    o_str << GF::CombineDescUnit(s_DriftConductivity) << " = " << m_DriftConductivity << std::endl;
    o_str << GF::CombineDescUnit(s_DriftMobility) << " = " << m_DriftMobility << std::endl;
    o_str << GF::CombineDescUnit(s_DiffusionCoefficient) << " = " << m_DiffusionCoefficient << std::endl;
    o_str << GF::CombineDescUnit(s_DiffusionCoefficientParallel) << " = " << m_DiffusionCoefficientParallel << std::endl;
    o_str << GF::CombineDescUnit(s_DiffusionCoefficientTransverse) << " = " << m_DiffusionCoefficientTransverse << std::endl;
    o_str << GF::CombineDescUnit(s_HavenRatio) << " = " << m_HavenRatio << std::endl;
    o_str << GF::CombineDescUnit(s_HavenRatioParallel) << " = " << m_HavenRatioParallel << std::endl;
    o_str << GF::CombineDescUnit(s_HavenRatioTransverse) << " = " << m_HavenRatioTransverse << std::endl;
    o_str << GF::CombineDescUnit(s_PartialEntropy) << " = " << m_PartialEntropy << std::endl;
    o_str << GF::CombineDescUnit(s_EffChemPot) << " = " << m_EffChemPot << std::endl;
    o_str << GF::CombineDescUnit(s_EffTemp) << " = " << m_EffTemp << std::endl;
    o_str << GF::CombineDescUnit(s_EffCarriers) << " = " << m_EffCarriers << std::endl;
    o_str << GF::CombineDescUnit(s_EffCarrierDensity) << " = " << m_EffCarrierDensity << std::endl;
    o_str << GF::CombineDescUnit(s_ElectronCount) << " = " << m_ElectronCount << std::endl;
    o_str << GF::CombineDescUnit(s_MobileElectrons) << " = " << m_MobileElectrons << std::endl;
    o_str << GF::CombineDescUnit(s_ZeroHopElectrons) << " = " << m_ZeroHopElectrons << std::endl;
    o_str << GF::CombineDescUnit(s_OscElectrons) << " = " << m_OscElectrons << std::endl;
    o_str << GF::CombineDescUnit(s_TotalTime) << " = " << m_TotalTime << std::endl;
    o_str << GF::CombineDescUnit(s_MeanDisp) << " = " << m_MeanDisp << std::endl;
    o_str << GF::CombineDescUnit(s_MeanDisp_x) << " = " << m_MeanDisp_x << std::endl;
    o_str << GF::CombineDescUnit(s_MeanDisp_y) << " = " << m_MeanDisp_y << std::endl;
    o_str << GF::CombineDescUnit(s_MeanDisp_z) << " = " << m_MeanDisp_z << std::endl;
    o_str << GF::CombineDescUnit(s_MeanDispVariance_x) << " = " << m_MeanDispVariance_x << std::endl;
    o_str << GF::CombineDescUnit(s_MeanDispVariance_y) << " = " << m_MeanDispVariance_y << std::endl;
    o_str << GF::CombineDescUnit(s_MeanDispVariance_z) << " = " << m_MeanDispVariance_z << std::endl;
    o_str << GF::CombineDescUnit(s_MeanSquaredDisp_x) << " = " << m_MeanSquaredDisp_x << std::endl;
    o_str << GF::CombineDescUnit(s_MeanSquaredDisp_y) << " = " << m_MeanSquaredDisp_y << std::endl;
    o_str << GF::CombineDescUnit(s_MeanSquaredDisp_z) << " = " << m_MeanSquaredDisp_z << std::endl;
    o_str << GF::CombineDescUnit(s_NonOscHops) << " = " << m_NonOscHops << std::endl;
    o_str << GF::CombineDescUnit(s_NonOscHopRatio) << " = " << m_NonOscHopRatio << std::endl;
    o_str << GF::CombineDescUnit(s_InXDirNonOscRatio) << " = " << m_InXDirNonOscRatio << std::endl;
    o_str << GF::CombineDescUnit(s_MeanFieldContribution) << " = " << m_MeanFieldContribution << std::endl;
    o_str << GF::CombineDescUnit(s_TotalEnergy) << " = " << m_TotalEnergy << std::endl;
    o_str << GF::CombineDescUnit(s_CellSize) << " = " << m_CellSize << std::endl;
    o_str << GF::CombineDescUnit(s_MaxPathDist) << " = " << m_MaxPathDist << std::endl;
    o_str << GF::CombineDescUnit(s_MaxUsedPathDist) << " = " << m_MaxUsedPathDist << std::endl;
    o_str << GF::CombineDescUnit(s_MaxPathEdiff) << " = " << m_MaxPathEdiff << std::endl;
    o_str << GF::CombineDescUnit(s_MaxUsedPathEdiff) << " = " << m_MaxUsedPathEdiff << std::endl;
    o_str << GF::CombineDescUnit(s_MaxPathCount) << " = " << m_MaxPathCount << std::endl;
    o_str << GF::CombineDescUnit(s_MeanPathCount) << " = " << m_MeanPathCount << std::endl;
    o_str << GF::CombineDescUnit(s_MinUsedStateEnergy) << " = " << m_MinUsedStateEnergy << std::endl;
    o_str << GF::CombineDescUnit(s_MaxUsedStateEnergy) << " = " << m_MaxUsedStateEnergy << std::endl;
    o_str << "</Results>" << std::endl;
    
    if (!m_EqProgress.empty())
    {
        o_str << std::endl;
        o_str << "<EquilibrationConvergence>" << std::endl;
        WriteConvergenceTable(o_str,true);
        o_str << "</EquilibrationConvergence>" << std::endl;
        o_str << "<!-- averages refer to effective carriers -->" << std::endl;
        o_str << "<!-- respective values are adjusted to the effective carrier density at the end of the complete simulation -->" << std::endl;
    }
    
    if (!m_Progress.empty())
    {
        o_str << std::endl;
        o_str << "<SimulationConvergence>" << std::endl;
        WriteConvergenceTable(o_str,false);
        o_str << "</SimulationConvergence>" << std::endl;
        o_str << "<!-- averages refer to effective carriers -->" << std::endl;
        o_str << "<!-- respective values are adjusted to the effective carrier density at the end of the complete simulation -->" << std::endl;
    }

    if (m_HStateEnergy.m_IsConfigured)
    {
        o_str << std::endl;
        m_HStateEnergy.Write(o_str,"StateEnergies");
    }

    if (m_HPathTime.m_IsConfigured)
    {
        o_str << std::endl;
        m_HPathTime.Write(o_str,"PathTimes");
    }

    if (m_HHopTime.m_IsConfigured)
    {
        o_str << std::endl;
        m_HHopTime.Write(o_str,"HopTimes");
    }

    if (m_HDisp.m_IsConfigured)
    {
        o_str << std::endl;
        m_HDisp.Write(o_str,"Displacements");
    }

    if (m_HDisp_x.m_IsConfigured)
    {
        o_str << std::endl;
        m_HDisp_x.Write(o_str,"xDisplacements");
    }

    if (m_HDisp_y.m_IsConfigured)
    {
        o_str << std::endl;
        m_HDisp_y.Write(o_str,"yDisplacements");
    }

    if (m_HDisp_z.m_IsConfigured)
    {
        o_str << std::endl;
        m_HDisp_z.Write(o_str,"zDisplacements");
    }

    if (m_HHopCount.m_IsConfigured)
    {
        o_str << std::endl;
        m_HHopCount.Write(o_str,"HopCounts");
    }

    if (m_HOscHopCount.m_IsConfigured)
    {
        o_str << std::endl;
        m_HOscHopCount.Write(o_str,"OscHopCounts");
    }

    if (m_HNonOscHopCount.m_IsConfigured)
    {
        o_str << std::endl;
        m_HNonOscHopCount.Write(o_str,"NonOscHopCounts");
    }

    if (m_HStateEnergyDifference.m_IsConfigured)
    {
        o_str << std::endl;
        m_HStateEnergyDifference.Write(o_str,"StateEnergyDifferences");
    }

    if (m_HFieldEnergyContribution.m_IsConfigured)
    {
        o_str << std::endl;
        m_HFieldEnergyContribution.Write(o_str,"FieldEnergyContributions");
    }

    if (m_HDistance.m_IsConfigured)
    {
        o_str << std::endl;
        m_HDistance.Write(o_str,"Distances");
    }

    if (m_HNextTime.m_IsConfigured)
    {
        o_str << std::endl;
        m_HNextTime.Write(o_str,"NextTransitionTimes");
    }

    if (m_HRelOccTime.m_IsConfigured)
    {
        o_str << std::endl;
        m_HRelOccTime.Write(o_str,"RelativeOccupationTimes");
    }

    o_str << std::endl;
}

// Write table header (as individual elements)
std::vector<std::array<std::string,2>> MC::TResult::WriteTableHeader()
{
    std::vector<std::array<std::string,2>> header;
    header.push_back(s_RepID);
    header.push_back(s_DriftConductivity);
    header.push_back(s_DriftMobility);
    header.push_back(s_DiffusionCoefficient);
    header.push_back(s_DiffusionCoefficientParallel);
    header.push_back(s_DiffusionCoefficientTransverse);
    header.push_back(s_HavenRatio);
    header.push_back(s_HavenRatioParallel);
    header.push_back(s_HavenRatioTransverse);
    header.push_back(s_PartialEntropy);
    header.push_back(s_EffChemPot);
    header.push_back(s_EffTemp);
    header.push_back(s_EffCarriers);
    header.push_back(s_EffCarrierDensity);
    header.push_back(s_ElectronCount);
    header.push_back(s_MobileElectrons);
    header.push_back(s_ZeroHopElectrons);
    header.push_back(s_OscElectrons);
    header.push_back(s_TotalTime);
    header.push_back(s_MeanDisp);
    header.push_back(s_MeanDisp_x);
    header.push_back(s_MeanDisp_y);
    header.push_back(s_MeanDisp_z);
    header.push_back(s_MeanDispVariance_x);
    header.push_back(s_MeanDispVariance_y);
    header.push_back(s_MeanDispVariance_z);
    header.push_back(s_MeanSquaredDisp_x);
    header.push_back(s_MeanSquaredDisp_y);
    header.push_back(s_MeanSquaredDisp_z);
    header.push_back(s_NonOscHops);
    header.push_back(s_NonOscHopRatio);
    header.push_back(s_InXDirNonOscRatio);
    header.push_back(s_MeanFieldContribution);
    header.push_back(s_TotalEnergy);
    header.push_back(s_CellSize);
    header.push_back(s_MaxPathDist);
    header.push_back(s_MaxUsedPathDist);
    header.push_back(s_MaxPathEdiff);
    header.push_back(s_MaxUsedPathEdiff);
    header.push_back(s_MaxPathCount);
    header.push_back(s_MeanPathCount);
    header.push_back(s_MinUsedStateEnergy);
    header.push_back(s_MaxUsedStateEnergy);

    for (auto& item : header)
    {
        if (item[1].empty()) item[1] = "_";
    }
    
    return header;
}

// Write table line (as individual elements)
std::vector<std::string> MC::TResult::WriteTableLine() const
{
    std::vector<std::string> line;
    std::stringstream sstr;
    
    sstr << m_RepID;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_DriftConductivity;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_DriftMobility;
    line.push_back(sstr.str());
    sstr.str("");
    
    sstr << m_DiffusionCoefficient;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_DiffusionCoefficientParallel;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_DiffusionCoefficientTransverse;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_HavenRatio;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_HavenRatioParallel;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_HavenRatioTransverse;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_PartialEntropy;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_EffChemPot;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_EffTemp;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_EffCarriers;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_EffCarrierDensity;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_ElectronCount;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MobileElectrons;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_ZeroHopElectrons;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_OscElectrons;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_TotalTime;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanDisp;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanDisp_x;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanDisp_y;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanDisp_z;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanDispVariance_x;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanDispVariance_y;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanDispVariance_z;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanSquaredDisp_x;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanSquaredDisp_y;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanSquaredDisp_z;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_NonOscHops;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_NonOscHopRatio;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_InXDirNonOscRatio;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanFieldContribution;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_TotalEnergy;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_CellSize;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MaxPathDist;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MaxUsedPathDist;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MaxPathEdiff;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MaxUsedPathEdiff;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MaxPathCount;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MeanPathCount;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MinUsedStateEnergy;
    line.push_back(sstr.str());
    sstr.str("");

    sstr << m_MaxUsedStateEnergy;
    line.push_back(sstr.str());
    sstr.str("");
    
    return line;
}

// Read results (true = successful read)
bool MC::TResult::Read(const std::string& str, bool raise_errors)
{
    // Define lambda functions to get results
    auto get_dbl = [&] (const std::array<std::string,2>& desc, double& val) -> bool
    {
        std::smatch match;
        if ((std::regex_search(str, match, 
            std::regex(GF::DescRegex(desc) + "\\s*=\\s*(" + Constant::dblex + ")"))) &&
            (GF::StringToDouble(match[1],val)))
        {
            return true;
        } else {
            if (raise_errors) throw EX::TInvalidInput("Result \"" + desc[0] + "\" is missing.");
            return false;
        }
    };
    auto get_uint32 = [&] (const std::array<std::string,2>& desc, std::uint32_t& val) -> bool
    {
        std::smatch match;
        if ((std::regex_search(str, match, 
            std::regex(GF::DescRegex(desc) + "\\s*=\\s*(" + Constant::uintex + ")"))) &&
            (GF::StringToUInt32(match[1],val)))
        {
            return true;
        } else {
            if (raise_errors) throw EX::TInvalidInput("Result \"" + desc[0] + "\" is missing.");
            return false;
        }
    };
    auto get_uint64 = [&] (const std::array<std::string,2>& desc, std::uint64_t& val) -> bool
    {
        std::smatch match;
        if ((std::regex_search(str, match, 
            std::regex(GF::DescRegex(desc) + "\\s*=\\s*(" + Constant::uintex + ")"))) &&
            (GF::StringToUInt64(match[1],val)))
        {
            return true;
        } else {
            if (raise_errors) throw EX::TInvalidInput("Result \"" + desc[0] + "\" is missing.");
            return false;
        }
    };
    
    // Retrieve results
    m_TotalHops = 0;
    if (get_uint32(s_RepID,m_RepID) == false) return false;
    if (m_RepID == 0) return false;
    if (get_dbl(s_DriftConductivity,m_DriftConductivity) == false) return false;
    if (get_dbl(s_DriftMobility,m_DriftMobility) == false) return false;
    if (get_dbl(s_DiffusionCoefficient,m_DiffusionCoefficient) == false) return false;
    if (get_dbl(s_DiffusionCoefficientParallel,m_DiffusionCoefficientParallel) == false) return false;
    if (get_dbl(s_DiffusionCoefficientTransverse,m_DiffusionCoefficientTransverse) == false) return false;
	if (get_dbl(s_HavenRatio,m_HavenRatio) == false) return false;
    if (get_dbl(s_HavenRatioParallel,m_HavenRatioParallel) == false) return false;
    if (get_dbl(s_HavenRatioTransverse,m_HavenRatioTransverse) == false) return false;
    if (get_dbl(s_PartialEntropy,m_PartialEntropy) == false) return false;
    if (get_dbl(s_EffChemPot,m_EffChemPot) == false) return false;
    if (get_dbl(s_EffTemp,m_EffTemp) == false) return false;
    if (get_dbl(s_EffCarriers,m_EffCarriers) == false) return false;
    if (get_dbl(s_EffCarrierDensity,m_EffCarrierDensity) == false) return false;
    if (get_uint32(s_ElectronCount,m_ElectronCount) == false) return false;
    if (get_uint32(s_MobileElectrons,m_MobileElectrons) == false) return false;
	if (get_uint32(s_ZeroHopElectrons,m_ZeroHopElectrons) == false) return false;
    if (get_uint32(s_OscElectrons,m_OscElectrons) == false) return false;
    if (get_dbl(s_TotalTime,m_TotalTime) == false) return false;
    if (get_dbl(s_MeanDisp,m_MeanDisp) == false) return false;
	if (get_dbl(s_MeanDisp_x,m_MeanDisp_x) == false) return false;
	if (get_dbl(s_MeanDisp_y,m_MeanDisp_y) == false) return false;
	if (get_dbl(s_MeanDisp_z,m_MeanDisp_z) == false) return false;
    if (get_dbl(s_MeanDispVariance_x,m_MeanDispVariance_x) == false) return false;
    if (get_dbl(s_MeanDispVariance_y,m_MeanDispVariance_y) == false) return false;
    if (get_dbl(s_MeanDispVariance_z,m_MeanDispVariance_z) == false) return false;
	if (get_dbl(s_MeanSquaredDisp_x,m_MeanSquaredDisp_x) == false) return false;
	if (get_dbl(s_MeanSquaredDisp_y,m_MeanSquaredDisp_y) == false) return false;
	if (get_dbl(s_MeanSquaredDisp_z,m_MeanSquaredDisp_z) == false) return false;
    if (get_uint64(s_NonOscHops,m_NonOscHops) == false) return false;
    if (get_dbl(s_NonOscHopRatio,m_NonOscHopRatio) == false) return false;
    if (get_dbl(s_InXDirNonOscRatio,m_InXDirNonOscRatio) == false) return false;
    if (get_dbl(s_MeanFieldContribution,m_MeanFieldContribution) == false) return false;
    if (get_dbl(s_TotalEnergy,m_TotalEnergy) == false) return false;
    if (get_dbl(s_CellSize,m_CellSize) == false) return false;
    if (get_dbl(s_MaxPathDist,m_MaxPathDist) == false) return false;
    if (get_dbl(s_MaxUsedPathDist,m_MaxUsedPathDist) == false) return false;
    if (get_dbl(s_MaxPathEdiff,m_MaxPathEdiff) == false) return false;
    if (get_dbl(s_MaxUsedPathEdiff,m_MaxUsedPathEdiff) == false) return false;
	if (get_uint32(s_MaxPathCount,m_MaxPathCount) == false) return false;
	if (get_dbl(s_MeanPathCount,m_MeanPathCount) == false) return false;
    if (get_dbl(s_MinUsedStateEnergy,m_MinUsedStateEnergy) == false) return false;
    if (get_dbl(s_MaxUsedStateEnergy,m_MaxUsedStateEnergy) == false) return false;
	
    return true;
}

// Write analysis for results of multiple simulations
void MC::TResult::WriteMultiAnalysis(std::ostream& o_str, const std::vector<std::unique_ptr<const TResult>>& results)
{
    o_str << "Combined analysis: mean +/- stddev [min, max]" << std::endl;

    // Search index of the first defined result object
    auto first_it = std::find_if(results.begin(),results.end(),[](const std::unique_ptr<const TResult>& ptr){return ptr.get() != nullptr;});
    if (first_it == results.end())
    {
        o_str << "Finished simulations: 0" << std::endl;
        return;
    }
    std::size_t result_counter = std::count_if(results.begin(),results.end(),[](const std::unique_ptr<const TResult>& ptr){return ptr.get() != nullptr;});
    o_str << "Finished simulations: " << result_counter << std::endl;
    if (result_counter <= 1)
    {
        return;
    }

    // Define lambda functions to calculate mean, stddev, min and max for a specified property
    auto write_dbl_analysis = [&] (const std::array<std::string,2>& desc,double TResult::* prop)
    {
        double min_val = first_it->get()->*prop;
        double max_val = first_it->get()->*prop;
        double stddev_val = 0.0;
        double mean_val = 0.0;
        double prev_mean_val = 0.0;
        std::size_t counter = 0;
        for (auto it = first_it; it != results.end(); it++)
        {
            if (*it)
            {
                if (it->get()->*prop < min_val) min_val = it->get()->*prop;
                if (it->get()->*prop > max_val) max_val = it->get()->*prop;
                counter++;
                prev_mean_val = mean_val;
                mean_val += (it->get()->*prop - mean_val)/static_cast<double>(counter);
                stddev_val += (it->get()->*prop - mean_val) * (it->get()->*prop - prev_mean_val);
            }
        }
        if ((counter > 1) && (stddev_val != 0.0))
        {
            stddev_val = sqrt(fabs(stddev_val) / static_cast<double>(counter - 1));
        } else stddev_val = 0.0;
        o_str << GF::CombineDescUnit(desc) << " = " << mean_val << " +/- " << stddev_val << " [" << min_val << ", " << max_val << "]" << std::endl;
    };
    auto write_uint32_analysis = [&] (const std::array<std::string,2>& desc, std::uint32_t TResult::* prop)
    {
        std::uint32_t min_val = first_it->get()->*prop;
        std::uint32_t max_val = first_it->get()->*prop;
        double mean_val = 0.0;
        double stddev_val = 0.0;
        double prev_mean_val = 0.0;
        std::size_t counter = 0;
        for (auto it = first_it; it != results.end(); it++)
        {
            if (*it)
            {
                if (it->get()->*prop < min_val) min_val = it->get()->*prop;
                if (it->get()->*prop > max_val) max_val = it->get()->*prop;
                counter++;
                prev_mean_val = mean_val;
                mean_val += (static_cast<double>(it->get()->*prop) - mean_val)/static_cast<double>(counter);
                stddev_val += (static_cast<double>(it->get()->*prop) - mean_val) * 
                    (static_cast<double>(it->get()->*prop) - prev_mean_val);
            }
        }
        if ((counter > 1) && (stddev_val != 0.0))
        {
            stddev_val = sqrt(fabs(stddev_val) / static_cast<double>(counter - 1));
        } else stddev_val = 0.0;
        o_str << GF::CombineDescUnit(desc) << " = " << mean_val << " +/- " << stddev_val << " [" << min_val << ", " << max_val << "]" << std::endl;
    };
    auto write_uint64_analysis = [&] (const std::array<std::string,2>& desc, std::uint64_t TResult::* prop)
    {
        std::uint64_t min_val = first_it->get()->*prop;
        std::uint64_t max_val = first_it->get()->*prop;
        double mean_val = 0.0;
        double stddev_val = 0.0;
        double prev_mean_val = 0.0;
        std::size_t counter = 0;
        for (auto it = first_it; it != results.end(); it++)
        {
            if (*it)
            {
                if (it->get()->*prop < min_val) min_val = it->get()->*prop;
                if (it->get()->*prop > max_val) max_val = it->get()->*prop;
                counter++;
                prev_mean_val = mean_val;
                mean_val += (static_cast<double>(it->get()->*prop) - mean_val)/static_cast<double>(counter);
                stddev_val += (static_cast<double>(it->get()->*prop) - mean_val) * 
                    (static_cast<double>(it->get()->*prop) - prev_mean_val);
            }
        }
        if ((counter > 1) && (stddev_val != 0.0))
        {
            stddev_val = sqrt(fabs(stddev_val) / static_cast<double>(counter - 1));
        } else stddev_val = 0.0;
        o_str << GF::CombineDescUnit(desc) << " = " << mean_val << " +/- " << stddev_val << " [" << min_val << ", " << max_val << "]" << std::endl;
    };

    // Analysis
    write_dbl_analysis(s_DriftConductivity,&TResult::m_DriftConductivity);
    write_dbl_analysis(s_DriftMobility,&TResult::m_DriftMobility);
    write_dbl_analysis(s_DiffusionCoefficient,&TResult::m_DiffusionCoefficient);
    write_dbl_analysis(s_DiffusionCoefficientParallel,&TResult::m_DiffusionCoefficientParallel);
    write_dbl_analysis(s_DiffusionCoefficientTransverse,&TResult::m_DiffusionCoefficientTransverse);
    write_dbl_analysis(s_HavenRatio,&TResult::m_HavenRatio);
    write_dbl_analysis(s_HavenRatioParallel,&TResult::m_HavenRatioParallel);
    write_dbl_analysis(s_HavenRatioTransverse,&TResult::m_HavenRatioTransverse);
    write_dbl_analysis(s_PartialEntropy,&TResult::m_PartialEntropy);
    write_dbl_analysis(s_EffChemPot,&TResult::m_EffChemPot);
    write_dbl_analysis(s_EffTemp,&TResult::m_EffTemp);
    write_dbl_analysis(s_EffCarriers,&TResult::m_EffCarriers);
    write_dbl_analysis(s_EffCarrierDensity,&TResult::m_EffCarrierDensity);
    write_uint32_analysis(s_ElectronCount,&TResult::m_ElectronCount);
    write_uint32_analysis(s_MobileElectrons,&TResult::m_MobileElectrons);
	write_uint32_analysis(s_ZeroHopElectrons,&TResult::m_ZeroHopElectrons);
    write_uint32_analysis(s_OscElectrons,&TResult::m_OscElectrons);
	write_dbl_analysis(s_TotalTime,&TResult::m_TotalTime);
    write_dbl_analysis(s_MeanDisp,&TResult::m_MeanDisp);
	write_dbl_analysis(s_MeanDisp_x,&TResult::m_MeanDisp_x);
	write_dbl_analysis(s_MeanDisp_y,&TResult::m_MeanDisp_y);
	write_dbl_analysis(s_MeanDisp_z,&TResult::m_MeanDisp_z);
    write_dbl_analysis(s_MeanDispVariance_x,&TResult::m_MeanDispVariance_x);
    write_dbl_analysis(s_MeanDispVariance_y,&TResult::m_MeanDispVariance_y);
    write_dbl_analysis(s_MeanDispVariance_z,&TResult::m_MeanDispVariance_z);
	write_dbl_analysis(s_MeanSquaredDisp_x,&TResult::m_MeanSquaredDisp_x);
	write_dbl_analysis(s_MeanSquaredDisp_y,&TResult::m_MeanSquaredDisp_y);
	write_dbl_analysis(s_MeanSquaredDisp_z,&TResult::m_MeanSquaredDisp_z);
    write_uint64_analysis(s_NonOscHops,&TResult::m_NonOscHops);
    write_dbl_analysis(s_NonOscHopRatio,&TResult::m_NonOscHopRatio);
    write_dbl_analysis(s_InXDirNonOscRatio,&TResult::m_InXDirNonOscRatio);
    write_dbl_analysis(s_MeanFieldContribution,&TResult::m_MeanFieldContribution);
    write_dbl_analysis(s_TotalEnergy,&TResult::m_TotalEnergy);
    write_dbl_analysis(s_CellSize,&TResult::m_CellSize);
    write_dbl_analysis(s_MaxPathDist,&TResult::m_MaxPathDist);
    write_dbl_analysis(s_MaxUsedPathDist,&TResult::m_MaxUsedPathDist);
    write_dbl_analysis(s_MaxPathEdiff,&TResult::m_MaxPathEdiff);
    write_dbl_analysis(s_MaxUsedPathEdiff,&TResult::m_MaxUsedPathEdiff);
	write_uint32_analysis(s_MaxPathCount,&TResult::m_MaxPathCount);
	write_dbl_analysis(s_MeanPathCount,&TResult::m_MeanPathCount);
    write_dbl_analysis(s_MinUsedStateEnergy,&TResult::m_MinUsedStateEnergy);
    write_dbl_analysis(s_MaxUsedStateEnergy,&TResult::m_MaxUsedStateEnergy);
}

// Write analysis for results of multiple simulations
void MC::TResult::WriteMultiAnalysis(std::ostream& o_str, const std::vector<std::vector<std::unique_ptr<const TResult>>>& results)
{
    o_str << "Combined analysis: mean +/- stddev [min, max]" << std::endl;
    if (results.empty())
    {
        o_str << "Finished simulations: none" << std::endl;
        return;
    }

    // Search index of the first defined result object
    auto first_it = results.begin()->begin();
    std::size_t result_counter = 0;
    for (auto& vec : results)
    {
        if (result_counter == 0)
        {
            first_it = std::find_if(vec.begin(),vec.end(),[](const std::unique_ptr<const TResult>& ptr){return ptr.get() != nullptr;});
        }
        result_counter += std::count_if(vec.begin(),vec.end(),[](const std::unique_ptr<const TResult>& ptr){return ptr.get() != nullptr;});
    }
    o_str << "Finished simulations: " << result_counter << std::endl;
    if (result_counter <= 1)
    {
        return;
    }

    // Define lambda functions to calculate mean, stddev, min and max for a specified property
    auto write_dbl_analysis = [&] (const std::array<std::string,2>& desc,double TResult::* prop)
    {
        double min_val = first_it->get()->*prop;
        double max_val = first_it->get()->*prop;
        double stddev_val = 0.0;
        double mean_val = 0.0;
        double prev_mean_val = 0.0;
        std::size_t counter = 0;
        for (auto& vec : results)
        {
            for (auto it = vec.begin(); it != vec.end(); it++)
            {
                if (*it)
                {
                    if (it->get()->*prop < min_val) min_val = it->get()->*prop;
                    if (it->get()->*prop > max_val) max_val = it->get()->*prop;
                    counter++;
                    prev_mean_val = mean_val;
                    mean_val += (it->get()->*prop - mean_val)/static_cast<double>(counter);
                    stddev_val += (it->get()->*prop - mean_val) * (it->get()->*prop - prev_mean_val);
                }
            }
        }
        if ((counter > 1) && (stddev_val != 0.0))
        {
            stddev_val = sqrt(fabs(stddev_val) / static_cast<double>(counter - 1));
        } else stddev_val = 0.0;
        o_str << GF::CombineDescUnit(desc) << " = " << mean_val << " +/- " << stddev_val << " [" << min_val << ", " << max_val << "]" << std::endl;
    };
    auto write_uint32_analysis = [&] (const std::array<std::string,2>& desc, std::uint32_t TResult::* prop)
    {
        std::uint32_t min_val = first_it->get()->*prop;
        std::uint32_t max_val = first_it->get()->*prop;
        double mean_val = 0.0;
        double stddev_val = 0.0;
        double prev_mean_val = 0.0;
        std::size_t counter = 0;
        for (auto& vec : results)
        {
            for (auto it = vec.begin(); it != vec.end(); it++)
            {
                if (*it)
                {
                    if (it->get()->*prop < min_val) min_val = it->get()->*prop;
                    if (it->get()->*prop > max_val) max_val = it->get()->*prop;
                    counter++;
                    prev_mean_val = mean_val;
                    mean_val += (static_cast<double>(it->get()->*prop) - mean_val)/static_cast<double>(counter);
                    stddev_val += (static_cast<double>(it->get()->*prop) - mean_val) * 
                        (static_cast<double>(it->get()->*prop) - prev_mean_val);
                }
            }
        }
        if ((counter > 1) && (stddev_val != 0.0))
        {
            stddev_val = sqrt(fabs(stddev_val) / static_cast<double>(counter - 1));
        } else stddev_val = 0.0;
        o_str << GF::CombineDescUnit(desc) << " = " << mean_val << " +/- " << stddev_val << " [" << min_val << ", " << max_val << "]" << std::endl;
    };
    auto write_uint64_analysis = [&] (const std::array<std::string,2>& desc, std::uint64_t TResult::* prop)
    {
        std::uint64_t min_val = first_it->get()->*prop;
        std::uint64_t max_val = first_it->get()->*prop;
        double mean_val = 0.0;
        double stddev_val = 0.0;
        double prev_mean_val = 0.0;
        std::size_t counter = 0;
        for (auto& vec : results)
        {
            for (auto it = vec.begin(); it != vec.end(); it++)
            {
                if (*it)
                {
                    if (it->get()->*prop < min_val) min_val = it->get()->*prop;
                    if (it->get()->*prop > max_val) max_val = it->get()->*prop;
                    counter++;
                    prev_mean_val = mean_val;
                    mean_val += (static_cast<double>(it->get()->*prop) - mean_val)/static_cast<double>(counter);
                    stddev_val += (static_cast<double>(it->get()->*prop) - mean_val) * 
                        (static_cast<double>(it->get()->*prop) - prev_mean_val);
                }
            }
        }
        if ((counter > 1) && (stddev_val != 0.0))
        {
            stddev_val = sqrt(fabs(stddev_val) / static_cast<double>(counter - 1));
        } else stddev_val = 0.0;
        o_str << GF::CombineDescUnit(desc) << " = " << mean_val << " +/- " << stddev_val << " [" << min_val << ", " << max_val << "]" << std::endl;
    };

    // Analysis
    write_dbl_analysis(s_DriftConductivity,&TResult::m_DriftConductivity);
    write_dbl_analysis(s_DriftMobility,&TResult::m_DriftMobility);
    write_dbl_analysis(s_DiffusionCoefficient,&TResult::m_DiffusionCoefficient);
    write_dbl_analysis(s_DiffusionCoefficientParallel,&TResult::m_DiffusionCoefficientParallel);
    write_dbl_analysis(s_DiffusionCoefficientTransverse,&TResult::m_DiffusionCoefficientTransverse);
    write_dbl_analysis(s_HavenRatio,&TResult::m_HavenRatio);
    write_dbl_analysis(s_HavenRatioParallel,&TResult::m_HavenRatioParallel);
    write_dbl_analysis(s_HavenRatioTransverse,&TResult::m_HavenRatioTransverse);
    write_dbl_analysis(s_PartialEntropy,&TResult::m_PartialEntropy);
    write_dbl_analysis(s_EffChemPot,&TResult::m_EffChemPot);
    write_dbl_analysis(s_EffTemp,&TResult::m_EffTemp);
    write_dbl_analysis(s_EffCarriers,&TResult::m_EffCarriers);
    write_dbl_analysis(s_EffCarrierDensity,&TResult::m_EffCarrierDensity);
    write_uint32_analysis(s_ElectronCount,&TResult::m_ElectronCount);
    write_uint32_analysis(s_MobileElectrons,&TResult::m_MobileElectrons);
	write_uint32_analysis(s_ZeroHopElectrons,&TResult::m_ZeroHopElectrons);
    write_uint32_analysis(s_OscElectrons,&TResult::m_OscElectrons);
	write_dbl_analysis(s_TotalTime,&TResult::m_TotalTime);
    write_dbl_analysis(s_MeanDisp,&TResult::m_MeanDisp);
	write_dbl_analysis(s_MeanDisp_x,&TResult::m_MeanDisp_x);
	write_dbl_analysis(s_MeanDisp_y,&TResult::m_MeanDisp_y);
	write_dbl_analysis(s_MeanDisp_z,&TResult::m_MeanDisp_z);
    write_dbl_analysis(s_MeanDispVariance_x,&TResult::m_MeanDispVariance_x);
    write_dbl_analysis(s_MeanDispVariance_y,&TResult::m_MeanDispVariance_y);
    write_dbl_analysis(s_MeanDispVariance_z,&TResult::m_MeanDispVariance_z);
	write_dbl_analysis(s_MeanSquaredDisp_x,&TResult::m_MeanSquaredDisp_x);
	write_dbl_analysis(s_MeanSquaredDisp_y,&TResult::m_MeanSquaredDisp_y);
	write_dbl_analysis(s_MeanSquaredDisp_z,&TResult::m_MeanSquaredDisp_z);
    write_uint64_analysis(s_NonOscHops,&TResult::m_NonOscHops);
    write_dbl_analysis(s_NonOscHopRatio,&TResult::m_NonOscHopRatio);
    write_dbl_analysis(s_InXDirNonOscRatio,&TResult::m_InXDirNonOscRatio);
    write_dbl_analysis(s_MeanFieldContribution,&TResult::m_MeanFieldContribution);
    write_dbl_analysis(s_TotalEnergy,&TResult::m_TotalEnergy);
    write_dbl_analysis(s_CellSize,&TResult::m_CellSize);
    write_dbl_analysis(s_MaxPathDist,&TResult::m_MaxPathDist);
    write_dbl_analysis(s_MaxUsedPathDist,&TResult::m_MaxUsedPathDist);
    write_dbl_analysis(s_MaxPathEdiff,&TResult::m_MaxPathEdiff);
    write_dbl_analysis(s_MaxUsedPathEdiff,&TResult::m_MaxUsedPathEdiff);
	write_uint32_analysis(s_MaxPathCount,&TResult::m_MaxPathCount);
	write_dbl_analysis(s_MeanPathCount,&TResult::m_MeanPathCount);
    write_dbl_analysis(s_MinUsedStateEnergy,&TResult::m_MinUsedStateEnergy);
    write_dbl_analysis(s_MaxUsedStateEnergy,&TResult::m_MaxUsedStateEnergy);
}

// Write table header for mean and stddev of multiple simulations
std::vector<std::array<std::string,2>> MC::TResult::WriteMultiTableHeader()
{
    std::vector<std::array<std::string,2>> header;
    header.push_back(s_DriftConductivity);
    header.push_back(GF::StdDevDescUnit(s_DriftConductivity));
    header.push_back(s_DriftMobility);
    header.push_back(GF::StdDevDescUnit(s_DriftMobility));
    header.push_back(s_DiffusionCoefficient);
    header.push_back(GF::StdDevDescUnit(s_DiffusionCoefficient));
    header.push_back(s_DiffusionCoefficientParallel);
    header.push_back(GF::StdDevDescUnit(s_DiffusionCoefficientParallel));
    header.push_back(s_DiffusionCoefficientTransverse);
    header.push_back(GF::StdDevDescUnit(s_DiffusionCoefficientTransverse));
    header.push_back(s_HavenRatio);
    header.push_back(GF::StdDevDescUnit(s_HavenRatio));
    header.push_back(s_HavenRatioParallel);
    header.push_back(GF::StdDevDescUnit(s_HavenRatioParallel));
    header.push_back(s_HavenRatioTransverse);
    header.push_back(GF::StdDevDescUnit(s_HavenRatioTransverse));
    header.push_back(s_PartialEntropy);
    header.push_back(GF::StdDevDescUnit(s_PartialEntropy));
    header.push_back(s_EffChemPot);
    header.push_back(GF::StdDevDescUnit(s_EffChemPot));
    header.push_back(s_EffTemp);
    header.push_back(GF::StdDevDescUnit(s_EffTemp));
    header.push_back(s_EffCarriers);
    header.push_back(GF::StdDevDescUnit(s_EffCarriers));
    header.push_back(s_EffCarrierDensity);
    header.push_back(GF::StdDevDescUnit(s_EffCarrierDensity));
    header.push_back(s_ElectronCount);
    header.push_back(GF::StdDevDescUnit(s_ElectronCount));
    header.push_back(s_MobileElectrons);
    header.push_back(GF::StdDevDescUnit(s_MobileElectrons));
    header.push_back(s_ZeroHopElectrons);
    header.push_back(GF::StdDevDescUnit(s_ZeroHopElectrons));
    header.push_back(s_OscElectrons);
    header.push_back(GF::StdDevDescUnit(s_OscElectrons));
    header.push_back(s_TotalTime);
    header.push_back(GF::StdDevDescUnit(s_TotalTime));
    header.push_back(s_MeanDisp);
    header.push_back(GF::StdDevDescUnit(s_MeanDisp));
    header.push_back(s_MeanDisp_x);
    header.push_back(GF::StdDevDescUnit(s_MeanDisp_x));
    header.push_back(s_MeanDisp_y);
    header.push_back(GF::StdDevDescUnit(s_MeanDisp_y));
    header.push_back(s_MeanDisp_z);
    header.push_back(GF::StdDevDescUnit(s_MeanDisp_z));
    header.push_back(s_MeanDispVariance_x);
    header.push_back(GF::StdDevDescUnit(s_MeanDispVariance_x));
    header.push_back(s_MeanDispVariance_y);
    header.push_back(GF::StdDevDescUnit(s_MeanDispVariance_y));
    header.push_back(s_MeanDispVariance_z);
    header.push_back(GF::StdDevDescUnit(s_MeanDispVariance_z));
    header.push_back(s_MeanSquaredDisp_x);
    header.push_back(GF::StdDevDescUnit(s_MeanSquaredDisp_x));
    header.push_back(s_MeanSquaredDisp_y);
    header.push_back(GF::StdDevDescUnit(s_MeanSquaredDisp_y));
    header.push_back(s_MeanSquaredDisp_z);
    header.push_back(GF::StdDevDescUnit(s_MeanSquaredDisp_z));
    header.push_back(s_NonOscHops);
    header.push_back(GF::StdDevDescUnit(s_NonOscHops));
    header.push_back(s_NonOscHopRatio);
    header.push_back(GF::StdDevDescUnit(s_NonOscHopRatio));
    header.push_back(s_InXDirNonOscRatio);
    header.push_back(GF::StdDevDescUnit(s_InXDirNonOscRatio));
    header.push_back(s_MeanFieldContribution);
    header.push_back(GF::StdDevDescUnit(s_MeanFieldContribution));
    header.push_back(s_TotalEnergy);
    header.push_back(GF::StdDevDescUnit(s_TotalEnergy));
    header.push_back(s_CellSize);
    header.push_back(GF::StdDevDescUnit(s_CellSize));
    header.push_back(s_MaxPathDist);
    header.push_back(GF::StdDevDescUnit(s_MaxPathDist));
    header.push_back(GF::MaxDescUnit(s_MaxUsedPathDist));
    header.push_back(s_MaxPathEdiff);
    header.push_back(GF::StdDevDescUnit(s_MaxPathEdiff));
    header.push_back(GF::MaxDescUnit(s_MaxUsedPathEdiff));
    header.push_back(s_MaxPathCount);
    header.push_back(GF::StdDevDescUnit(s_MaxPathCount));
    header.push_back(s_MeanPathCount);
    header.push_back(GF::StdDevDescUnit(s_MeanPathCount));
    header.push_back(GF::MinDescUnit(s_MinUsedStateEnergy));
    header.push_back(GF::MaxDescUnit(s_MaxUsedStateEnergy));

    for (auto& item : header)
    {
        if (item[1].empty()) item[1] = "_";
    }

    return header;
}

// Write table line for mean and stddev of multiple simulations
std::vector<std::string> MC::TResult::WriteMultiTableLine(const std::vector<std::unique_ptr<const TResult>>& results)
{
    // Return empty vector if there are no defined results
    if (std::none_of(results.begin(),results.end(),[](const std::unique_ptr<const TResult>& ptr){return ptr.get() != nullptr;}))
    {
        return std::vector<std::string>();
    }

    std::vector<std::string> line;
    std::stringstream sstr;
    
    // Define lambda functions to calculate mean and stddev for a specified property
    auto write_dbl_analysis = [&] (double TResult::* prop)
    {
        double stddev_val = 0.0;
        double mean_val = 0.0;
        double prev_mean_val = 0.0;
        std::size_t counter = 0;
        for (auto it = results.begin(); it != results.end(); it++)
        {
            if (*it)
            {
                counter++;
                prev_mean_val = mean_val;
                mean_val += (it->get()->*prop - mean_val)/static_cast<double>(counter);
                stddev_val += (it->get()->*prop - mean_val) * (it->get()->*prop - prev_mean_val);
            }
        }
        if ((counter > 1) && (stddev_val != 0.0))
        {
            stddev_val = sqrt(fabs(stddev_val) / static_cast<double>(counter - 1));
        } else stddev_val = 0.0;

        sstr << mean_val;
        line.push_back(sstr.str());
        sstr.str("");

        sstr << stddev_val;
        line.push_back(sstr.str());
        sstr.str("");
    };
    auto write_uint32_analysis = [&] (std::uint32_t TResult::* prop)
    {
        double stddev_val = 0.0;
        double mean_val = 0.0;
        double prev_mean_val = 0.0;
        std::size_t counter = 0;
        for (auto it = results.begin(); it != results.end(); it++)
        {
            if (*it)
            {
                counter++;
                prev_mean_val = mean_val;
                mean_val += (static_cast<double>(it->get()->*prop) - mean_val)/static_cast<double>(counter);
                stddev_val += (static_cast<double>(it->get()->*prop) - mean_val) * 
                    (static_cast<double>(it->get()->*prop) - prev_mean_val);
            }
        }
        if ((counter > 1) && (stddev_val != 0.0))
        {
            stddev_val = sqrt(fabs(stddev_val) / static_cast<double>(counter - 1));
        } else stddev_val = 0.0;
        
        sstr << mean_val;
        line.push_back(sstr.str());
        sstr.str("");

        sstr << stddev_val;
        line.push_back(sstr.str());
        sstr.str("");
    };
    auto write_uint64_analysis = [&] (std::uint64_t TResult::* prop)
    {
        double stddev_val = 0.0;
        double mean_val = 0.0;
        double prev_mean_val = 0.0;
        std::size_t counter = 0;
        for (auto it = results.begin(); it != results.end(); it++)
        {
            if (*it)
            {
                counter++;
                prev_mean_val = mean_val;
                mean_val += (static_cast<double>(it->get()->*prop) - mean_val)/static_cast<double>(counter);
                stddev_val += (static_cast<double>(it->get()->*prop) - mean_val) * 
                    (static_cast<double>(it->get()->*prop) - prev_mean_val);
            }
        }
        if ((counter > 1) && (stddev_val != 0.0))
        {
            stddev_val = sqrt(fabs(stddev_val) / static_cast<double>(counter - 1));
        } else stddev_val = 0.0;
        
        sstr << mean_val;
        line.push_back(sstr.str());
        sstr.str("");

        sstr << stddev_val;
        line.push_back(sstr.str());
        sstr.str("");
    };
    auto write_dbl_min = [&] (double TResult::* prop)
    {
        double min_val = 0.0;
        bool has_data = false;
        for (auto it = results.begin(); it != results.end(); it++)
        {
            if (*it)
            {
                if (has_data)
                {
                    if (it->get()->*prop < min_val) min_val = it->get()->*prop;
                }
                else
                {
                    min_val = it->get()->*prop;
                    has_data = true;
                }
            }
        }

        sstr << min_val;
        line.push_back(sstr.str());
        sstr.str("");
    };
    auto write_dbl_max = [&] (double TResult::* prop)
    {
        double max_val = 0.0;
        bool has_data = false;
        for (auto it = results.begin(); it != results.end(); it++)
        {
            if (*it)
            {
                if (has_data)
                {
                    if (it->get()->*prop > max_val) max_val = it->get()->*prop;
                }
                else
                {
                    max_val = it->get()->*prop;
                    has_data = true;
                }
            }
        }

        sstr << max_val;
        line.push_back(sstr.str());
        sstr.str("");
    };

    // Analysis
    write_dbl_analysis(&TResult::m_DriftConductivity);
    write_dbl_analysis(&TResult::m_DriftMobility);
    write_dbl_analysis(&TResult::m_DiffusionCoefficient);
    write_dbl_analysis(&TResult::m_DiffusionCoefficientParallel);
    write_dbl_analysis(&TResult::m_DiffusionCoefficientTransverse);
    write_dbl_analysis(&TResult::m_HavenRatio);
    write_dbl_analysis(&TResult::m_HavenRatioParallel);
    write_dbl_analysis(&TResult::m_HavenRatioTransverse);
    write_dbl_analysis(&TResult::m_PartialEntropy);
    write_dbl_analysis(&TResult::m_EffChemPot);
    write_dbl_analysis(&TResult::m_EffTemp);
    write_dbl_analysis(&TResult::m_EffCarriers);
    write_dbl_analysis(&TResult::m_EffCarrierDensity);
    write_uint32_analysis(&TResult::m_ElectronCount);
    write_uint32_analysis(&TResult::m_MobileElectrons);
	write_uint32_analysis(&TResult::m_ZeroHopElectrons);
    write_uint32_analysis(&TResult::m_OscElectrons);
	write_dbl_analysis(&TResult::m_TotalTime);
    write_dbl_analysis(&TResult::m_MeanDisp);
	write_dbl_analysis(&TResult::m_MeanDisp_x);
	write_dbl_analysis(&TResult::m_MeanDisp_y);
	write_dbl_analysis(&TResult::m_MeanDisp_z);
    write_dbl_analysis(&TResult::m_MeanDispVariance_x);
    write_dbl_analysis(&TResult::m_MeanDispVariance_y);
    write_dbl_analysis(&TResult::m_MeanDispVariance_z);
	write_dbl_analysis(&TResult::m_MeanSquaredDisp_x);
	write_dbl_analysis(&TResult::m_MeanSquaredDisp_y);
	write_dbl_analysis(&TResult::m_MeanSquaredDisp_z);
    write_uint64_analysis(&TResult::m_NonOscHops);
    write_dbl_analysis(&TResult::m_NonOscHopRatio);
    write_dbl_analysis(&TResult::m_InXDirNonOscRatio);
    write_dbl_analysis(&TResult::m_MeanFieldContribution);
    write_dbl_analysis(&TResult::m_TotalEnergy);
    write_dbl_analysis(&TResult::m_CellSize);
    write_dbl_analysis(&TResult::m_MaxPathDist);
    write_dbl_max(&TResult::m_MaxUsedPathDist);
    write_dbl_analysis(&TResult::m_MaxPathEdiff);
    write_dbl_max(&TResult::m_MaxUsedPathEdiff);
	write_uint32_analysis(&TResult::m_MaxPathCount);
	write_dbl_analysis(&TResult::m_MeanPathCount);
    write_dbl_min(&TResult::m_MinUsedStateEnergy);
    write_dbl_max(&TResult::m_MaxUsedStateEnergy);

    return line;
}

// Save progress of equilibration or simulation
void MC::TResult::SaveProgress(double percentage, bool is_eq)
{
    std::vector<TProgress>::iterator new_progress;
    if (is_eq)
    {
        m_EqProgress.emplace_back();
        new_progress = m_EqProgress.end() - 1;
    }
    else
    {
        m_Progress.emplace_back();
        new_progress = m_Progress.end() - 1;
    }

    new_progress->m_Percentage = percentage;
    new_progress->m_TotalHops = m_TotalHops;
    new_progress->m_NonOscHops = m_NonOscHops;
    new_progress->m_NonOscHopRatio = m_NonOscHopRatio;
    new_progress->m_TotalTime = m_TotalTime;
	new_progress->m_DriftConductivity = m_DriftConductivity;
	new_progress->m_DriftMobility = m_DriftMobility;
    new_progress->m_DiffusionCoefficient = m_DiffusionCoefficient;
	new_progress->m_DiffusionCoefficientParallel = m_DiffusionCoefficientParallel;
	new_progress->m_DiffusionCoefficientTransverse = m_DiffusionCoefficientTransverse;
	new_progress->m_PartialEntropy = m_PartialEntropy;
    new_progress->m_TotalEnergy = m_TotalEnergy;
    new_progress->m_EffChemPot = m_EffChemPot;
    new_progress->m_EffTemp = m_EffTemp;
    new_progress->m_EffCarriers = m_EffCarriers;
	new_progress->m_EffCarrierDensity = m_EffCarrierDensity;
    new_progress->m_MobileElectrons = m_MobileElectrons;
	new_progress->m_ZeroHopElectrons = m_ZeroHopElectrons;
	new_progress->m_OscElectrons = m_OscElectrons;
	new_progress->m_MeanDisp = m_MeanDisp;
	new_progress->m_MeanDisp_x = m_MeanDisp_x;
	new_progress->m_MeanDisp_y = m_MeanDisp_y;
	new_progress->m_MeanDisp_z = m_MeanDisp_z;
	new_progress->m_MeanDispVariance_x = m_MeanDispVariance_x;
    new_progress->m_MeanDispVariance_y = m_MeanDispVariance_y;
    new_progress->m_MeanDispVariance_z = m_MeanDispVariance_z;
	new_progress->m_MeanSquaredDisp_x = m_MeanSquaredDisp_x;
	new_progress->m_MeanSquaredDisp_y = m_MeanSquaredDisp_y;
	new_progress->m_MeanSquaredDisp_z = m_MeanSquaredDisp_z;
    new_progress->m_MaxUsedPathDist = m_MaxUsedPathDist;
    new_progress->m_MaxUsedPathEdiff = m_MaxUsedPathEdiff;
    new_progress->m_MinUsedStateEnergy = m_MinUsedStateEnergy;
    new_progress->m_MaxUsedStateEnergy = m_MaxUsedStateEnergy;
}

// Write progress table header
void MC::TResult::WriteProgressHeader(std::ostream& o_str, std::uint64_t maxhops, bool has_field) const
{
    o_str << std::setfill(' ') << std::left;
    int pos_dbl_width = static_cast<int>(o_str.precision()) + 5;
    int hops_width = static_cast<int>(log10(maxhops)) + 1;
    int el_width = static_cast<int>(log10(m_ElectronCount)) + 1;

    // Write header
    o_str << "Progress %: ";
    o_str << std::setw(std::max(hops_width,4)) << "Hops" << " ";
    o_str << std::setw(std::max(pos_dbl_width,9)) << "Hops/Neff" << " ";
    o_str << std::setw(std::max(hops_width,10)) << "NonOscHops" << " ";
    o_str << std::setw(std::max(pos_dbl_width,15)) << "NonOscHops/Neff" << " ";
    o_str << std::setw(std::max(pos_dbl_width,11)) << "Timespan(s)" << " -> ";
    
    if (has_field)
    {
        o_str << std::setw(std::max(pos_dbl_width,10)) << "Sigma(S/m)" << " ";
        o_str << std::setw(std::max(pos_dbl_width,9)) << "Dx(cm2/s)" << " ";
        o_str << std::setw(std::max(pos_dbl_width,10)) << "Dyz(cm2/s)" << " ";
    }
    else
    {
        o_str << std::setw(std::max(pos_dbl_width,8)) << "D(cm2/s)" << " ";
    }

    o_str << std::setw(std::max(pos_dbl_width + 1,10)) << "Etotal(eV)" << " ";
    o_str << std::setw(std::max(pos_dbl_width,4)) << "Neff" << " ";
    o_str << std::setw(std::max(el_width,4)) << "Nmob" << " ";
    o_str << std::setw(std::max(el_width,5)) << "Nzero" << " ";
    o_str << std::setw(std::max(el_width,4)) << "Nosc" << " ";

    if (has_field)
    {
        o_str << std::setw(std::max(pos_dbl_width + 1,7)) << "<x>(nm)" << " ";
        o_str << std::setw(std::max(pos_dbl_width,19)) << "<(x-<x>mob)^2>(nm2)" << " ";
        o_str << std::setw(std::max(pos_dbl_width,20)) << "(<y^2>+<z^2>)/2(nm2)";
    }
    else
    {
        o_str << std::setw(std::max(pos_dbl_width,26)) << "(<x^2>+<y^2>+<z^2>)/3(nm2)";
    }
    o_str << std::endl;
}

// Write last progress line of equilibration or simulation
void MC::TResult::WriteProgressLine(std::ostream& o_str, std::uint64_t maxhops, bool has_field, bool is_eq) const
{
    o_str << std::setfill(' ') << std::left;
    int pos_dbl_width = static_cast<int>(o_str.precision()) + 5;
    int hops_width = static_cast<int>(log10(maxhops)) + 1;
    int el_width = static_cast<int>(log10(m_ElectronCount)) + 1;

    std::vector<TProgress>::const_iterator last;
    if (is_eq)
    {
        if (m_EqProgress.empty()) return;
        last = m_EqProgress.cend() - 1;
    }
    else
    {
        if (m_Progress.empty()) return;
        last = m_Progress.cend() - 1;
    }
    
    // Write progress line
    o_str << std::setw(8) << last->m_Percentage << " %: ";
    o_str << std::setw(std::max(hops_width,4)) << last->m_TotalHops << " ";
    o_str << std::setw(std::max(pos_dbl_width,9)) 
        << static_cast<double>(last->m_TotalHops)/last->m_EffCarriers << " ";
    o_str << std::setw(std::max(hops_width,10)) << last->m_NonOscHops << " ";
    o_str << std::setw(std::max(pos_dbl_width,15)) 
        << static_cast<double>(last->m_NonOscHops)/last->m_EffCarriers << " ";
    o_str << std::setw(std::max(pos_dbl_width,11)) << last->m_TotalTime << " -> ";
    
    if (has_field)
    {
        o_str << std::setw(std::max(pos_dbl_width,10)) << last->m_DriftConductivity << " ";
        o_str << std::setw(std::max(pos_dbl_width,9)) << last->m_DiffusionCoefficientParallel << " ";
        o_str << std::setw(std::max(pos_dbl_width,10)) << last->m_DiffusionCoefficientTransverse << " ";
    }
    else
    {
        o_str << std::setw(std::max(pos_dbl_width,8)) << last->m_DiffusionCoefficient << " ";
    }

    o_str << std::setw(std::max(pos_dbl_width + 1,10)) << last->m_TotalEnergy << " ";
    o_str << std::setw(std::max(pos_dbl_width,4)) << last->m_EffCarriers << " ";
    o_str << std::setw(std::max(el_width,4)) << last->m_MobileElectrons << " ";
    o_str << std::setw(std::max(el_width,5)) << last->m_ZeroHopElectrons << " ";
    o_str << std::setw(std::max(el_width,4)) << last->m_OscElectrons << " ";

    if (has_field)
    {
        o_str << std::setw(std::max(pos_dbl_width + 1,7)) << last->m_MeanDisp_x << " ";
        o_str << std::setw(std::max(pos_dbl_width,19)) << last->m_MeanDispVariance_x << " ";
        o_str << std::setw(std::max(pos_dbl_width,20)) 
            << (last->m_MeanSquaredDisp_y + last->m_MeanSquaredDisp_z)/2.0;
    }
    else
    {
        o_str << std::setw(std::max(pos_dbl_width,26)) 
            << (last->m_MeanSquaredDisp_x + last->m_MeanSquaredDisp_y + last->m_MeanSquaredDisp_z)/3.0;
    }
    o_str << std::endl;
}

// Write convergence table of equilibration or simulation
void MC::TResult::WriteConvergenceTable(std::ostream& o_str, bool is_eq) const
{
    if (is_eq)
    {
        if (m_EqProgress.empty()) return;
    }
    else
    {
        if (m_Progress.empty()) return;
    }

    // Create header
    std::vector<std::array<std::string,2>> header;
    header.push_back(std::array<std::string,2>{"Progress","%"});
    header.push_back(std::array<std::string,2>{"Hops",""});
    header.push_back(std::array<std::string,2>{"MeanHops",""});
    header.push_back(std::array<std::string,2>{"OscHops",""});
    header.push_back(std::array<std::string,2>{"MeanOscHops",""});
    header.push_back(s_NonOscHops);
    header.push_back(std::array<std::string,2>{"MeanNonOscHops",""});
    header.push_back(s_NonOscHopRatio);
    header.push_back(s_TotalTime);
    header.push_back(s_DriftConductivity);
    header.push_back(s_DriftMobility);
    header.push_back(s_DiffusionCoefficient);
    header.push_back(s_DiffusionCoefficientParallel);
    header.push_back(s_DiffusionCoefficientTransverse);
    header.push_back(s_PartialEntropy);
    header.push_back(s_TotalEnergy);
    header.push_back(s_EffChemPot);
    header.push_back(s_EffTemp);
    header.push_back(s_EffCarriers);
    header.push_back(s_EffCarrierDensity);
    header.push_back(s_MobileElectrons);
    header.push_back(s_ZeroHopElectrons);
    header.push_back(s_OscElectrons);
    header.push_back(s_MeanDisp);
    header.push_back(s_MeanDisp_x);
    header.push_back(s_MeanDisp_y);
    header.push_back(s_MeanDisp_z);
    header.push_back(s_MeanDispVariance_x);
    header.push_back(s_MeanDispVariance_y);
    header.push_back(s_MeanDispVariance_z);
    header.push_back(s_MeanSquaredDisp_x);
    header.push_back(s_MeanSquaredDisp_y);
    header.push_back(s_MeanSquaredDisp_z);
    header.push_back(std::array<std::string,2>{"MeanSqDispyz","nm2"});
    header.push_back(std::array<std::string,2>{"MeanSqDispxyz","nm2"});
    header.push_back(std::array<std::string,2>{"MeanDispVarxSqyz","nm2"});
    header.push_back(std::array<std::string,2>{"MeanDispVaryz","nm2"});
    header.push_back(std::array<std::string,2>{"MeanDispVarxyz","nm2"});
    header.push_back(s_MaxUsedPathDist);
    header.push_back(s_MaxUsedPathEdiff);
    header.push_back(s_MinUsedStateEnergy);
    header.push_back(s_MaxUsedStateEnergy);

    for (auto& item : header)
    {
        if (item[1].empty()) item[1] = "_";
    }

    // Create table
    std::vector<std::vector<std::string>> table;
    std::vector<TProgress>::const_iterator cbegin = m_Progress.cbegin();
    std::vector<TProgress>::const_iterator cend = m_Progress.cend();
    if (is_eq)
    {
        cbegin = m_EqProgress.cbegin();
        cend = m_EqProgress.cend();
    }
    std::stringstream sstr;
    for (auto it = cbegin; it != cend; it++)
    {
        table.emplace_back();

        double density_factor = it->m_EffCarriers/m_EffCarriers;

        sstr << it->m_Percentage;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_TotalHops;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << static_cast<double>(it->m_TotalHops)/it->m_EffCarriers;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_TotalHops - it->m_NonOscHops;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << static_cast<double>(it->m_TotalHops - it->m_NonOscHops)/it->m_EffCarriers;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_NonOscHops;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << static_cast<double>(it->m_NonOscHops)/it->m_EffCarriers;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_NonOscHopRatio;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_TotalTime;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_DriftConductivity;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_DriftMobility * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_DiffusionCoefficient * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_DiffusionCoefficientParallel * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_DiffusionCoefficientTransverse * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_PartialEntropy;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_TotalEnergy;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_EffChemPot;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_EffTemp;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_EffCarriers;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_EffCarrierDensity;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MobileElectrons;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_ZeroHopElectrons;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_OscElectrons;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MeanDisp * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MeanDisp_x * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MeanDisp_y * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MeanDisp_z * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MeanDispVariance_x * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MeanDispVariance_y * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MeanDispVariance_z * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MeanSquaredDisp_x * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MeanSquaredDisp_y * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MeanSquaredDisp_z * density_factor;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << (it->m_MeanSquaredDisp_y + it->m_MeanSquaredDisp_z) * density_factor/2.0;
        table.back().push_back(sstr.str());
        sstr.str("");
 
        sstr << (it->m_MeanSquaredDisp_x + it->m_MeanSquaredDisp_y + it->m_MeanSquaredDisp_z) * density_factor/3.0;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << (it->m_MeanDispVariance_x + it->m_MeanSquaredDisp_y + it->m_MeanSquaredDisp_z) * density_factor/3.0;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << (it->m_MeanDispVariance_y + it->m_MeanDispVariance_z) * density_factor/2.0;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << (it->m_MeanDispVariance_x + it->m_MeanDispVariance_y + it->m_MeanDispVariance_z) * density_factor/3.0;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MaxUsedPathDist;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MaxUsedPathEdiff;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MinUsedStateEnergy;
        table.back().push_back(sstr.str());
        sstr.str("");

        sstr << it->m_MaxUsedStateEnergy;
        table.back().push_back(sstr.str());
        sstr.str("");
    }

    GF::WriteTable(o_str, header, table);
}