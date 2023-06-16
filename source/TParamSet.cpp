#include "TParamSet.hpp"

#include <sstream>
#include <regex>
#include <iostream>
#include <algorithm>

#include "Constants.hpp"
#include "GlobalFunctions.hpp"
#include "CustomExceptions.hpp"

// Parameter descriptors and units (without white-spaces)
const std::array<std::string,2> MC::TParamSet::s_SimID = {"SimID",""};
const std::array<std::string,2> MC::TParamSet::s_Repetitions = {"Repetitions",""};
const std::array<std::string,2> MC::TParamSet::s_MinStateEnergy = {"Emin","eV"};
const std::array<std::string,2> MC::TParamSet::s_MaxStateEnergy = {"Emax","eV"};
const std::array<std::string,2> MC::TParamSet::s_ChemPot = {"ChemPot","eV"};
const std::array<std::string,2> MC::TParamSet::s_StateCount = {"States",""};
const std::array<std::string,2> MC::TParamSet::s_MinPathCount = {"MinPaths",""};
const std::array<std::string,2> MC::TParamSet::s_DistCutoff = {"DistCutoff","nm"};
const std::array<std::string,2> MC::TParamSet::s_EdiffCutoff = {"EdiffCutoff","eV"};
const std::array<std::string,2> MC::TParamSet::s_PreHopLimit = {"PreHops",""};
const std::array<std::string,2> MC::TParamSet::s_EqHopLimit = {"EqHops",""};
const std::array<std::string,2> MC::TParamSet::s_HopLimit = {"HopLimit",""};
const std::array<std::string,2> MC::TParamSet::s_RndSeed = {"Seed",""};
const std::array<std::string,2> MC::TParamSet::s_AttemptTime = {"AttTime","s"};
const std::array<std::string,2> MC::TParamSet::s_Temperature = {"Temp","K"};
const std::array<std::string,2> MC::TParamSet::s_PhiGradient = {"GradPhi","V/cm"};
const std::array<std::string,2> MC::TParamSet::s_LocRadius = {"LocRadius","nm"};

// Default constructor
MC::TParamSet::TParamSet()
    : m_SimID(0),
    m_EFTAdjust(true), m_InitialFDDistrib(true), 
    m_TeffFit(true), m_EnforceECount(true), m_CutoffAutoAdjust(false),
    m_DistCutoffAdjustPercentage(0.0), m_EdiffCutoffAdjustPercentage(0.0),
    m_UseYZVariance(false),
    m_Repetitions(1), c_Repetitions(true),
    m_MinStateEnergy(0.0), c_MinStateEnergy(true),
    m_MaxStateEnergy(0.0), c_MaxStateEnergy(true),
	m_ChemPot(0.0), c_ChemPot(true),
	mI_StateCount(0), c_StateCount(true),
	m_MinPathCount(0), c_MinPathCount(true),
    m_DistCutoff(0.0), c_DistCutoff(true),
    m_EdiffCutoff(0.0), c_EdiffCutoff(true),
    mI_PreHopLimit(0), c_PreHopLimit(true),
    mI_EqHopLimit(0), c_EqHopLimit(true),
	mI_HopLimit(0), c_HopLimit(true),
    m_RndSeed(0), c_RndSeed(true),
	m_AttemptTime(0.0), c_AttemptTime(true),
	m_Temperature(0.0), c_Temperature(true),
	m_PhiGradient(0.0), c_PhiGradient(true),
	m_LocRadius(0.0), c_LocRadius(true)
{

}

// Equality comparison
bool MC::TParamSet::operator==(const TParamSet& rhs) const
{
    return ((m_SimID == rhs.m_SimID) &&
        (GF::AlmostEqual(m_MinStateEnergy, rhs.m_MinStateEnergy)) &&
        (GF::AlmostEqual(m_MaxStateEnergy, rhs.m_MaxStateEnergy)) &&
        (GF::AlmostEqual(m_ChemPot, rhs.m_ChemPot)) &&
        (mI_StateCount == rhs.mI_StateCount) &&
        (m_MinPathCount == rhs.m_MinPathCount) &&
        (GF::AlmostEqual(m_DistCutoff, rhs.m_DistCutoff)) &&
        (GF::AlmostEqual(m_EdiffCutoff, rhs.m_EdiffCutoff)) &&
        (mI_PreHopLimit == rhs.mI_PreHopLimit) &&
        (mI_EqHopLimit == rhs.mI_EqHopLimit) &&
        (mI_HopLimit == rhs.mI_HopLimit) &&
        (GF::AlmostEqual(m_AttemptTime, rhs.m_AttemptTime)) &&
        (GF::AlmostEqual(m_Temperature, rhs.m_Temperature)) &&
        (GF::AlmostEqual(m_PhiGradient, rhs.m_PhiGradient)) &&
        (GF::AlmostEqual(m_LocRadius, rhs.m_LocRadius)));
}

// Write all parameters
std::string MC::TParamSet::Write(bool is_output) const
{
    std::stringstream sstr;
    if (is_output) sstr << GF::CombineDescUnit(s_SimID) << " = " << m_SimID << std::endl;
    sstr << GF::CombineDescUnit(s_Repetitions) << " = " << m_Repetitions << std::endl;
    sstr << GF::CombineDescUnit(s_MinStateEnergy) << " = " << m_MinStateEnergy << std::endl;
    sstr << GF::CombineDescUnit(s_MaxStateEnergy) << " = " << m_MaxStateEnergy << std::endl;
    sstr << GF::CombineDescUnit(s_ChemPot) << " = " << m_ChemPot << std::endl;
    sstr << GF::CombineDescUnit(s_StateCount) << " = " << mO_StateCount() << std::endl;
    sstr << GF::CombineDescUnit(s_MinPathCount) << " = " << m_MinPathCount << std::endl;
    sstr << GF::CombineDescUnit(s_DistCutoff) << " = " << m_DistCutoff << std::endl;
    sstr << GF::CombineDescUnit(s_EdiffCutoff) << " = " << m_EdiffCutoff << std::endl;
    sstr << GF::CombineDescUnit(s_PreHopLimit) << " = " << mO_PreHopLimit() << std::endl;
    sstr << GF::CombineDescUnit(s_EqHopLimit) << " = " << mO_EqHopLimit() << std::endl;
    sstr << GF::CombineDescUnit(s_HopLimit) << " = " << mO_HopLimit() << std::endl;
    sstr << GF::CombineDescUnit(s_RndSeed) << " = " << m_RndSeed << std::endl;
    sstr << GF::CombineDescUnit(s_AttemptTime) << " = " << m_AttemptTime << std::endl;
    sstr << GF::CombineDescUnit(s_Temperature) << " = " << m_Temperature << std::endl;
    sstr << GF::CombineDescUnit(s_PhiGradient) << " = " << m_PhiGradient << std::endl;
    sstr << GF::CombineDescUnit(s_LocRadius) << " = " << m_LocRadius << std::endl;
    return sstr.str();
}

// Write non-varied parameters
std::string MC::TParamSet::WriteConstant() const
{
    std::stringstream sstr;
    if (c_Repetitions) sstr << GF::CombineDescUnit(s_Repetitions) << " = " << m_Repetitions << std::endl;
    if (c_MinStateEnergy) sstr << GF::CombineDescUnit(s_MinStateEnergy) << " = " << m_MinStateEnergy << std::endl;
    if (c_MaxStateEnergy) sstr << GF::CombineDescUnit(s_MaxStateEnergy) << " = " << m_MaxStateEnergy << std::endl;
    if (c_ChemPot) sstr << GF::CombineDescUnit(s_ChemPot) << " = " << m_ChemPot << std::endl;
    if (c_StateCount) sstr << GF::CombineDescUnit(s_StateCount) << " = " << mO_StateCount() << std::endl;
    if (c_MinPathCount) sstr << GF::CombineDescUnit(s_MinPathCount) << " = " << m_MinPathCount << std::endl;
    if (c_DistCutoff) sstr << GF::CombineDescUnit(s_DistCutoff) << " = " << m_DistCutoff << std::endl;
    if (c_EdiffCutoff) sstr << GF::CombineDescUnit(s_EdiffCutoff) << " = " << m_EdiffCutoff << std::endl;
    if (c_PreHopLimit) sstr << GF::CombineDescUnit(s_PreHopLimit) << " = " << mO_PreHopLimit() << std::endl;
    if (c_EqHopLimit) sstr << GF::CombineDescUnit(s_EqHopLimit) << " = " << mO_EqHopLimit() << std::endl;
    if (c_HopLimit) sstr << GF::CombineDescUnit(s_HopLimit) << " = " << mO_HopLimit() << std::endl;
    if (c_RndSeed) sstr << GF::CombineDescUnit(s_RndSeed) << " = " << m_RndSeed << std::endl;
    if (c_AttemptTime) sstr << GF::CombineDescUnit(s_AttemptTime) << " = " << m_AttemptTime << std::endl;
    if (c_Temperature) sstr << GF::CombineDescUnit(s_Temperature) << " = " << m_Temperature << std::endl;
    if (c_PhiGradient) sstr << GF::CombineDescUnit(s_PhiGradient) << " = " << m_PhiGradient << std::endl;
    if (c_LocRadius) sstr << GF::CombineDescUnit(s_LocRadius) << " = " << m_LocRadius << std::endl;
    return sstr.str();
}

// Write table header (description and unit of varied parameters as individual elements)
std::vector<std::array<std::string,2>> MC::TParamSet::WriteTableHeader(bool is_output) const
{
    std::vector<std::array<std::string,2>> header;
    if (is_output) header.push_back(s_SimID);
    if (!c_Repetitions) header.push_back(s_Repetitions);
    if (!c_MinStateEnergy) header.push_back(s_MinStateEnergy);
    if (!c_MaxStateEnergy) header.push_back(s_MaxStateEnergy);
    if (!c_ChemPot) header.push_back(s_ChemPot);
    if (!c_StateCount) header.push_back(s_StateCount);
    if (!c_MinPathCount) header.push_back(s_MinPathCount);
    if (!c_DistCutoff) header.push_back(s_DistCutoff);
    if (!c_EdiffCutoff) header.push_back(s_EdiffCutoff);
    if (!c_PreHopLimit) header.push_back(s_PreHopLimit);
    if (!c_EqHopLimit) header.push_back(s_EqHopLimit);
    if (!c_HopLimit) header.push_back(s_HopLimit);
    if (!c_RndSeed) header.push_back(s_RndSeed);
    if (!c_AttemptTime) header.push_back(s_AttemptTime);
    if (!c_Temperature) header.push_back(s_Temperature);
    if (!c_PhiGradient) header.push_back(s_PhiGradient);
    if (!c_LocRadius) header.push_back(s_LocRadius);

    for (auto& item : header)
    {
        if (item[1].empty()) item[1] = "_";
    }

    return header;
}

// Write table line (values of varied parameters as individual elements)
std::vector<std::string> MC::TParamSet::WriteTableLine(bool is_output) const
{
    std::vector<std::string> line;
    std::stringstream sstr;
    if (is_output)
    {
        sstr << m_SimID;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_Repetitions)
    {
        sstr << m_Repetitions;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_MinStateEnergy)
    {
        sstr << m_MinStateEnergy;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_MaxStateEnergy)
    {
        sstr << m_MaxStateEnergy;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_ChemPot)
    {
        sstr << m_ChemPot;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_StateCount)
    {
        sstr << mO_StateCount();
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_MinPathCount)
    {
        sstr << m_MinPathCount;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_DistCutoff)
    {
        sstr << m_DistCutoff;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_EdiffCutoff)
    {
        sstr << m_EdiffCutoff;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_PreHopLimit)
    {
        sstr << mO_PreHopLimit();
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_EqHopLimit)
    {
        sstr << mO_EqHopLimit();
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_HopLimit)
    {
        sstr << mO_HopLimit();
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_RndSeed)
    {
        sstr << m_RndSeed;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_AttemptTime)
    {
        sstr << m_AttemptTime;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_Temperature)
    {
        sstr << m_Temperature;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_PhiGradient)
    {
        sstr << m_PhiGradient;
        line.push_back(sstr.str());
        sstr.str("");
    }
    if (!c_LocRadius)
    {
        sstr << m_LocRadius;
        line.push_back(sstr.str());
        sstr.str("");
    }
    return line;
}

// Read parameters (true = successful read)
bool MC::TParamSet::Read(const std::string& general, const std::string& header, const std::string& line, bool raise_errors)
{
    // Generate lists of header and line items
    auto get_list = [] (const std::string& str) -> std::vector<std::string>
    {
        std::vector<std::string> items;
        if (!str.empty())
        {
            std::smatch match;
            std::string::const_iterator start(str.cbegin());
            while (std::regex_search(start,str.cend(),match,std::regex("\\S+")))
            {
                items.push_back(match[0]);  
                start = match.suffix().first;
            }
        }
        return items;
    };
    std::vector<std::string> header_items = get_list(header);
    std::vector<std::string> line_items = get_list(line);
    if (header_items.size() != line_items.size())
    {
        if (raise_errors) throw EX::TInvalidInput("Unequal number of header items and table entries.");
        return false;
    }

    // Define lambda functions to get parameters
    bool is_essential = false;
    auto get_dbl = [&] (const std::array<std::string,2>& desc, double& val, bool& cflag) -> bool
    {
        cflag = true;
        auto head_it = std::find(header_items.cbegin(),header_items.cend(),desc[0]);
        if (head_it == header_items.cend())
        {
            head_it = std::find(header_items.cbegin(),header_items.cend(),GF::CombineDescUnit(desc));
        }
        if (head_it != header_items.cend())
        {
            cflag = false;
            auto line_it = line_items.cbegin() + std::distance(header_items.cbegin(),head_it);
            if (!((std::regex_match(*line_it,std::regex(Constant::dblex))) &&
                (GF::StringToDouble(*line_it,val))))
            {
                if ((is_essential) && (raise_errors)) 
                    throw EX::TInvalidInput(desc[0] + " value \"" + *line_it + "\" is not a real number.");
                return false;
            }
        } else {
            std::smatch match;
            if (!((std::regex_search(general, match, 
                std::regex(GF::DescRegex(desc) + "\\s*=\\s*(" + Constant::dblex + ")"))) &&
                (GF::StringToDouble(match[1],val))))
            {
                if ((is_essential) && (raise_errors)) 
                    throw EX::TInvalidInput("Parameter \"" + desc[0] + "\" is missing.");
                return false;
            }
        }
        return true;
    };
    auto get_uint32 = [&] (const std::array<std::string,2>& desc, std::uint32_t& val, bool& cflag) -> bool
    {
        cflag = true;
        auto head_it = std::find(header_items.cbegin(),header_items.cend(),desc[0]);
        if (head_it == header_items.cend())
        {
            head_it = std::find(header_items.cbegin(),header_items.cend(),GF::CombineDescUnit(desc));
        }
        if (head_it != header_items.cend())
        {
            cflag = false;
            auto line_it = line_items.cbegin() + std::distance(header_items.cbegin(),head_it);
            if (!((std::regex_match(*line_it,std::regex(Constant::uintex))) &&
                (GF::StringToUInt32(*line_it,val))))
            {
                if ((is_essential) && (raise_errors)) 
                    throw EX::TInvalidInput(desc[0] + " value \"" + *line_it + "\" is not an unsigned integer.");
                return false;
            }
        } else {
            std::smatch match;
            if (!((std::regex_search(general, match, 
                std::regex(GF::DescRegex(desc) + "\\s*=\\s*(" + Constant::uintex + ")"))) &&
                (GF::StringToUInt32(match[1],val))))
            {
                if ((is_essential) && (raise_errors)) 
                    throw EX::TInvalidInput("Parameter \"" + desc[0] + "\" is missing.");
                return false;
            }
        }
        return true;
    };
    auto get_uint64 = [&] (const std::array<std::string,2>& desc, std::uint64_t& val, bool& cflag) -> bool
    {
        cflag = true;
        auto head_it = std::find(header_items.cbegin(),header_items.cend(),desc[0]);
        if (head_it == header_items.cend())
        {
            head_it = std::find(header_items.cbegin(),header_items.cend(),GF::CombineDescUnit(desc));
        }
        if (head_it != header_items.cend())
        {
            cflag = false;
            auto line_it = line_items.cbegin() + std::distance(header_items.cbegin(),head_it);
            if (!((std::regex_match(*line_it,std::regex(Constant::uintex))) &&
                (GF::StringToUInt64(*line_it,val))))
            {
                if ((is_essential) && (raise_errors)) 
                    throw EX::TInvalidInput(desc[0] + " value \"" + *line_it + "\" is not an unsigned integer.");
                return false;
            }
        } else {
            std::smatch match;
            if (!((std::regex_search(general, match, 
                std::regex(GF::DescRegex(desc) + "\\s*=\\s*(" + Constant::uintex + ")"))) &&
                (GF::StringToUInt64(match[1],val))))
            {
                if ((is_essential) && (raise_errors)) 
                    throw EX::TInvalidInput("Parameter \"" + desc[0] + "\" is missing.");
                return false;
            }
        }
        return true;
    };
    auto get_int64 = [&] (const std::array<std::string,2>& desc, std::int64_t& val, bool& cflag) -> bool
    {
        cflag = true;
        auto head_it = std::find(header_items.cbegin(),header_items.cend(),desc[0]);
        if (head_it == header_items.cend())
        {
            head_it = std::find(header_items.cbegin(),header_items.cend(),GF::CombineDescUnit(desc));
        }
        if (head_it != header_items.cend())
        {
            cflag = false;
            auto line_it = line_items.cbegin() + std::distance(header_items.cbegin(),head_it);
            if (!((std::regex_match(*line_it,std::regex(Constant::intex))) &&
                (GF::StringToInt64(*line_it,val))))
            {
                if ((is_essential) && (raise_errors)) 
                    throw EX::TInvalidInput(desc[0] + " value \"" + *line_it + "\" is not an integer.");
                return false;
            }
        } else {
            std::smatch match;
            if (!((std::regex_search(general, match, 
                std::regex(GF::DescRegex(desc) + "\\s*=\\s*(" + Constant::intex + ")"))) &&
                (GF::StringToInt64(match[1],val))))
            {
                if ((is_essential) && (raise_errors)) 
                    throw EX::TInvalidInput("Parameter \"" + desc[0] + "\" is missing.");
                return false;
            }
        }
        return true;
    };
    
    // Retrieve parameters from table line (first priority) or general parameters (second priority)
    is_essential = false;
    bool c_SimID = false;
    if (get_uint32(s_SimID,m_SimID,c_SimID) == false) m_SimID = 0;
    if (get_uint32(s_Repetitions,m_Repetitions,c_Repetitions) == false) m_Repetitions = 1;
    
    is_essential = true;
    if (get_dbl(s_MinStateEnergy,m_MinStateEnergy,c_MinStateEnergy) == false) return false;
    if (get_dbl(s_MaxStateEnergy,m_MaxStateEnergy,c_MaxStateEnergy) == false) return false;
    if (get_dbl(s_ChemPot,m_ChemPot,c_ChemPot) == false) return false;
    if (get_uint32(s_StateCount,mI_StateCount,c_StateCount) == false) return false;

    is_essential = false;
    if (get_uint32(s_MinPathCount,m_MinPathCount,c_MinPathCount) == false) m_MinPathCount = 0;
    if (get_dbl(s_DistCutoff,m_DistCutoff,c_DistCutoff) == false) m_DistCutoff = 0.0;
    if ((m_MinPathCount == 0) && (m_DistCutoff == 0.0))
    {
        if (raise_errors) throw EX::TInvalidInput("Parameters \"" + s_MinPathCount[0] + "\" and \"" + s_DistCutoff[0] + "\" are both missing or zero.");
        return false;
    }
    if (get_dbl(s_EdiffCutoff,m_EdiffCutoff,c_EdiffCutoff) == false) m_EdiffCutoff = 0.0;
    if (get_uint64(s_PreHopLimit,mI_PreHopLimit,c_PreHopLimit) == false) mI_PreHopLimit = 0;
    if (get_uint64(s_EqHopLimit,mI_EqHopLimit,c_EqHopLimit) == false) mI_EqHopLimit = 0;

    is_essential = true;
    if (get_uint64(s_HopLimit,mI_HopLimit,c_HopLimit) == false) return false;
    if (get_int64(s_RndSeed,m_RndSeed,c_RndSeed) == false) return false;
    if (get_dbl(s_AttemptTime,m_AttemptTime,c_AttemptTime) == false) return false;
    if (get_dbl(s_Temperature,m_Temperature,c_Temperature) == false) return false;
    if (get_dbl(s_PhiGradient,m_PhiGradient,c_PhiGradient) == false) return false;
    if (get_dbl(s_LocRadius,m_LocRadius,c_LocRadius) == false) return false;

    // Convert states and hop counts to internal values (for just one spin-type)
    // (rounded down to avoid potential overflow in output values)
    mI_StateCount /= 2U;
    mI_PreHopLimit /= 2U;
    mI_EqHopLimit /= 2U;
    mI_HopLimit /= 2U;

    return true;
}