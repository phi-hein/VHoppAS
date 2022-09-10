#include "THistogram.hpp"

#include <stdexcept>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>

#include "GlobalFunctions.hpp"
#include "CustomExceptions.hpp"

// Default constructor
MC::THistogram::THistogram()
    : m_IsConfigured(false), m_IsLog(false), m_xLabel(""), m_xUnit(""), m_xMinimum(0.0), m_xMaximum(0.0)
{

}

// Configure x-axis with linear scale
void MC::THistogram::ConfigureLinear(const std::string& xlabel, const std::string& xunit,
    double xmin, double xmax, std::size_t bincount)
{
    if ((xmin >= xmax) || (bincount == 0))
    {
        throw EX::TInvalidStatus("Invalid histogram specifications.",__func__);
    }

    m_xLabel = xlabel;
    m_xUnit = xunit;
    m_IsLog = false;
    m_xMinimum = xmin;
    m_xMaximum = xmax;
    try
    {
        m_xValues = std::vector<double>(bincount,0.0);
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot create histogram axis.",__func__,e.what());
    }
    m_yValues = std::vector<THistogramProperty>();
    for (std::size_t i = 0; i < bincount; ++i)
    {
        m_xValues[i] = xmin + (i+1)*(xmax - xmin)/bincount;
    }
    m_xValues[bincount-1] = xmax;
    m_IsConfigured = true;
}

// Configure x-axis with logarithmic scale
void MC::THistogram::ConfigureLogarithmic(const std::string& xlabel, const std::string& xunit,
    double xmin, double xmax, std::size_t bincount)
{
    if ((xmin >= xmax) || (bincount == 0) || (xmin <= 0.0))
    {
        throw EX::TInvalidStatus("Invalid histogram specifications.",__func__);
    }

    m_xLabel = xlabel;
    m_xUnit = xunit;
    m_IsLog = true;
    m_xMinimum = xmin;
    m_xMaximum = xmax;
    try
    {
        m_xValues = std::vector<double>(bincount,0.0);
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot create histogram axis.",__func__,e.what());
    }
    m_yValues = std::vector<THistogramProperty>();
    for (std::size_t i = 0; i < bincount; ++i)
    {
        m_xValues[i] = exp(log(xmin) + (i+1)*(log(xmax) - log(xmin))/bincount);
    }
    m_xValues[bincount-1] = xmax;
    m_IsConfigured = true;
}

// Add property to histogram
void MC::THistogram::AddProperty(HP prop, HT type, const std::string& label, 
    const std::string& unit, const std::string& comment)
{
    if (m_IsConfigured != true)
    {
        throw EX::TInvalidStatus("Histogram not configured.",__func__);
    }
    for (const auto& property : m_yValues)
    {
        if (prop == property.m_ID)
        {
            throw EX::TInvalidStatus("Histogram property ID already defined.",__func__);
        }
    }

    try
    {
        m_yValues.emplace_back();
        m_yValues.back().m_ID = prop;
        m_yValues.back().m_Type = type;
        m_yValues.back().m_Label = label;
        m_yValues.back().m_Unit = unit;
        m_yValues.back().m_Comment = comment;

        if (type == HT::COUNT)
        {
            m_yValues.back().m_Count = std::vector<std::uint64_t>(m_xValues.size(),0);
        }
        if ((type == HT::SUM) || (type == HT::VALUE))
        {
            m_yValues.back().m_Avg = std::vector<double>(m_xValues.size(),0.0);
        }
        if (type == HT::AVG)
        {
            m_yValues.back().m_Avg = std::vector<double>(m_xValues.size(),0.0);
            m_yValues.back().m_Count = std::vector<std::uint64_t>(m_xValues.size(),0);
        }
        if (type == HT::AVG_STDDEV)
        {
            m_yValues.back().m_Avg = std::vector<double>(m_xValues.size(),0.0);
            m_yValues.back().m_Variance = std::vector<double>(m_xValues.size(),0.0);
            m_yValues.back().m_Count = std::vector<std::uint64_t>(m_xValues.size(),0);
        }
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot create histogram property.",__func__,e.what());
    }
}

// Add value to the histogram (add to count of the appropriate bin)
void MC::THistogram::IncCount(HP prop, double xvalue, std::uint64_t inc)
{
    if (m_IsConfigured != true)
    {
        throw EX::TInvalidStatus("Histogram not configured.",__func__);
    }
    if ((xvalue < m_xMinimum) || (xvalue > m_xMaximum)) return;
    if (inc == 0) return;
    
    auto xit = std::lower_bound(m_xValues.cbegin(),m_xValues.cend(),xvalue);
    auto offset = std::distance(m_xValues.cbegin(),xit);

    auto yit = std::find_if(m_yValues.begin(),m_yValues.end(),
        [&] (const THistogramProperty& property) { return property.m_ID == prop; });
    if (yit == m_yValues.end())
    {
        throw EX::TInvalidStatus("Histogram property not found.",__func__);
    }
    if (yit->m_Type != HT::COUNT)
    {
        throw EX::TInvalidStatus("Inconsistent histogram property type.",__func__);
    }

    *(yit->m_Count.begin() + offset) += inc;
}

// Add value to the histogram (add to sum/avg/stddev or set value of the appropriate bin)
void MC::THistogram::AddValue(HP prop, double xvalue, double yvalue, std::uint64_t inc)
{
    if (m_IsConfigured != true)
    {
        throw EX::TInvalidStatus("Histogram not configured.",__func__);
    }
    if ((xvalue < m_xMinimum) || (xvalue > m_xMaximum)) return;
    if (inc == 0) return;
    
    auto xit = std::lower_bound(m_xValues.cbegin(),m_xValues.cend(),xvalue);
    auto offset = std::distance(m_xValues.cbegin(),xit);

    auto yit = std::find_if(m_yValues.begin(),m_yValues.end(),
        [&] (const THistogramProperty& property) { return property.m_ID == prop; });
    if (yit == m_yValues.end())
    {
        throw EX::TInvalidStatus("Histogram property not found.",__func__);
    }
    if (yit->m_Type == HT::COUNT)
    {
        throw EX::TInvalidStatus("Inconsistent histogram property type.",__func__);
    }

    if (yit->m_Type == HT::SUM)
    {
        *(yit->m_Avg.begin() + offset) += yvalue * static_cast<double>(inc);
    }
    if (yit->m_Type == HT::AVG)
    {
        auto avg_it = yit->m_Avg.begin() + offset;
        auto count_it = yit->m_Count.begin() + offset;

        *count_it += inc;
        *avg_it += (yvalue - *avg_it) * static_cast<double>(inc) / static_cast<double>(*count_it);
    }
    if (yit->m_Type == HT::AVG_STDDEV)
    {
        auto avg_it = yit->m_Avg.begin() + offset;
        auto var_it = yit->m_Variance.begin() + offset;
        auto count_it = yit->m_Count.begin() + offset;

        double prev_avg = *avg_it;
        *count_it += inc;
        *avg_it += (yvalue - *avg_it) * static_cast<double>(inc) / static_cast<double>(*count_it);
        *var_it += (yvalue - *avg_it) * (yvalue - prev_avg) * static_cast<double>(inc);
    }
    if (yit->m_Type == HT::VALUE)
    {
        *(yit->m_Avg.begin() + offset) = yvalue;
    }
}

// Calculate central values of bins
std::vector<double> MC::THistogram::GetBinCenters() const
{
    if (m_IsConfigured != true)
    {
        throw EX::TInvalidStatus("Histogram not configured.",__func__);
    }

    std::vector<double> xcenters (m_xValues.size(), 0.0);
    for (std::size_t i = 0; i < m_xValues.size(); ++i)
    {
        if (m_IsLog == true)
        {
            if (i == 0)
            {
                xcenters[i] = exp(0.5*(log(m_xMinimum) + log(m_xValues[i])));
            }
            else
            {
                xcenters[i] = exp(0.5*(log(m_xValues[i-1]) + log(m_xValues[i])));
            }
        }
        else
        {
            if (i == 0)
            {
                xcenters[i] = 0.5*(m_xMinimum + m_xValues[i]);
            }
            else
            {
                xcenters[i] = 0.5*(m_xValues[i-1] + m_xValues[i]);
            }
        }
    }

    return xcenters;
}

// Write histogram
void MC::THistogram::Write(std::ostream& o_str, const std::string& name) const
{
    if (m_IsConfigured != true)
    {
        throw EX::TInvalidStatus("Histogram not configured.",__func__);
    }

    try
    {
        std::vector<std::array<std::string,2>> header{
            {"min("+m_xLabel+")", (m_xUnit == "") ? "_" : m_xUnit},
            {"max("+m_xLabel+")", (m_xUnit == "") ? "_" : m_xUnit},
            {m_xLabel, (m_xUnit == "") ? "_" : m_xUnit}};
        for (const auto& prop : m_yValues)
        {
            if (prop.m_Type == HT::AVG_STDDEV)
            {
                header.push_back(std::array<std::string,2>{prop.m_Label, (prop.m_Unit == "") ? "_" : prop.m_Unit});
                header.push_back(std::array<std::string,2>{"d"+prop.m_Label, (prop.m_Unit == "") ? "_" : prop.m_Unit});
            }
            else
            {
                header.push_back(std::array<std::string,2>{prop.m_Label, (prop.m_Unit == "") ? "_" : prop.m_Unit});
            }
        }
        std::vector<std::vector<std::string>> table;
        std::stringstream sstr;
        std::vector<double> xcenters = GetBinCenters();
        for (std::size_t i = 0; i < m_xValues.size(); ++i)
        {
            table.push_back(std::vector<std::string>());
            if (i == 0)
            {
                sstr << m_xMinimum;
            }
            else
            {
                sstr << m_xValues[i-1];
            }
            table[i].push_back(sstr.str());
            sstr.str("");
            sstr << m_xValues[i];
            table[i].push_back(sstr.str());
            sstr.str("");
            sstr << xcenters[i];
            table[i].push_back(sstr.str());

            for (const auto& prop : m_yValues)
            {
                sstr.str("");
                if (prop.m_Type == HT::COUNT)
                {
                    sstr << prop.m_Count[i];
                    table[i].push_back(sstr.str());
                }
                if ((prop.m_Type == HT::SUM) || (prop.m_Type == HT::VALUE))
                {
                    sstr << prop.m_Avg[i];
                    table[i].push_back(sstr.str());
                }
                if (prop.m_Type == HT::AVG)
                {
                    if (prop.m_Count[i] != 0)
                    {
                        sstr << prop.m_Avg[i];
                    }
                    else sstr << "--";
                    table[i].push_back(sstr.str());
                }
                if (prop.m_Type == HT::AVG_STDDEV)
                {
                    if (prop.m_Count[i] != 0)
                    {
                        sstr << prop.m_Avg[i];
                    }
                    else sstr << "--";
                    table[i].push_back(sstr.str());
                    sstr.str("");
                    if (prop.m_Count[i] > 1)
                    {
                        sstr << sqrt(fabs(prop.m_Variance[i]) / static_cast<double>(prop.m_Count[i] - 1));
                    } 
                    else if (prop.m_Count[i] == 1)
                    {
                        sstr << 0.0;
                    } else sstr << "--";
                    table[i].push_back(sstr.str());
                }
            }
            sstr.str("");
        }

        o_str << "<Histogram:" << name << ">" << std::endl;
        GF::WriteTable(o_str,header,table);
        o_str << "</Histogram:" << name << ">" << std::endl;
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Cannot create histogram output.",__func__,e.what());
    }

    for (const auto& prop : m_yValues)
    {
        if (prop.m_Comment != "")
        {
            o_str << "<!-- " << prop.m_Comment << " -->" << std::endl;
        }
    }
}