#ifndef MC_THistogram_H_
#define MC_THistogram_H_

#include <vector>
#include <string>
#include <ostream>
#include <cstdint>

#include "Constants.hpp"

namespace MC
{

struct THistogramProperty
{
    // Identification number
    HP m_ID;

    // Type
    HT m_Type;

    // Label (without white-spaces)
    std::string m_Label;

    // Unit (without white-spaces)
    std::string m_Unit;

    // Comment (written below histogram)
    std::string m_Comment;

    // Sum or average data
    std::vector<double> m_Avg;

    // Variance data
    std::vector<double> m_Variance;

    // Count data
    std::vector<std::uint64_t> m_Count;
};

class THistogram 
{
public:
    // Default constructor
    THistogram();

    // Configure x-axis with linear scale
    void ConfigureLinear(const std::string& xlabel, const std::string& xunit,
        double xmin, double xmax, std::size_t bincount);

    // Configure x-axis with logarithmic scale
    void ConfigureLogarithmic(const std::string& xlabel, const std::string& xunit,
        double xmin, double xmax, std::size_t bincount);

    // Add property to histogram
    void AddProperty(HP prop, HT type, const std::string& label, 
        const std::string& unit = "", const std::string& comment = "");

	// Add value to the histogram (add to count of the appropriate bin)
    void IncCount(HP prop, double xvalue, std::uint64_t inc = 1U);

    // Add value to the histogram (add to sum/avg/stddev or set value of the appropriate bin)
    void AddValue(HP prop, double xvalue, double yvalue, std::uint64_t inc = 1U);
    
    // Calculate central values of bins
    std::vector<double> GetBinCenters() const;

    // Write histogram
	void Write(std::ostream& o_str, const std::string& name) const;

    // Configured switch (true = scale is defined)
    bool m_IsConfigured;

    // Logarithmic switch (true = logarithmic scale)
    bool m_IsLog;

    // x axis label (without white-spaces)
	std::string m_xLabel;

    // x axis unit (without white-spaces)
	std::string m_xUnit;
    
    // x extrema
    double m_xMinimum;
    double m_xMaximum;

	// x axis values (= upper bin boundary)
	std::vector<double> m_xValues;

    // y axis values
    std::vector<THistogramProperty> m_yValues;
};

} // MC namespace
#endif  // MC_THistogram_H_