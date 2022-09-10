#ifndef MC_GlobalFunctions_H_
#define MC_GlobalFunctions_H_

#include <cstdint>
#include <string>
#include <vector>
#include <array>
#include <ostream>
#include <chrono>
#include <sstream>

namespace MC::GF
{

// Convert value to string (by using a stringstream)
template <class T>
std::string ToString(T val)
{
    std::stringstream sstr;
    sstr << val;
    return sstr.str();
}

// Trim string from leading and trailing spaces, quotes and line terminators
void TrimString(std::string& str);

// Convert string to double
bool StringToDouble(const std::string& str, double& val);

// Convert string to bool
bool StringToBool(const std::string& str, bool& val);

// Convert string to int32_t
bool StringToInt32(const std::string& str, std::int32_t& val);

// Convert string to int64_t
bool StringToInt64(const std::string& str, std::int64_t& val);

// Convert string to uint32_t
bool StringToUInt32(const std::string& str, std::uint32_t& val);

// Convert string to uint64_t
bool StringToUInt64(const std::string& str, std::uint64_t& val);

// Compare two double values with tolerance
bool AlmostEqual(double x, double y, std::uint32_t ulp = 10);

// Combine property descriptor and unit
std::string CombineDescUnit(const std::array<std::string,2>& desc);

// Contruct descriptor and unit for standard deviation
std::array<std::string,2> StdDevDescUnit(const std::array<std::string,2>& desc);

// Contruct descriptor and unit for minimum
std::array<std::string,2> MinDescUnit(const std::array<std::string,2>& desc);

// Contruct descriptor and unit for maximum
std::array<std::string,2> MaxDescUnit(const std::array<std::string,2>& desc);

// Write a formatted table (modifies vectors with padding)
void WriteTable(std::ostream& o_str, std::vector<std::array<std::string,2>>& header, 
    std::vector<std::vector<std::string>>& table);

// Write a formatted table (modifies vectors with padding)
void WriteTable(std::ostream& o_str, std::vector<std::string>& header, 
    std::vector<std::vector<std::string>>& table);

// Write duration formatted as hours, minutes, seconds
void WriteDuration(std::ostream& o_str, const std::chrono::steady_clock::duration& duration);

// Template class for min/max object
template <class T>
class TMinMax {
    public:
        T min, max;
        bool has_data;
        TMinMax()
            : has_data(false)
        { }
        void check(T val)
        {
            if (has_data)
            {
                if (val < min) min = val;
                if (val > max) max = val;
            }
            else
            {
                min = val;
                max = val;
                has_data = true;
            }
        }
        friend std::ostream& operator<< (std::ostream& o_str, const TMinMax& obj)
        {
            if (obj.has_data)
                if (obj.min != obj.max)
                    o_str << "[" << obj.min << ", " << obj.max << "]";
                else
                    o_str << obj.min;
            else 
                o_str << "[-- no data --]";
            return o_str;
        }
};

// Template class for min/max object with holding indices
template <class T>
class TMinMaxIdx {
    public:
        T min, max;
        std::size_t min_idx, max_idx;
        bool has_data;
        TMinMaxIdx()
            : has_data(false)
        { }
        void check(T val, std::size_t val_idx)
        {
            if (has_data)
            {
                if (val < min) 
                {
                    min = val;
                    min_idx = val_idx;
                }
                if (val > max) 
                {
                    max = val;
                    max_idx = val_idx;
                }
            }
            else
            {
                min = val;
                min_idx = val_idx;
                max = val;
                max_idx = val_idx;
                has_data = true;
            }
        }
        friend std::ostream& operator<< (std::ostream& o_str, const TMinMaxIdx& obj)
        {
            if (obj.has_data)
                if (obj.min != obj.max)
                    o_str << "[" << obj.min << ", " << obj.max << "]";
                else
                    o_str << obj.min;
            else 
                o_str << "[-- no data --]";
            return o_str;
        }
};

// Template class for min/max/mean object
template <class T>
class TMinMaxMean {
    public:
        T min, max;
        double mean;
        std::uint64_t count;
        TMinMaxMean()
            : mean(0.0), count(0)
        { }
        void check(T val, std::uint64_t inc = 1)
        {
            if (inc == 0) return;
            if (count != 0)
            {
                if (val < min) min = val;
                if (val > max) max = val;
                count += inc;
                mean += (static_cast<double>(val) - mean)*static_cast<double>(inc)/static_cast<double>(count);
            }
            else
            {
                min = val;
                max = val;
                count = inc;
                mean = static_cast<double>(val);    
            }
        }
        friend std::ostream& operator<< (std::ostream& o_str, const TMinMaxMean& obj)
        {
            if (obj.count != 0)
                if (obj.min != obj.max)
                    o_str << obj.mean << " [" << obj.min << ", " << obj.max << "]";
                else
                    o_str << obj.min;
            else 
                o_str << "[-- no data --]";
            return o_str;
        }
};

// Template class for counting values on outer parts of ranges
template <class T>
class TOuterCounter {
    public:
        double lower_bound, upper_bound;
        T lower, upper;
        TOuterCounter(double lower_boundary, double upper_boundary)
            : lower_bound(lower_boundary), upper_bound(upper_boundary), lower(0), upper(0)
        { }
        void check(double val, T inc = 1)
        {
            if (val < lower_bound) lower += inc;
            if (val > upper_bound) upper += inc;
        }
        friend std::ostream& operator<< (std::ostream& o_str, const TOuterCounter& obj)
        {
            o_str << obj.lower << " <--> " << obj.upper;
            return o_str;
        }
};

} // MC::GF namespace
#endif  // MC_GlobalFunctions_H_