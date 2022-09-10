#include "GlobalFunctions.hpp"

#include <cstdlib>
#include <cmath>
#include <limits>
#include <regex>
#include <algorithm>
#include <stdexcept>
#include <iostream>

#include "Constants.hpp"

// Trim string from leading and trailing spaces, quotes and line terminators
void MC::GF::TrimString(std::string& str)
{
    size_t start = str.find_first_not_of(" \n\r\t\f\v\"\'");
    size_t end = str.find_last_not_of(" \n\r\t\f\v\"\'");
    if ((start == std::string::npos) || (end == std::string::npos))
    {
        str = "";
    }
    else
    {
        str = str.substr(start, end + 1 - start);
    }
}

// Convert string to double
bool MC::GF::StringToDouble(const std::string& str, double& val)
{
    try
    {
        val = std::stod(str);
    }
    catch(const std::exception& e)
    {
        #ifndef NDEBUG
        std::cerr << "Error: Conversion of string (" << str << ") to double: " << e.what() << std::endl;
        #endif
        return false;
    }
    return true;
}

// Convert string to bool
bool MC::GF::StringToBool(const std::string& str, bool& val)
{
    if (std::regex_match(str,std::regex("\\s*"+ Constant::boolex_true + "\\s*",std::regex::icase)))
    {
        val = true;
    }
    else if (std::regex_match(str,std::regex("\\s*"+ Constant::boolex_false + "\\s*",std::regex::icase)))
    {
        val = false;
    }
    else
    {
        #ifndef NDEBUG
        std::cerr << "Error: Conversion of string (" << str << ") to bool: not true or false." << std::endl;
        #endif
        return false;
    }
    return true;
}

// Convert string to int32_t
bool MC::GF::StringToInt32(const std::string& str, std::int32_t& val)
{
    try
    {
        long long_val = std::stol(str);
        std::int32_t int32_val = static_cast<std::int32_t>(long_val);
        if (long_val != int32_val)
        {
            #ifndef NDEBUG
            std::cerr << "Error: Conversion of string (" << str << ") to int32: out of range." << std::endl;
            #endif
            return false;
        }
        val = int32_val;
    }
    catch(const std::exception& e)
    {
        #ifndef NDEBUG
        std::cerr << "Error: Conversion of string (" << str << ") to int32: " << e.what() << std::endl;
        #endif
        return false;
    }
    return true;
}

// Convert string to int64_t
bool MC::GF::StringToInt64(const std::string& str, std::int64_t& val)
{
    try
    {
        long long llong_val = std::stoll(str);
        std::int64_t int64_val = static_cast<std::int64_t>(llong_val);
        if (llong_val != int64_val)
        {
            #ifndef NDEBUG
            std::cerr << "Error: Conversion of string (" << str << ") to int64: out of range." << std::endl;
            #endif
            return false;
        }
        val = int64_val;
    }
    catch(const std::exception& e)
    {
        #ifndef NDEBUG
        std::cerr << "Error: Conversion of string (" << str << ") to int64: " << e.what() << std::endl;
        #endif
        return false;
    }
    return true;
}

// Convert string to uint32_t
bool MC::GF::StringToUInt32(const std::string& str, std::uint32_t& val)
{
    try
    {
        unsigned long ulong_val = std::stoul(str);
        std::uint32_t uint32_val = static_cast<std::uint32_t>(ulong_val);
        if (ulong_val != uint32_val)
        {
            #ifndef NDEBUG
            std::cerr << "Error: Conversion of string (" << str << ") to uint32: out of range." << std::endl;
            #endif
            return false;
        }
        val = uint32_val;
    }
    catch(const std::exception& e)
    {
        #ifndef NDEBUG
        std::cerr << "Error: Conversion of string (" << str << ") to uint32: " << e.what() << std::endl;
        #endif
        return false;
    }
    return true;
}

// Convert string to uint64_t
bool MC::GF::StringToUInt64(const std::string& str, std::uint64_t& val)
{
    try
    {
        unsigned long long ullong_val = std::stoull(str);
        std::uint64_t uint64_val = static_cast<std::uint64_t>(ullong_val);
        if (ullong_val != uint64_val)
        {
            #ifndef NDEBUG
            std::cerr << "Error: Conversion of string (" << str << ") to uint64: out of range." << std::endl;
            #endif
            return false;
        }
        val = uint64_val;
    }
    catch(const std::exception& e)
    {
        #ifndef NDEBUG
        std::cerr << "Error: Conversion of string (" << str << ") to uint64: " << e.what() << std::endl;
        #endif
        return false;
    }
    return true;
}

// Compare two double values with tolerance
bool MC::GF::AlmostEqual(double x, double y, std::uint32_t ulp)
{
    if (x == y) return true;
    return (std::fabs(x-y) <= std::numeric_limits<double>::epsilon() * std::fabs(x+y) * ulp)
        || (std::fabs(x-y) < std::numeric_limits<double>::min());
}

// Combine property descriptor and unit
std::string MC::GF::CombineDescUnit(const std::array<std::string,2>& desc)
{
    if (desc[1].empty()) return desc[0];
    return desc[0] + "(" + desc[1] + ")";
}

// Contruct descriptor and unit for standard deviation
std::array<std::string,2> MC::GF::StdDevDescUnit(const std::array<std::string,2>& desc)
{
    return std::array<std::string,2>{"d" + desc[0], desc[1]};
}

// Contruct descriptor and unit for minimum
std::array<std::string,2> MC::GF::MinDescUnit(const std::array<std::string,2>& desc)
{
    return std::array<std::string,2>{"min" + desc[0], desc[1]};
}

// Contruct descriptor and unit for maximum
std::array<std::string,2> MC::GF::MaxDescUnit(const std::array<std::string,2>& desc)
{
    return std::array<std::string,2>{"max" + desc[0], desc[1]};
}

// Write a formatted table (modifies vectors with padding)
void MC::GF::WriteTable(std::ostream& o_str, std::vector<std::array<std::string,2>>& header, 
    std::vector<std::vector<std::string>>& table)
{
    if (header.empty())
    {
        return;
    }

    // Add padding to table
    std::vector<std::size_t> max_length;;
    for (auto& str : header)
    {
        max_length.push_back(std::max(str[0].length(),str[1].length())); 
    }
    for (auto& line : table)
    {
        for (std::size_t i = 0; i < line.size(); ++i)
        {
            if (line[i].length() > max_length[i]) max_length[i] = line[i].length();
        }
    }
    for (auto& line : table)
    {
        for (std::size_t i = 0; i < line.size(); ++i)
        {
            if (line[i].length() != max_length[i])
            {
                line[i].append(max_length[i] - line[i].length(),' ');
            }
        }
    }
    for (std::size_t i = 0; i < header.size(); ++i)
    {
        if (header[i][0].length() != max_length[i])
        {
            header[i][0].append(max_length[i] - header[i][0].length(),' ');
        }
        if (header[i][1].length() != max_length[i])
        {
            header[i][1].append(max_length[i] - header[i][1].length(),' ');
        }
    }
    
    // Write table
    for (auto& str : header)
    {
        o_str << str[0] << " ";
    }
    o_str << std::endl;
    for (auto& str : header)
    {
        o_str << str[1] << " ";
    }
    o_str << std::endl;
    for (auto& line : table)
    {
        for (auto& str : line)
        {
            o_str << str << " ";
        }
        o_str << std::endl;
    }
}

// Write a formatted table (modifies vectors with padding)
void MC::GF::WriteTable(std::ostream& o_str, std::vector<std::string>& header, 
    std::vector<std::vector<std::string>>& table)
{
    if (header.empty())
    {
        return;
    }

    // Add padding to table
    std::vector<std::size_t> max_length;;
    for (auto& str : header)
    {
        max_length.push_back(str.length()); 
    }
    for (auto& line : table)
    {
        for (std::size_t i = 0; i < line.size(); ++i)
        {
            if (line[i].length() > max_length[i]) max_length[i] = line[i].length();
        }
    }
    for (auto& line : table)
    {
        for (std::size_t i = 0; i < line.size(); ++i)
        {
            if (line[i].length() != max_length[i])
            {
                line[i].append(max_length[i] - line[i].length(),' ');
            }
        }
    }
    for (std::size_t i = 0; i < header.size(); ++i)
    {
        if (header[i].length() != max_length[i])
        {
            header[i].append(max_length[i] - header[i].length(),' ');
        }
    }
    
    // Write table
    for (auto& str : header)
    {
        o_str << str << " ";
    }
    o_str << std::endl;
    for (auto& line : table)
    {
        for (auto& str : line)
        {
            o_str << str << " ";
        }
        o_str << std::endl;
    }
}

// Write duration formatted as hours, minutes, seconds
void MC::GF::WriteDuration(std::ostream& o_str, const std::chrono::steady_clock::duration& duration)
{
    const auto hrs = std::chrono::duration_cast<std::chrono::hours>(duration);
    const auto mins = std::chrono::duration_cast<std::chrono::minutes>(duration - hrs);
    const auto secs = std::chrono::duration_cast<std::chrono::seconds>(duration - hrs - mins);

    if (hrs.count() != 0)
    {
        o_str << hrs.count() << " h, " << mins.count() << " min, " << secs.count() << " s";
    }
    else
    {
        if (mins.count() != 0)
        {
            o_str << mins.count() << " min, " << secs.count() << " s";
        }
        else
        {
            o_str << secs.count() << " s";
        }
    }
}