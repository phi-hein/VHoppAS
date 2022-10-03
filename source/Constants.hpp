#ifndef MC_Constants_H_
#define MC_Constants_H_

#include <string>
#include <regex>

namespace MC
{

// Output verbosity
enum class Verbosity : int
{
    MINIMUM = 0,
    MEDIUM = 1,
    MAXIMUM = 2
};

// Histogram property types
enum class HT : int
{
    COUNT,
    SUM,
    AVG,
    AVG_STDDEV,
    VALUE
};

// Histogram property IDs
enum class HP : int 
{
    ST, USED_ST, OCC_ST, OSC_ST,
    PT, USED_PT,
    EL, MOB_EL, ZEROHOP_EL, OSC_EL, ZERODISP_EL, INCOSC_EL,
    NEW_EL, NEWMOB_EL, NEWOSC_EL, NEW_INCOSC_EL,
    EFF_EL_CALC, EFF_EL_A, EFF_EL_B, EFF_EL_C,

    USED_R, OCC_R, MOB_R, OSC_R, NONOSC_R,
    OCC_FIT, OCC_RTIME,

    ST_EGY,

    HOPS, OSC, NONOSC,
    OUT_HOPS, OUT_OSC, OUT_NONOSC, IN_HOPS, IN_OSC, IN_NONOSC,

    HOPS_PER_PT, OSC_PER_PT,
    HOPS_PER_EL, OSC_PER_EL, NONOSC_PER_EL, 

    DIST_PER_PT, EDIFF_PER_PT,
    DIST_PER_HOP, EDIFF_PER_HOP,
    DIST_PER_NONOSC, EDIFF_PER_NONOSC,

    DISP_PER_EL, DISP_PER_INCOSC_EL
};

namespace Constant
{
    // Boltzmann constant (in eV/K)
    constexpr double kboltz = 8.617333E-5;

    // Elementary charge (in C)
    constexpr double echarge = 1.6021766E-19;

    // Pi
    constexpr double pi = 3.1415926535;

    // Minimum integration delta on energy axis (for interpolation; in units of kBT)
    constexpr double deltaE = 0.1;

    // Auto-adjustment of cut-offs after equilibration: Increment for used ranges in percent
    constexpr double autocutoff_inc = 5.0;

    // Regex that represents all regex meta-characters
    const std::regex regex_metachars ("\\^|\\$|\\\\|\\.|\\*|\\+|\\?|\\(|\\)|\\[|\\]|\\{|\\}|\\|");

    // Regex string that represents a double value
    const std::string dblex = "[\\+-]?(?:[1-9]\\d*|0)(?:\\.\\d+)?(?:[Ee][\\+-]?\\d+)?";

    // Regex string that represents an uint32 or uint64 value
    const std::string uintex = "\\+?\\d+";

    // Regex string that represents an int32 or int64 value
    const std::string intex = "[\\+-]?\\d+";

    // Regex string that represents a boolean value: true
    const std::string boolex_true = "(?:true|t|yes|y)";

    // Regex string that represents a boolean value: false
    const std::string boolex_false = "(?:false|f|no|n)";
}

namespace XMLSection
{
    // XML identifier: DOS
    const std::string DOS = "DOS";

    // XML identifier: Project
    const std::string Project = "MC-Project";

    // XML identifier: Input parameters
    const std::string Params = "Parameters";

    // XML identifier: Input parameters
    const std::string VariedParams = "VariedParameters";

    // XML identifier: Results
    const std::string Results = "Results";

    // XML identifier: Convergence of equilibration
    const std::string EqConv = "EquilibrationConvergence";

    // XML identifier: Convergence of simulation
    const std::string SimConv = "SimulationConvergence";

    // XML identifier: Histogram
    const std::string Histogram = "Histogram";

    // XML identifier: Table with results of individual repetitions
    const std::string SummaryTable = "Summary";

    // XML identifier: Table with average results of parameter sets (SimIDs)
    const std::string MeanTable = "Mean/StdDev-Summary";
}

} // MC namespace
#endif  // MC_Constants_H_