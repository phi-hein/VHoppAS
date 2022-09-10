#include "TController.hpp"

#include <fstream>
#include <iostream>
#include <regex>
#include <algorithm>
#include <sstream>
#include <chrono>
#include <stdexcept>

#include "TEngine.hpp"
#include "Constants.hpp"
#include "TPiecewiseLinearDOS.hpp"
#include "GlobalFunctions.hpp"
#include "CustomExceptions.hpp"

// Parameter descriptors
const std::string MC::TController::s_ProjectID = "ProjectID";
const std::string MC::TController::s_ProjectName = "Name";
const std::string MC::TController::s_DOSFile = "DOS-File";
const std::string MC::TController::s_OutputFile = "Output-File";
const std::string MC::TController::s_VL = "Verbosity(0-2)";
const std::string MC::TController::s_EFTAdjust = "EFTAdjust(y/n)";
const std::string MC::TController::s_InitialFDDistrib = "InitialFDDistrib(y/n)";
const std::string MC::TController::s_TeffFit = "TeffFit(y/n)";
const std::string MC::TController::s_EnforceECount = "EnforceECount(y/n)";
const std::string MC::TController::s_CutoffAutoAdjust = "CutoffAutoAdjust(y/n)";
const std::string MC::TController::s_DistCutoffAdjustPercentage = "DistAdjustPercentage";
const std::string MC::TController::s_EdiffCutoffAdjustPercentage = "EdiffAdjustPercentage";
const std::string MC::TController::s_OnlyCompareSimID = "OnlyCompareSimID(y/n)";
const std::string MC::TController::s_UseYZVariance = "UseYZVariance(y/n)";
const std::string MC::TController::s_ProjectDescription = "Description";

// Default constructor
MC::TController::TController()
    : m_IsReady(false), m_IsFinished(false), m_InputFile(""), m_SelectedSimID(0),
    m_ProjectID(0), m_ProjectName(""), m_DOSFile(""), m_OutputFile(""), 
    m_VL(Verbosity::MAXIMUM), m_EFTAdjust(true), m_InitialFDDistrib(true), 
    m_TeffFit(true), m_EnforceECount(true), m_CutoffAutoAdjust(false), 
    m_DistCutoffAdjustPercentage(0.0), m_EdiffCutoffAdjustPercentage(0.0),
    m_OnlyCompareSimID(false), m_UseYZVariance(false),
    m_ProjectDescription("")
{

}

// Create example input files in current working directory
void MC::TController::GenerateExampleInputFiles()
{
    std::cout << "Creating input file examples ..." << std::endl;

    // DOS file
    if (m_VL >= Verbosity::MEDIUM) std::cout << "Writing DOS file (ExampleDOS.txt): " << std::flush;
    std::ofstream dos_file ("ExampleDOS.txt", std::ofstream::trunc);
  	if (dos_file.is_open())
	{
		dos_file << "<DOS>" << std::endl;
		dos_file << "Type = PiecewiseLinear" << std::endl;
        dos_file << "RefTemp(K) = 300" << std::endl;
        dos_file <<  std::endl;
        dos_file << "E-EF(eV) DOS(1/cm3eV)" << std::endl;
        dos_file << "-0.25    1.5E23" << std::endl;
        dos_file << "0.0      1.75E23" << std::endl;
        dos_file << "0.25     2.0E23" << std::endl;
		dos_file.close();
        if (m_VL >= Verbosity::MEDIUM) std::cout << "done" << std::endl;
	}
    else throw EX::TFileAccess("Cannot open/create DOS file (ExampleDOS.txt).");

    // Set example parameters
    m_IsReady = true;
    m_IsFinished = false;
    m_ProjectID = 1;
    m_ProjectName = "SingleExample";
    m_DOSFile = "ExampleDOS.txt";
    m_OutputFile = "ExampleResult.txt";  
    m_VL = Verbosity::MAXIMUM;
    m_EFTAdjust = true;
    m_InitialFDDistrib = true;
    m_TeffFit = false;
    m_EnforceECount = true;
    m_CutoffAutoAdjust = false;
    m_DistCutoffAdjustPercentage = 5.0;
    m_EdiffCutoffAdjustPercentage = 10.0;
    m_OnlyCompareSimID = false;
    m_UseYZVariance = false;
    m_ProjectDescription = "This is an example project for\n a single MC hopping simulation.";
    m_ParamSets.clear();
    m_Results.clear();
    m_ParamSets.push_back(std::unique_ptr<TParamSet>(new TParamSet()));
    m_ParamSets[0]->m_SimID = 1;
    m_ParamSets[0]->m_Repetitions = 3;
    m_ParamSets[0]->m_MinStateEnergy = -0.2;
    m_ParamSets[0]->m_MaxStateEnergy = 0.2;
	m_ParamSets[0]->m_ChemPot = 0.0;
	m_ParamSets[0]->m_StateCount = 32768U;
	m_ParamSets[0]->m_MinPathCount = 100U;
    m_ParamSets[0]->m_DistCutoff = 0.0;
    m_ParamSets[0]->m_EdiffCutoff = 0.25;
    m_ParamSets[0]->m_PreHopLimit = 1000000U;
    m_ParamSets[0]->m_EqHopLimit = 10000000U;
	m_ParamSets[0]->m_HopLimit = 10000000U;
    m_ParamSets[0]->m_RndSeed = -6526958;
	m_ParamSets[0]->m_AttemptTime = 1.0E-13;
	m_ParamSets[0]->m_Temperature = 273.15;
	m_ParamSets[0]->m_PhiGradient = 0.0;
	m_ParamSets[0]->m_LocRadius = 0.25;

    // Write single simulation file
    WriteInputFile("ExampleSingle.txt");

    // Set multi-example parameters
    m_ProjectID = 2;
    m_ProjectName = "MultiExample";
    m_ProjectDescription = "This is an example project for\n multiple MC hopping simulations.";
    m_ParamSets[0]->c_PhiGradient = false;
    m_ParamSets[0]->c_Temperature = false;
    m_ParamSets[0]->m_PhiGradient = 0.125;
    for (std::uint32_t i = 1; i < 3; ++i)
    {
        m_ParamSets.push_back(std::unique_ptr<TParamSet>(new TParamSet(*(m_ParamSets[0]))));
        m_ParamSets[i]->m_SimID = i + 1;
        m_ParamSets[i]->m_PhiGradient *= pow(10.0,i);
        m_ParamSets[i]->m_Temperature -= 20.0*i;
    }

    // Write multi simulation file
    WriteInputFile("ExampleMulti.txt");
    
    // Clear example input
    m_IsReady = false;
    m_IsFinished = false;
    m_ProjectID = 0;
    m_ProjectName = "";
    m_DOSFile = "";
    m_OutputFile = "";  
    m_VL = Verbosity::MAXIMUM;
    m_EFTAdjust = true;
    m_InitialFDDistrib = true;
    m_TeffFit = true;
    m_EnforceECount = true;
    m_CutoffAutoAdjust = false;
    m_DistCutoffAdjustPercentage = 0.0;
    m_EdiffCutoffAdjustPercentage = 0.0;
    m_OnlyCompareSimID = false;
    m_UseYZVariance = false;
    m_ProjectDescription = "";
    m_ParamSets.clear();
    m_Results.clear();
}

// Read input file (optional: select parameter set, 0 = all)
void MC::TController::ReadInputFile(const std::string& filename, const std::uint32_t sim_id)
{
    std::cout << "Reading input file: " << filename << std::endl;
    if (!std::filesystem::exists(filename))
    {
        throw EX::TInvalidInput("Input file not found.");
    }

    // Separate input file
    std::string general_str;
    std::string param_str;
    std::string varied_str;
    {
        std::string file_content;
        std::ifstream input_file (filename);
        if (input_file.is_open())
        {
            try
            {
                file_content = std::string{std::istreambuf_iterator<char>(input_file), std::istreambuf_iterator<char>()};
                input_file.close();
            }
            catch(const std::bad_alloc& e)
            {
                throw EX::TOutOfMemory("Cannot load input file.",__func__,e.what());
            }
        }
        else throw EX::TFileAccess("Cannot open input file.");
        
        std::smatch general_match;
        if (std::regex_search(file_content,general_match,std::regex("<MC-Project>((?:(?!</MC-Project>)(?:.|[\\n\\r]))+)")))
        {
            general_str = general_match[1];
        }
        std::smatch param_match;
        if (std::regex_search(file_content,param_match,std::regex("<Parameters>((?:(?!</Parameters>)(?:.|[\\n\\r]))+)")))
        {
            param_str = param_match[1];
        }
        std::smatch varied_match;
        if (std::regex_search(file_content,varied_match,std::regex("<VariedParameters>((?:.|[\\n\\r])+)")))
        {
            varied_str = varied_match[1];
        }
    }
    if ((general_str.empty()) || ((param_str.empty()) && (varied_str.empty())))
    {
        throw EX::TInvalidInput("No appropriate content detected in input file.");
    }

    // Split varied parameters into lines
    std::vector<std::string> param_lines;
    if (!varied_str.empty())
    {
        std::smatch match;
        std::string::const_iterator start(varied_str.cbegin());
        while (std::regex_search(start,varied_str.cend(),match,std::regex("[\\n\\r\\s]*([\\w\\+-].*)")))
        {
            param_lines.push_back(match[1]);  
            start = match.suffix().first;
        }
    }

    // Reset parameters
    m_IsReady = false;
    m_IsFinished = false;
    m_SelectedSimID = sim_id;
    m_DOS.reset(nullptr);
    m_ParamSets.clear();
    m_Results.clear();
    m_InputFile = filename;
    if (m_SelectedSimID == 0)
    {
        std::cout << "Parameter set selection: all" << std::endl;
    }
    else
    {
        std::cout << "Parameter set selection: " << m_SelectedSimID << std::endl;
    }

    // Read general parameters
    auto makeex = [] (const std::string& str) -> std::string
    {
        return std::regex_replace(std::regex_replace(str, std::regex("\\)"), "\\)"), std::regex("\\("), "\\(");
    };
    auto get_optional_bool = [&] (const std::string& desc, const std::string& name, bool& val, bool def_val)
    {
        std::smatch match;
        if (std::regex_search(general_str,match,std::regex(makeex(desc) + "\\s*=\\s*(\\w+)")))
        {
            if (!GF::StringToBool(match[1], val))
            {
                val = def_val;
                std::cout << name << " setting has invalid value -> set to default." << std::endl;
            }
        }
        else
        {
            val = def_val;
            std::cout << name << " setting not found -> set to default." << std::endl;
        }
    };
    auto get_optional_perc = [&] (const std::string& desc, const std::string& name, double& val, bool show_msg)
    {
        std::smatch match;
        if (std::regex_search(general_str,match,std::regex(makeex(desc) + "\\s*=\\s*(" + Constant::dblex + ")")))
        {
            if (!GF::StringToDouble(match[1], val))
            {
                val = 0.0;
                if (show_msg) std::cout << name << " setting has invalid value -> set to default." << std::endl;
            }
            else if (val <= 0.0)
            {
                val = 0.0;
                if (show_msg) std::cout << name << " setting is zero or negative -> default is used." << std::endl;
            }
        }
        else
        {
            val = 0.0;
            if (show_msg) std::cout << name << " setting not set -> default is used." << std::endl;
        }
    };
    {
        std::smatch match;
        if (std::regex_search(general_str,match,std::regex(makeex(s_ProjectID) + "\\s*=\\s*(\\d+)")))
        {
            if (!GF::StringToUInt32(match[1], m_ProjectID))
            {
                throw EX::TInvalidInput("Project-ID is no unsigned integer.");
            }
        }
        else
        {
            throw EX::TInvalidInput("Project-ID is not specified.");
        }
    }
    {
        std::smatch match;
        if (std::regex_search(general_str,match,std::regex(makeex(s_ProjectName) + "\\s*=\\s*(.+)")))
        {
            m_ProjectName = match[1];
            GF::TrimString(m_ProjectName);
        }
        else
        {
            m_ProjectName = "Default";
            std::cout << "Project name not found -> set to default." << std::endl;
        }
    }
    {
        std::smatch match;
        if (std::regex_search(general_str,match,std::regex(makeex(s_DOSFile) + "\\s*=\\s*(.+)")))
        {
            std::string str = match[1];
            GF::TrimString(str);
            m_DOSFile = str;
            if (!std::filesystem::exists(m_DOSFile))
            {
                throw EX::TInvalidInput("DOS file not found.");
            }
        }
        else
        {
            throw EX::TInvalidInput("DOS file is not specified.");
        }
    }
    {
        std::smatch match;
        if (std::regex_search(general_str,match,std::regex(makeex(s_OutputFile) + "\\s*=\\s*(.+)")))
        {
            std::string str = match[1];
            GF::TrimString(str);
            m_OutputFile = str;
            if (!m_OutputFile.has_filename())
            {
                m_OutputFile = "DefaultOutput.txt";
                std::cout << "Invalid output file name -> set to default." << std::endl;
            }
            if (!m_OutputFile.has_extension())
            {
                m_OutputFile += ".txt";
            }
        }
        else
        {
            m_OutputFile = "DefaultOutput.txt";
            std::cout << "Output file name not specified -> set to default." << std::endl;
        }
    }
    {
        std::smatch match;
        if (std::regex_search(general_str,match,std::regex(makeex(s_VL) + "\\s*=\\s*(\\d+)")))
        {
            std::uint32_t vl = 0;
            if (GF::StringToUInt32(match[1], vl))
            {
                switch (vl)
                {
                case 0:
                    m_VL = Verbosity::MINIMUM;
                    break;
                case 1:
                    m_VL = Verbosity::MEDIUM;
                    break;
                case 2:
                    m_VL = Verbosity::MAXIMUM;
                    break;
                default:
                    m_VL = Verbosity::MEDIUM;
                    std::cout << "Verbosity setting has invalid value -> set to default." << std::endl;
                    break;
                }
            }
            else
            {
                m_VL = Verbosity::MEDIUM;
                std::cout << "Verbosity setting has invalid value -> set to default." << std::endl;
            }
        }
        else
        {
            m_VL = Verbosity::MEDIUM;
            std::cout << "Verbosity setting not found -> set to default." << std::endl;
        }
    }

    get_optional_bool(s_EFTAdjust, "EFTAdjust", m_EFTAdjust, true);
    get_optional_bool(s_InitialFDDistrib, "InitialFDDistrib", m_InitialFDDistrib, true);
    get_optional_bool(s_TeffFit, "TeffFit", m_TeffFit, true);
    get_optional_bool(s_EnforceECount, "EnforceECount", m_EnforceECount, true);
    get_optional_bool(s_CutoffAutoAdjust, "CutoffAutoAdjust", m_CutoffAutoAdjust, false);
    get_optional_perc(s_DistCutoffAdjustPercentage, "DistAdjustPercentage", m_DistCutoffAdjustPercentage, m_CutoffAutoAdjust);
    get_optional_perc(s_EdiffCutoffAdjustPercentage, "EdiffAdjustPercentage", m_EdiffCutoffAdjustPercentage, m_CutoffAutoAdjust);
    get_optional_bool(s_OnlyCompareSimID, "OnlyCompareSimID", m_OnlyCompareSimID, false);
    get_optional_bool(s_UseYZVariance, "UseYZVariance", m_UseYZVariance, false);

    {
        std::smatch match;
        if (std::regex_search(general_str,match,std::regex(makeex(s_ProjectDescription) + "\\s*=\\s*((?:.|[\\n\\r])+)")))
        {
            m_ProjectDescription = match[1];
            GF::TrimString(m_ProjectDescription);
        }
        else
        {
            m_ProjectDescription = "";
        }
    }

    if (m_VL >= Verbosity::MAXIMUM) 
    {
        WriteHeader(std::cout);
        std::cout << std::endl;

        std::cout << "Clarifications:" << std::endl;
        std::cout << "- Energies refer to the energy axis of the DOS and do not include electric potential energy." << std::endl;
        std::cout << "- Energy reference: E = 0 is the position of the Fermi level in the charge neutral material." << std::endl;
        std::cout << "- Chemical potential = difference betw. the Fermi level in the simulated cell and the charge neutral Fermi level (on the energy axis of the DOS)." << std::endl;
        std::cout << "- Mobile electrons = electrons with > 0 non-oscillating hops." << std::endl;
        std::cout << std::endl;
    }

    // Read DOS file
    std::cout << "Reading DOS file: " << m_DOSFile << std::endl;
    {
        std::string dos_content;
        std::ifstream dos_file (m_DOSFile);
        if (dos_file.is_open())
        {
            try
            {
                dos_content = std::string{std::istreambuf_iterator<char>(dos_file), std::istreambuf_iterator<char>()};
                dos_file.close();
            }
            catch(const std::bad_alloc& e)
            {
                throw EX::TOutOfMemory("Cannot load DOS file.",__func__,e.what());
            }
        }
        else throw EX::TFileAccess("Cannot open DOS file.");

        if (std::regex_search(dos_content,std::regex("<DOS>[\\s\\n\\r]*Type\\s*=\\s*PiecewiseLinear")))
        {
            std::unique_ptr<MC::TPiecewiseLinearDOS> dos (new MC::TPiecewiseLinearDOS());
            dos->SpecifyDOS(dos_content);
            if (dos->HasDOS() == false) throw EX::TInvalidStatus("DOS creation failed.",__func__);
            m_DOS = std::move(dos);
            std::cout << "DOS is ready." << std::endl;
        }
        // Hint: <- insert here other TDOS sub-classes
    }
    if (!m_DOS) throw EX::TInvalidInput("No valid DOS defined.");

    // Check if unit line contains parameter set
    std::size_t first_param_line = 2;
    if (param_lines.size() >= 2)
    {
        if (m_VL >= Verbosity::MINIMUM) std::cout << "Check unit line: " << std::flush;
        std::unique_ptr<TParamSet> pset (new TParamSet());
        if (pset->Read(param_str,param_lines[0],param_lines[1],false))
        {
            first_param_line = 1;
            if (m_VL >= Verbosity::MINIMUM) std::cout << "contains parameter line" << std::endl;
        }
        else
        {
            if (m_VL >= Verbosity::MINIMUM) std::cout << "OK (not parameter line)" << std::endl;
        }
    }
    
    // Case: single parameter set
    if (param_lines.size() <= first_param_line + 1)
    {
        if (m_VL >= Verbosity::MINIMUM) std::cout << "Validating single parameter set: " << std::flush;
        std::unique_ptr<TParamSet> pset;
        try
        {
            pset = std::make_unique<TParamSet>();
        }
        catch(const std::bad_alloc& e)
        {
            throw EX::TOutOfMemory("Cannot create single parameter set",__func__,e.what());
        }
        pset->m_EFTAdjust = m_EFTAdjust;
        pset->m_InitialFDDistrib = m_InitialFDDistrib;
        pset->m_TeffFit = m_TeffFit;
        pset->m_EnforceECount = m_EnforceECount;
        pset->m_CutoffAutoAdjust = m_CutoffAutoAdjust;
        pset->m_DistCutoffAdjustPercentage = m_DistCutoffAdjustPercentage;
        pset->m_EdiffCutoffAdjustPercentage = m_EdiffCutoffAdjustPercentage;
        pset->m_OnlyCompareSimID = m_OnlyCompareSimID;
        pset->m_UseYZVariance = m_UseYZVariance;
        if (((param_lines.size() < first_param_line + 1) && (pset->Read(param_str))) ||
            ((param_lines.size() == first_param_line + 1) && (pset->Read(param_str,param_lines[0],param_lines[first_param_line]))))
        {
            TEngineData::ValidateParameters(*pset);
            m_DOS->ValidateParameters(*pset);
            m_ParamSets.push_back(std::move(pset));
            if (m_VL >= Verbosity::MINIMUM) std::cout << "OK" << std::endl;
        }
        else throw EX::TInvalidInput("Could not read single parameter set.");
    }
    // Case: selection from multiple parameter sets
    else if (m_SelectedSimID != 0)
    {
        if (m_SelectedSimID + first_param_line - 1 >= param_lines.size())
        {
            throw EX::TInvalidInput("Number of parameter sets is lower than the selected simulation ID.");
        }
          
        if (m_VL >= Verbosity::MINIMUM) std::cout << "Validating selected parameter set " 
            << m_SelectedSimID << " of " << param_lines.size() - first_param_line << ": " << std::flush;
        std::unique_ptr<TParamSet> pset;
        try
        {
            pset = std::make_unique<TParamSet>();
        }
        catch(const std::bad_alloc& e)
        {
            throw EX::TOutOfMemory("Cannot create selected parameter set",__func__,e.what());
        }
        pset->m_EFTAdjust = m_EFTAdjust;
        pset->m_InitialFDDistrib = m_InitialFDDistrib;
        pset->m_TeffFit = m_TeffFit;
        pset->m_EnforceECount = m_EnforceECount;
        pset->m_CutoffAutoAdjust = m_CutoffAutoAdjust;
        pset->m_DistCutoffAdjustPercentage = m_DistCutoffAdjustPercentage;
        pset->m_EdiffCutoffAdjustPercentage = m_EdiffCutoffAdjustPercentage;
        pset->m_OnlyCompareSimID = m_OnlyCompareSimID;
        pset->m_UseYZVariance = m_UseYZVariance;
        if (pset->Read(param_str,param_lines[0],param_lines[m_SelectedSimID + first_param_line - 1]))
        {
            TEngineData::ValidateParameters(*pset);
            m_DOS->ValidateParameters(*pset);
            m_ParamSets.push_back(std::move(pset));
            if (m_VL >= Verbosity::MINIMUM) std::cout << "OK" << std::endl;
        }
        else throw EX::TInvalidInput("Could not read selected parameter set.");
    }
    // Case: multiple parameter sets (load all)
    else
    {
        for (std::size_t i = first_param_line; i < param_lines.size(); ++i)
        {
            if (m_VL >= Verbosity::MEDIUM) std::cout << "Validating parameter set " << i - first_param_line + 1
                << " of " << param_lines.size() - first_param_line << ": " << std::flush;
            std::unique_ptr<TParamSet> pset;
            try
            {
                pset = std::make_unique<TParamSet>();
            }
            catch(const std::bad_alloc& e)
            {
                throw EX::TOutOfMemory("Cannot create parameter sets",__func__,e.what());
            }
            pset->m_EFTAdjust = m_EFTAdjust;
            pset->m_InitialFDDistrib = m_InitialFDDistrib;
            pset->m_TeffFit = m_TeffFit;
            pset->m_EnforceECount = m_EnforceECount;
            pset->m_CutoffAutoAdjust = m_CutoffAutoAdjust;
            pset->m_DistCutoffAdjustPercentage = m_DistCutoffAdjustPercentage;
            pset->m_EdiffCutoffAdjustPercentage = m_EdiffCutoffAdjustPercentage;
            pset->m_OnlyCompareSimID = m_OnlyCompareSimID;
            pset->m_UseYZVariance = m_UseYZVariance;
            if (pset->Read(param_str,param_lines[0],param_lines[i]))
            {
                TEngineData::ValidateParameters(*pset);
                m_DOS->ValidateParameters(*pset);
                m_ParamSets.push_back(std::move(pset));
                if (m_VL >= Verbosity::MINIMUM) std::cout << "OK" << std::endl;
            }
            else throw EX::TInvalidInput("Could not read parameter set.");
        }
    }
    if (m_ParamSets.empty())
    {
        throw EX::TInvalidInput("No valid simulation parameters defined.");
    }
    std::cout << std::endl;

    // Assign simulation IDs
    if (m_ParamSets.size() == 1)
    {
        if (m_SelectedSimID == 0)
        {
            m_ParamSets[0]->m_SimID = 1;
        }
        else
        {
            m_ParamSets[0]->m_SimID = m_SelectedSimID;
        }
    }
    else for (std::uint32_t i = 0; i < m_ParamSets.size(); ++i)
    {
        m_ParamSets[i]->m_SimID = i + 1;
    }

    // Create result lists for all parameter sets
    m_Results.resize(m_ParamSets.size());

    // Search already finished result files
    std::cout << "Search finished result files:" << std::endl;
    std::vector<std::filesystem::path> output_files;
    std::filesystem::path output_dir = m_OutputFile;
    output_dir.remove_filename();
    if (output_dir.empty()) output_dir = ".";
    for (auto const& dir_entry : std::filesystem::directory_iterator(output_dir)) 
    {
        if ((dir_entry.is_regular_file()) &&
            (dir_entry.path().extension() == m_OutputFile.extension()) &&
            (std::regex_match(dir_entry.path().stem().string(), std::regex(m_OutputFile.stem().string() + ".*"))))
        {
            output_files.push_back(dir_entry.path());
        }
    }
    if (output_files.empty())
        std::cout << "Output file candidates: none" << std::endl;
    else
        std::cout << "Output file candidates: " << output_files.size() << std::endl;

    // Generate and compare ParamSet and Result objects from output files
    std::size_t valid_results = 0;
    for (auto const& file_path : output_files)
    {
        // Read file content
        std::uint32_t proj_id = 0;
        std::string param_str;
        std::string result_str;
        {
            std::string file_content;
            std::ifstream output_file (file_path);
            if (output_file.is_open())
            {
                try
                {
                    file_content = std::string{std::istreambuf_iterator<char>(output_file), std::istreambuf_iterator<char>()};
                    output_file.close();
                }
                catch(const std::bad_alloc& e)
                {
                    throw EX::TOutOfMemory("Cannot load potential result file (" + file_path.string() + ").",__func__,e.what());
                }
            }
            else continue;
            
            std::smatch id_match;
            if ((!std::regex_search(file_content,id_match,std::regex(s_ProjectID + "\\s*=\\s*(\\d+)"))) ||
                (!GF::StringToUInt32(id_match[1],proj_id)) ||
                (proj_id != m_ProjectID))
            {
                continue;
            }

            std::smatch param_match;
            if (std::regex_search(file_content,param_match,std::regex("<Parameters>((?:(?!</Parameters>)(?:.|[\\n\\r]))+)")))
            {
                param_str = param_match[1];
            }
            else continue;
            
            std::smatch result_match;
            if (std::regex_search(file_content,result_match,std::regex("<Results>((?:(?!</Results>)(?:.|[\\n\\r]))+)")))
            {
                result_str = result_match[1];
            }
            else continue;
        }
        
        // Generate and compare parameter set object
        std::unique_ptr<TParamSet> pset;
        try
        {
            pset = std::make_unique<TParamSet>();
        }
        catch(const std::bad_alloc& e)
        {
            throw EX::TOutOfMemory("Cannot create parameter container.",__func__,e.what());
        }
        if (pset->Read(param_str,"","",false))
        {
            if (m_SelectedSimID != 0)
            {
                if (m_OnlyCompareSimID)
                {
                    if (pset->m_SimID != m_ParamSets[0]->m_SimID) continue;
                }
                else
                {
                    if (!(*pset == *(m_ParamSets[0]))) continue;
                }
            }
            else
            {
                if ((pset->m_SimID == 0) || (pset->m_SimID > m_ParamSets.size())) continue;
                if (m_OnlyCompareSimID)
                {
                    if (pset->m_SimID != m_ParamSets[pset->m_SimID - 1]->m_SimID) continue;
                }
                else
                {
                    if (!(*pset == *(m_ParamSets[pset->m_SimID - 1]))) continue;
                }
            }
        }
        else continue;

        // Generate result object
        std::unique_ptr<TResult> res;
        try
        {
            res = std::make_unique<TResult>();
        }
        catch(const std::bad_alloc& e)
        {
            throw EX::TOutOfMemory("Cannot create results container.",__func__,e.what());
        }
        if (!res->Read(result_str,false))
        {
            continue;
        }

        // Store the valid result
        std::cout << "Result found (Sim-ID: " << pset->m_SimID << ", Rep-ID: " << res->m_RepID 
            << "): " << file_path << std::endl;
        if (m_SelectedSimID != 0)
        {
            m_Results[0].push_back(std::move(res));
        }
        else
        {
            m_Results[pset->m_SimID - 1].push_back(std::move(res));
        }
        res = nullptr;
        valid_results++;
    }
    if (valid_results == 0)
    {
        std::cout << "Valid output files: none" << std::endl;
    }
    else
    {
        std::cout << "Valid output files: " << valid_results << std::endl;

        // Sort by Rep-ID
        for (auto& vec : m_Results)
        {
            if (vec.size() > 1)
            {
                std::sort(vec.begin(),vec.end(),
                    [](const std::unique_ptr<const TResult>& lhs,const std::unique_ptr<const TResult>& rhs) -> bool
                    {
                        return (lhs->m_RepID < rhs->m_RepID);
                    });
            }
        }

        // Write overview
        for (std::size_t i = 0; i < m_ParamSets.size(); ++i)
        {
            std::cout << "Sim-ID " << m_ParamSets[i]->m_SimID << ": " << m_Results[i].size() << " of "
                << m_ParamSets[i]->m_Repetitions << " repetitions";
            if (!m_Results[i].empty())
            {
                std::cout << " (Rep-IDs: " << m_Results[i][0]->m_RepID;
                for (std::size_t j = 1; j < m_Results[i].size(); ++j)
                {
                    std::cout << ", " << m_Results[i][j]->m_RepID;
                }
                std::cout << ")";
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;

    m_IsReady = true;
}

// Execute simulations
void MC::TController::ExecuteSimulations()
{
    if (!m_IsReady)
        throw EX::TInvalidStatus("Not ready for simulation.",__func__);
    if (m_IsFinished)
        throw EX::TInvalidStatus("Simulation is already finished.",__func__);

    // Calculate total number of simulations
    std::size_t total_sim = 0;
    for (std::size_t i = 0; i < m_ParamSets.size(); ++i)
    {
        if (m_Results[i].size() < m_ParamSets[i]->m_Repetitions)
        {
            total_sim += m_ParamSets[i]->m_Repetitions - m_Results[i].size();
        }
    }
    if (m_VL >= Verbosity::MEDIUM) std::cout << "Number of pending simulations: " << total_sim << std::endl;
    if (total_sim == 0)
    {
        std::cout << "All specified repetitions of all parameter sets are already finished." << std::endl;
        std::cout << "Increase number of repetitions to conduct more simulations." << std::endl;
        m_IsFinished = true;
        return;
    }

    // Construct summary output path
    std::filesystem::path summary_path = m_OutputFile;
    summary_path.replace_extension();
    summary_path += "_Summary";
    summary_path += m_OutputFile.extension();

    // Construct mean output path
    std::filesystem::path mean_path = m_OutputFile;
    mean_path.replace_extension();
    mean_path += "_Mean";
    mean_path += m_OutputFile.extension();
    bool write_mean_file = false;

    // Create engine object
    TEngine simulation;

    // Set verbosity
    simulation.SetVerbosity(m_VL);

	// Set DOS
    if (m_VL >= Verbosity::MEDIUM) std::cout << "Applying DOS: " << std::flush;
	simulation.SetDOS(m_DOS->Copy());
    if (m_VL >= Verbosity::MEDIUM) std::cout << "done" << std::endl << std::endl;

	// Simulations-loop
    std::size_t current_sim = 0;
    for (std::size_t i = 0; i < m_ParamSets.size(); ++i)
    {
        if ((m_SelectedSimID == 0) && (m_ParamSets[i]->m_Repetitions > 1))
        {
            write_mean_file = true;
        }
        if (m_Results[i].size() >= m_ParamSets[i]->m_Repetitions)
        {
            continue;
        }

        // Repetitions-loop
        for (std::size_t j = m_Results[i].size(); j < m_ParamSets[i]->m_Repetitions; ++j)
        {
            // Construct output file path and assign repetition ID
            std::filesystem::path output_path;
            std::filesystem::path output_path_base = m_OutputFile;
            output_path_base.replace_extension();
            output_path_base += "(" + std::to_string(m_ParamSets[i]->m_SimID) + ".";
            std::uint32_t rep_id = 0;
            if (j > 0) rep_id = m_Results[i][j - 1]->m_RepID;
            do
            {
                rep_id++;
                output_path = output_path_base;
                output_path += std::to_string(rep_id) + ")";
                output_path += m_OutputFile.extension();
            } while (std::filesystem::exists(output_path));

            // Print simulation header
            current_sim++;
            if (m_VL >= Verbosity::MEDIUM)
            {
                std::cout << " ========== SIMULATION: " << current_sim << " of " << total_sim 
                    << " (Sim-ID: " << m_ParamSets[i]->m_SimID << ", Repetition: " << j + 1 << " of " 
                    << m_ParamSets[i]->m_Repetitions << ", Rep-ID: " << rep_id << ") ========== " << std::endl;
            }
            else
            {
                std::cout << "Sim-ID " << m_ParamSets[i]->m_SimID << ", Rep-ID " << rep_id << ": ";
            }

            // Save start runtime
            auto runtime_start = std::chrono::steady_clock::now();

            // Set parameters
            simulation.SetParameters(
                std::unique_ptr<TParamSet>(new TParamSet(*(m_ParamSets[i]))), rep_id,
                m_ParamSets[i]->m_SimID * 1000U + j + 1);

            // Create random structure
            simulation.GenerateStructure();

            // Create hopping paths
            simulation.GeneratePaths();

            // Place electrons in the structure
            simulation.GenerateElectrons();

            // Run the Monte-Carlo simulation
            simulation.RunSimulation();

            // GenerateResults
            simulation.GenerateResults();

            // GetResults
            std::unique_ptr<const TResult> result;
            simulation.GetResults(result);

            // Write individual output file
            WriteOutputFile(output_path,*(m_ParamSets[i]),*result);

            // Store memory-reduced results object
            m_Results[i].push_back(TResult::ValueCopy(result.get()));

            // Write summary output file
            if ((m_SelectedSimID == 0) && (total_sim > 1))
            {
                WriteOutputFileSummary(summary_path,"");
            }

            // Write summary output file
            if (write_mean_file)
            {
                WriteOutputFileMean(mean_path,"");
            }

            // Print tests and analytics
            #ifndef NDEBUG
            if (m_VL >= Verbosity::MAXIMUM)
            {
                simulation.PrintDiagnostics();
            }
            #endif

            if (m_VL >= Verbosity::MEDIUM)
            {
                std::cout << "(total repetition runtime: ";
                GF::WriteDuration(std::cout, std::chrono::steady_clock::now() - runtime_start);
                std::cout << ")" << std::endl;
            }
            if (m_VL >= Verbosity::MINIMUM) std::cout << std::endl;
        }

        // Print combined analysis for completed simulation ID
        if ((m_VL >= Verbosity::MEDIUM) && (m_Results[i].size() > 1))
        {
            std::cout << " ========== Collective results for simulation ID: " << m_ParamSets[i]->m_SimID << " ========== " << std::endl;
            TResult::WriteMultiAnalysis(std::cout, m_Results[i]);
            std::cout << std::endl;
        }
    }

    // Print combined analysis of all simulations
    if ((m_VL >= Verbosity::MEDIUM) && (m_Results.size() > 1))
    {
        std::cout << " ========== Collective results of all simulations ========== " << std::endl;
        TResult::WriteMultiAnalysis(std::cout, m_Results);
        std::cout << std::endl;
    }

    m_IsFinished = true;
}

// Collect results of already finished simulations
void MC::TController::CollectResults()
{
    if (!m_IsReady)
        throw EX::TInvalidStatus("No valid input available.",__func__);
    if ((m_ParamSets.empty()) || (m_ParamSets.size() != m_Results.size()))
        throw EX::TInvalidStatus("Inconsistent number of parameter sets and results.",__func__);

    // Print overview of results
    std::cout << std::endl;
    std::cout << " ============ RESULTS COLLECTION ============ " << std::endl;
    if ((m_ParamSets.size() == 1) && (m_Results[0].size() == 1))
    {
        std::cout << "No collection: Only one parameter set with one result." << std::endl;
        return;
    }
    std::cout << "Number of parameter sets: " << m_ParamSets.size() << std::endl;
    std::size_t result_counter = 0;
    std::vector<std::uint32_t> incomplete_sim_ids;
    bool write_mean_file = false;
    for (std::size_t i = 0; i < m_ParamSets.size(); ++i)
    {
        std::size_t valid_count = std::count_if(m_Results[i].begin(),m_Results[i].end(),[](std::unique_ptr<const TResult>& ptr){return ptr.get() != nullptr;});
        std::cout << "Sim-ID " << m_ParamSets[i]->m_SimID << ": " << valid_count << " results for "
            << m_ParamSets[i]->m_Repetitions << " repetitions -> ";
        if (valid_count >= m_ParamSets[i]->m_Repetitions)
        {
            std::cout << "complete" << std::endl;
        }
        else
        {
            std::cout << m_ParamSets[i]->m_Repetitions - valid_count << " missing" << std::endl;
            incomplete_sim_ids.push_back(m_ParamSets[i]->m_SimID);
        }
        if (valid_count > 1) write_mean_file = true;
        result_counter += valid_count;
    }
    std::cout << "Total number of results: " << result_counter << std::endl;
    if (result_counter <= 1)
    {
        std::cout << "No collection: Too few results (zero or one)." << std::endl;
        return;
    }
    std::cout << "Incomplete simulation IDs: ";
    std::stringstream incomplete_sstr;
    if (incomplete_sim_ids.empty())
    {
        incomplete_sstr << "none";
    }
    else
    {
        std::sort(incomplete_sim_ids.begin(),incomplete_sim_ids.end());
        bool in_range = false;
        for (std::size_t i = 0; i < incomplete_sim_ids.size(); i++)
        {
            if (in_range == false)
            {
                if (i != 0) incomplete_sstr << ",";
                incomplete_sstr << incomplete_sim_ids[i];
                if (i + 1 < incomplete_sim_ids.size())
                {
                    if (incomplete_sim_ids[i] + 1 == incomplete_sim_ids[i + 1]) in_range = true;
                }
            }
            else
            {
                if (i + 1 < incomplete_sim_ids.size())
                {
                    if (incomplete_sim_ids[i] + 1 != incomplete_sim_ids[i + 1])
                    {
                        incomplete_sstr << "-" << incomplete_sim_ids[i];
                        in_range = false;
                    }
                }
                else incomplete_sstr << "-" << incomplete_sim_ids[i];
            }
        }
    }
    std::cout << incomplete_sstr.str() << std::endl;
    std::cout << std::endl;

    // Write summary output file
    std::filesystem::path summary_path = m_OutputFile;
    summary_path.replace_extension();
    summary_path += "_Summary";
    summary_path += m_OutputFile.extension();
    WriteOutputFileSummary(summary_path,incomplete_sstr.str());

    // Write mean output file
    if (write_mean_file)
    {
        std::filesystem::path mean_path = m_OutputFile;
        mean_path.replace_extension();
        mean_path += "_Mean";
        mean_path += m_OutputFile.extension();
        WriteOutputFileMean(mean_path,incomplete_sstr.str());
    }
    else std::cout << "No mean/stddev file written because only up to one result per parameter set." << std::endl;

    // Print analysis of the complete results table
    std::cout << std::endl << "Collective results of all simulations:" << std::endl;
    TResult::WriteMultiAnalysis(std::cout, m_Results);
    std::cout << std::endl;
}

// Write input file (based on m_ParamSets)
void MC::TController::WriteInputFile(const std::filesystem::path& filename) const
{
    if (!m_IsReady)
        throw EX::TInvalidStatus("No valid input available.",__func__);
    if (m_SelectedSimID != 0)
        throw EX::TInvalidStatus("Cannot write input file when simulation is selected.",__func__);
    if (m_ParamSets.empty())
        throw EX::TInvalidStatus("No parameter sets defined.",__func__);

    if (m_ParamSets.size() > 1)
    {
        if (m_VL >= Verbosity::MEDIUM) std::cout << "Writing multi-type input file (" << filename << "): " << std::flush;
    }
    else
    {
        if (m_VL >= Verbosity::MEDIUM) std::cout << "Writing single-type input file (" << filename << "): " << std::flush;
    }

    std::ofstream file (filename, std::ofstream::trunc);
  	if (file.is_open())
	{
        file << "<MC-Project>" << std::endl;
        WriteHeader(file);
        file << "</MC-Project>" << std::endl;
        file << std::endl;

        if (m_ParamSets.size() > 1)
        {
            auto header = m_ParamSets[0]->WriteTableHeader(false);
            if (!header.empty())
            {
                file << "<Parameters>" << std::endl;
                file << m_ParamSets[0]->WriteConstant();
                file << "</Parameters>" << std::endl;
                file << std::endl;

                std::vector<std::vector<std::string>> table;
                try
                {
                    for (std::size_t i = 0; i < m_ParamSets.size(); ++i)
                    {
                        table.push_back(m_ParamSets[i]->WriteTableLine(false));
                    }
                }
                catch(const std::bad_alloc& e)
                {
                    throw EX::TOutOfMemory("Parameter table too large.",__func__,e.what());
                }

                file << "<VariedParameters>" << std::endl;
                GF::WriteTable(file,header,table);
            }
            else
            {
                file << "<Parameters>" << std::endl;
                file << m_ParamSets[0]->Write(false);
                file << "</Parameters>" << std::endl;
            }
        }
        else
        {
            file << "<Parameters>" << std::endl;
            file << m_ParamSets[0]->Write(false);
            file << "</Parameters>" << std::endl;
        }
		file.close();
	}
    else throw EX::TFileAccess("Cannot open/create input file.");

    if (m_VL >= Verbosity::MEDIUM) std::cout << "done" << std::endl;
}

// Write controller parameters
void MC::TController::WriteHeader(std::ostream& o_str) const
{
    o_str << s_ProjectID << " = " << m_ProjectID << std::endl;
    o_str << s_ProjectName << " = " << m_ProjectName << std::endl;
    o_str << s_DOSFile << " = " << m_DOSFile << std::endl;
    o_str << s_OutputFile << " = " << m_OutputFile << std::endl;
    o_str << s_VL << " = ";
    if (m_VL == Verbosity::MINIMUM) o_str << "0" << std::endl;
    if (m_VL == Verbosity::MEDIUM) o_str << "1" << std::endl;
    if (m_VL == Verbosity::MAXIMUM) o_str << "2" << std::endl;
    o_str << s_EFTAdjust << " = " << ((m_EFTAdjust) ? "y" : "n") << std::endl;
    o_str << s_InitialFDDistrib << " = " << ((m_InitialFDDistrib) ? "y" : "n") << std::endl;
    o_str << s_TeffFit << " = " << ((m_TeffFit) ? "y" : "n") << std::endl;
    o_str << s_EnforceECount << " = " << ((m_EnforceECount) ? "y" : "n") << std::endl;
    o_str << s_CutoffAutoAdjust << " = " << ((m_CutoffAutoAdjust) ? "y" : "n") << std::endl;
    if ((m_DistCutoffAdjustPercentage != 0.0) || (m_CutoffAutoAdjust))
        o_str << s_DistCutoffAdjustPercentage << " = " << m_DistCutoffAdjustPercentage << std::endl;
    if ((m_EdiffCutoffAdjustPercentage != 0.0) || (m_CutoffAutoAdjust))
        o_str << s_EdiffCutoffAdjustPercentage << " = " << m_EdiffCutoffAdjustPercentage << std::endl;
    o_str << s_OnlyCompareSimID << " = " << ((m_OnlyCompareSimID) ? "y" : "n") << std::endl;
    o_str << s_UseYZVariance << " = " << ((m_UseYZVariance) ? "y" : "n") << std::endl;
    if (!m_ProjectDescription.empty())
    {
        o_str << s_ProjectDescription << " = " << m_ProjectDescription << std::endl;
    }
}

// Write output file for single simulation
void MC::TController::WriteOutputFile(const std::filesystem::path& filename, 
    const TParamSet& params, const TResult& result) const
{
    if (!m_IsReady)
        throw EX::TInvalidStatus("No valid input available.",__func__);

    if (m_VL >= Verbosity::MEDIUM) std::cout << "Writing single-type output file (" << filename << "): " << std::flush;

    std::ofstream file (filename, std::ofstream::trunc);
  	if (file.is_open())
	{
        file << "<MC-Project>" << std::endl;
        WriteHeader(file);
        file << "</MC-Project>" << std::endl;
        file << std::endl;

		file << "<Parameters>" << std::endl;
		file << params.Write();
        file << "</Parameters>" << std::endl;

        file << std::endl;
        result.Write(file);

		file.close();
	}
    else throw EX::TFileAccess("Cannot open/create output file.");

    if (m_VL >= Verbosity::MEDIUM) std::cout << "done" << std::endl;
}

// Write summary output file for all simulations (based on m_ParamSets and m_Results)
void MC::TController::WriteOutputFileSummary(const std::filesystem::path& filename, 
    const std::string& incomplete_sim_ids) const
{
    if (!m_IsReady)
        throw EX::TInvalidStatus("No valid input available.",__func__);
    if ((m_ParamSets.empty()) || (m_Results.size() != m_ParamSets.size()))
        throw EX::TInvalidStatus("Inconsistent parameter sets and results.",__func__);

    if (m_VL >= Verbosity::MEDIUM) std::cout << "Writing summary output file (" << filename << "): " << std::flush;

    // Construct table header
    auto header = m_ParamSets[0]->WriteTableHeader();
    {
        auto res_vec = TResult::WriteTableHeader();
        header.insert(header.end(),res_vec.begin(),res_vec.end());
    }

    // Construct table
    std::vector<std::vector<std::string>> table;
    try
    {
        for (std::size_t i = 0; i < m_ParamSets.size(); ++i)
        {
            auto param_line = m_ParamSets[i]->WriteTableLine();
            for (const auto& result : m_Results[i])
            {
                if (result) 
                {
                    auto res_line = result->WriteTableLine();
                    if (!res_line.empty())
                    {
                        table.push_back(param_line);
                        table.back().insert(table.back().end(),res_line.begin(),res_line.end());
                    }
                }
            }
        }
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Summary table too large.",__func__,e.what());
    }

    // Write summary file
    std::ofstream file (filename, std::ofstream::trunc);
    if (file.is_open())
    {
        file << "<MC-Project>" << std::endl;
        WriteHeader(file);
        file << "</MC-Project>" << std::endl;
        if (incomplete_sim_ids != "")
            file << "<!-- Incomplete simulation IDs: " << incomplete_sim_ids << " -->" << std::endl;
        file << std::endl;

        file << "<Parameters>" << std::endl;
        file << m_ParamSets[0]->WriteConstant();
        file << "</Parameters>" << std::endl;
        file << std::endl;

        file << "<Summary>" << std::endl;
        GF::WriteTable(file,header,table);
    }
    else throw EX::TFileAccess("Cannot open/create summary file.");

    if (m_VL >= Verbosity::MEDIUM) std::cout << "done" << std::endl;
}

// Write output file with mean of repetitions (based on m_ParamSets and m_Results)
void MC::TController::WriteOutputFileMean(const std::filesystem::path& filename, 
    const std::string& incomplete_sim_ids) const
{
    if (!m_IsReady)
        throw EX::TInvalidStatus("No valid input available.",__func__);
    if ((m_ParamSets.empty()) || (m_Results.size() != m_ParamSets.size()))
        throw EX::TInvalidStatus("Inconsistent parameter sets and results.",__func__);

    if (m_VL >= Verbosity::MEDIUM) std::cout << "Writing mean/stddev output file (" << filename << "): " << std::flush;

    // Construct table header
    auto header = m_ParamSets[0]->WriteTableHeader();
    {
        auto res_vec = TResult::WriteMultiTableHeader();
        header.insert(header.end(),res_vec.begin(),res_vec.end());
    }

    // Construct table
    std::vector<std::vector<std::string>> table;
    try
    {
        for (std::size_t i = 0; i < m_ParamSets.size(); ++i)
        {
            if (!m_Results[i].empty())
            {
                auto res_line = TResult::WriteMultiTableLine(m_Results[i]);
                if (!res_line.empty())
                {
                    table.push_back(m_ParamSets[i]->WriteTableLine());
                    table.back().insert(table.back().end(),res_line.begin(),res_line.end());
                }
            }
        }
    }
    catch(const std::bad_alloc& e)
    {
        throw EX::TOutOfMemory("Mean/stddev table too large.",__func__,e.what());
    }

    // Write mean/stddev file
    std::ofstream file (filename, std::ofstream::trunc);
    if (file.is_open())
    {
        file << "<MC-Project>" << std::endl;
        WriteHeader(file);
        file << "</MC-Project>" << std::endl;
        if (incomplete_sim_ids != "")
            file << "<!-- Incomplete simulation IDs: " << incomplete_sim_ids << " -->" << std::endl;
        file << std::endl;

        file << "<Parameters>" << std::endl;
        file << m_ParamSets[0]->WriteConstant();
        file << "</Parameters>" << std::endl;
        file << std::endl;

        file << "<Mean/StdDev-Summary>" << std::endl;
        GF::WriteTable(file,header,table);
    }
    else throw EX::TFileAccess("Cannot open/create mean/stddev file.");

    if (m_VL >= Verbosity::MEDIUM) std::cout << "done" << std::endl;
}