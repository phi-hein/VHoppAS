#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <filesystem>
#include <cstdint>
#include <exception>

#include "source/TController.hpp"
#include "source/GlobalFunctions.hpp"
#include "source/CustomExceptions.hpp"

int main(int argc, char* argv[])
{
	std::cout << " -- VHoppAS: Variable-range Hopping in Amorphous Solids --" << std::endl; 
	std::cout << "    DOS-based kinetic Monte-Carlo simulation" << std::endl;
	std::cout << "    (by Philipp Hein, version 1.7.0";
	#ifndef NDEBUG
	std::cout << ", DEBUG build";
	#endif
	std::cout << ")" << std:: endl << std::endl;

	// Display available options when no command line arguments specified
	// (argv always contains executable name as argv[0])
	if (argc == 1)
	{
		std::cout << "Help: available command line arguments" << std::endl;
		std::cout << "-example : creates example input files in the working directory" << std::endl;
		std::cout << "-input <filename> : specify input file (relative to working directory)" << std::endl;
		std::cout << "  optionally combined with:" << std::endl;
		std::cout << "  -validate : validate complete input file (no simulation)" << std::endl;
		std::cout << "  -job <id> : execute only a certain part of all simulations, where <id> is a 1-based identification number with the following meaning:" << std::endl;
		std::cout << "              when \"ParallelizeReps = n\" (or not specified): execute all repetitions of the simulation <id> in the input parameter list" << std::endl;
		std::cout << "              when \"ParallelizeReps = y\": execute only the repetition <id> of the serialized input parameter list" << std::endl;
		std::cout << "  -collect : create summary file for all finished simulations of the input parameter list" << std::endl;
		std::cout << std::endl;
		std::cout << "Program finished." << std::endl;
		return 0;
	}

	// Parse command line arguments
	if (argc > 20)
	{
		std::cout << "Error: Too many command line arguments." << std::endl;
		return 1;
	}
	const std::vector<std::string> args(argv + 1, argv + argc);
	bool arg_example = std::find(args.begin(),args.end(),"-example") != args.end();
	bool arg_validate = std::find(args.begin(),args.end(),"-validate") != args.end();
	bool arg_collect = std::find(args.begin(),args.end(),"-collect") != args.end();
	auto it_input = std::find(args.begin(),args.end(),"-input");
	std::string arg_input_filename = "";
	if ((it_input != args.end()) && (it_input + 1 != args.end())) arg_input_filename = *(it_input + 1);
	auto it_job_id = std::find(args.begin(),args.end(),"-job");
	std::uint32_t arg_selected_job_id = 0;
	if ((it_job_id != args.end()) && ((it_job_id + 1 == args.end()) ||
		(!MC::GF::StringToUInt32(*(it_job_id + 1),arg_selected_job_id))))
	{
		std::cout << "Program aborted: INVALID JOB-ID SELECTION" << std::endl;
		return 1;
	}

	// Save runtime start
	auto runtime_start = std::chrono::steady_clock::now();

	try
	{
		// Create controller object
		MC::TController controller;

		// Write example files
		if (arg_example)
		{
			controller.GenerateExampleInputFiles();

			std::cout << "Program finished." << std::endl;
			return 0;
		}

		// Read input file
		if (arg_input_filename.empty())
		{
			std::cout << "Program aborted: NO INPUT FILE SPECIFIED" << std::endl;
			return 1;
		}
		if (!std::filesystem::exists(arg_input_filename))
		{
			std::cout << "Program aborted: INPUT FILE NOT FOUND (" << arg_input_filename << ")" << std::endl;
			return 1;
		}
		if ((arg_validate) || (arg_collect))
		{
			controller.ReadInputFile(arg_input_filename);
		}
		else
		{
			controller.ReadInputFile(arg_input_filename, arg_selected_job_id);
		}

		// Start no simulation after successful validation
		if (arg_validate)
		{
			std::cout << "Input file is valid." << std::endl;
			
			// Re-save validated input file
			controller.WriteInputFile(arg_input_filename);

			std::cout << "Program finished." << std::endl;
			return 0;
		}

		// Collect finished simulations
		if (arg_collect)
		{
			controller.CollectResults();

			std::cout << "Program finished." << std::endl;
			return 0;
		}

		// Run simulations
		controller.ExecuteSimulations();
	}
	catch (const MC::EX::TInvalidInput& e)
	{
		std::cout << "Invalid input: " << e.what() << std::endl;
		std::cout << "Program finished." << std::endl;
		return 1;
	}
	catch (const MC::EX::TFileAccess& e)
	{
		std::cout << "File access error: " << e.what() << std::endl;
		std::cout << "Program finished." << std::endl;
		return 1;
	}
	catch (const MC::EX::TOutOfMemory& e)
	{
		std::cout << "Out of memory: " << e.what() << std::endl;
		std::cout << "Program finished." << std::endl;
		return 1;
	}
	catch (const MC::EX::TInvalidStatus& e)
	{
		std::cout << "Invalid status: " << e.what() << std::endl;
		std::cout << "-- This should not happen. Report to developer! --" << std::endl;
		std::cout << "Program terminated due to inconsistency." << std::endl;
		throw;
	}
	catch (const std::exception& e)
	{
		std::cout << "Exception: " << e.what() << std::endl;
		std::cout << "-- This should not happen. Report to developer! --" << std::endl;
		std::cout << "Program terminated by exception." << std::endl;
		throw;
	}
	catch (...)
	{
		std::cout << "-- This should not happen. Report to developer! --" << std::endl;
		std::cout << "Program terminated by unknown exception." << std::endl;
		throw;
	}

	std::cout << "Program finished (total runtime: ";
	MC::GF::WriteDuration(std::cout, std::chrono::steady_clock::now() - runtime_start);
    std::cout << ")." << std::endl;
	return 0;
}