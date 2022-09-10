# CMake settings file

# How to build the project manually:
# (configure project by calling cmake from inside build directory with source directory (contains 
# CMakeLists.txt) as argument, then call cmake with --build option to compile the project)
# cd <build-directory>  
#   (e.g.: cd build )
# cmake -DCMAKE_BUILD_TYPE=<config> <source-directory>  
#   (e.g.: cmake -DCMAKE_BUILD_TYPE=Debug .. )
#   (e.g.: cmake -DCMAKE_BUILD_TYPE=Release .. ) 
# cmake --build <build-directory>  
#   (e.g.: cmake --build . )

# Require CMake 3.10
cmake_minimum_required(VERSION 3.10)

# Specify project name and version
set(PROJECT_NAME VHoppAS)
project(${PROJECT_NAME} VERSION 1.6.0 LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Create list of source files (additional to main.cpp; header files are optional = just for IDEs)
set(MY_SOURCE_FILES
  source/TController.hpp
  source/TController.cpp
  source/TEngine.hpp
  source/TEngine.cpp
  source/TEngineData.hpp
  source/TEngineData.cpp
  source/Constants.hpp
  source/GlobalFunctions.hpp
  source/GlobalFunctions.cpp
  source/CustomExceptions.hpp
  source/CustomExceptions.cpp
  source/TDOS.hpp
  source/TPiecewiseLinearDOS.hpp
  source/TPiecewiseLinearDOS.cpp
  source/TElectron.hpp
  source/TElectron.cpp
  source/TParamSet.hpp
  source/TParamSet.cpp
  source/TLocalState.hpp
  source/TLocalState.cpp
  source/TPath.hpp
  source/TPath.cpp
  source/TResult.hpp
  source/TResult.cpp
  source/THistogram.hpp
  source/THistogram.cpp
)

# Executable
add_executable(${PROJECT_NAME} 
  main.cpp
  ${MY_SOURCE_FILES}
)

# Specify include directory
target_include_directories(${PROJECT_NAME} PUBLIC "${CMAKE_SOURCE_DIR}/source")

# Specify tests
enable_testing()
add_test(
  NAME Remove_result_and_log_files
  COMMAND sh -c "rm -f TestResult* Test*.log Example*.txt" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )
add_test(
  NAME Test01_Single_simulation_without_repetitions 
  COMMAND sh -c "$<TARGET_FILE:${PROJECT_NAME}> -input TestInput_01.txt > Test_01.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )
add_test(
  NAME Test02_Continue_single_simulation_with_repetitions 
  COMMAND sh -c "$<TARGET_FILE:${PROJECT_NAME}> -input TestInput_02.txt > Test_02.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )
add_test(
  NAME Test03.1_Validate_multi_type_input_file 
  COMMAND sh -c "$<TARGET_FILE:${PROJECT_NAME}> -validate -input TestInput_03.txt > Test_03.1.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )
add_test(
  NAME Test03.2_Multiple_simulations_without_repetitions
  COMMAND sh -c "$<TARGET_FILE:${PROJECT_NAME}> -input TestInput_03.txt > Test_03.2.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )
add_test(
  NAME Test04.1_Adding_repetitions_for_selected_simulation
  COMMAND sh -c "$<TARGET_FILE:${PROJECT_NAME}> -sim 2 -input TestInput_04.txt > Test_04.1.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )
add_test(
  NAME Test04.2_Result_collection
  COMMAND sh -c "$<TARGET_FILE:${PROJECT_NAME}> -collect -input TestInput_04.txt > Test_04.2.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )
add_test(
  NAME Test05.1_Selected_simulation_with_repetitions 
  COMMAND sh -c "$<TARGET_FILE:${PROJECT_NAME}> -sim 3 -input TestInput_05.txt > Test_05.1.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )
add_test(
  NAME Test05.2_Simulate_missing_repetitions 
  COMMAND sh -c "$<TARGET_FILE:${PROJECT_NAME}> -input TestInput_05.txt > Test_05.2.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )
add_test(
  NAME Test06_Write_example_files 
  COMMAND sh -c "$<TARGET_FILE:${PROJECT_NAME}> -example > Test_06.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )

set_tests_properties(
  Remove_result_and_log_files
  Test01_Single_simulation_without_repetitions
  Test02_Continue_single_simulation_with_repetitions
  Test03.1_Validate_multi_type_input_file
  Test03.2_Multiple_simulations_without_repetitions
  Test04.1_Adding_repetitions_for_selected_simulation
  Test04.2_Result_collection
  Test05.1_Selected_simulation_with_repetitions
  Test05.2_Simulate_missing_repetitions 
  Test06_Write_example_files 
  PROPERTIES RUN_SERIAL TRUE)