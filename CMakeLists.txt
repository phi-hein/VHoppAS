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
project(${PROJECT_NAME} VERSION 2.0.0 LANGUAGES CXX)

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

# Enable tests
enable_testing()

# Test 01.1: Single simulation without repetitions
# Test 01.2: Continue single simulation with repetitions
add_test(
  NAME Test_01 
  COMMAND sh -c "mkdir -p output && cd output && rm -f Test_01_Result* Test_01*.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -input ../Test_01_Input_01.txt > Test_01.1.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -input ../Test_01_Input_02.txt > Test_01.2.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )

# Test 02.1: Validate multi-type input file
# Test 02.2: Multiple simulations without repetitions
# Test 02.3: Adding repetitions for selected simulation
# Test 02.4: Result collection
add_test(
  NAME Test_02 
  COMMAND sh -c "mkdir -p output && cd output && rm -f Test_02_Result* Test_02*.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -validate -input ../Test_02_Input_01.txt > Test_02.1.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -input ../Test_02_Input_01.txt > Test_02.2.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 2 -input ../Test_02_Input_02.txt > Test_02.3.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -collect -input ../Test_02_Input_02.txt > Test_02.4.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )

# Test 03.1: Selected simulation with repetitions
# Test 03.2: Simulate missing repetitions
add_test(
  NAME Test_03 
  COMMAND sh -c "mkdir -p output && cd output && rm -f Test_03_Result* Test_03*.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 3 -input ../Test_03_Input.txt > Test_03.1.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -input ../Test_03_Input.txt > Test_03.2.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )

# Test 04.1: Write example files
# Test 04.2: Run single-type example simulation
# Test 04.3: Run multi-type example simulation
add_test(
  NAME Test_04 
  COMMAND sh -c "mkdir -p output && cd output && rm -f Test_04*.log Example*.txt && \
  $<TARGET_FILE:${PROJECT_NAME}> -example > Test_04.1.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -input ExampleSingle.txt > Test_04.2.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -input ExampleMulti.txt > Test_04.3.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )

# Test 05.1: Multiple repetitions with maximum verbosity
# Test 05.2: Multiple repetitions with medium verbosity
# Test 05.3: Multiple repetitions with minimum verbosity
# Test 05.4: Result collection with minimum verbosity
add_test(
  NAME Test_05 
  COMMAND sh -c "mkdir -p output && cd output && rm -f Test_05_Result* Test_05*.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 2 -input ../Test_05_Input_01.txt > Test_05.1.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 1 -input ../Test_05_Input_02.txt > Test_05.2.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 3 -input ../Test_05_Input_03.txt > Test_05.3.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -collect -input ../Test_05_Input_03.txt > Test_05.4.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )

# Test 06.1: Parallelized repetitions -> repetition without any previous results
# Test 06.2: Parallelized repetitions -> first repetition with partial results
# Test 06.3: Parallelized repetitions -> missing repetition with partial results
# Test 06.4: Parallelized repetitions -> last repetition with partial results
# Test 06.5: Result collection with partial results
# Test 06.6: Parallelized repetitions -> all missing repetitions (no selection)
# Test 06.7: Parallelized repetitions -> new repetition adjacent to previous results
# Test 06.8: Parallelized repetitions -> new repetition not adjacent to previous results
# Test 06.9: Result collection with new repetitions
add_test(
  NAME Test_06 
  COMMAND sh -c "mkdir -p output && cd output && rm -f Test_06_Result* Test_06*.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 4 -input ../Test_06_Input_01.txt > Test_06.1.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 1 -input ../Test_06_Input_01.txt > Test_06.2.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 3 -input ../Test_06_Input_01.txt > Test_06.3.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 6 -input ../Test_06_Input_01.txt > Test_06.4.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -collect -input ../Test_06_Input_01.txt > Test_06.5.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -input ../Test_06_Input_01.txt > Test_06.6.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 11 -input ../Test_06_Input_02.txt > Test_06.7.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -job 4 -input ../Test_06_Input_02.txt > Test_06.8.log && \
  $<TARGET_FILE:${PROJECT_NAME}> -collect -input ../Test_06_Input_02.txt > Test_06.9.log" 
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/test"
  )
