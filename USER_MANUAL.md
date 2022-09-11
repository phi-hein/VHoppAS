# User manual
This document describes how the VHoppAS program can be build on the target system (where the simulations should run) and how to set up and execute simulations.

## How to build the VHoppAS program from source
Ensure that the following two components are installed on the target system (on cluster systems this may involve loading of the correct modules):

- **C++ compiler** that supports at least the **C++17** standard (e.g. GCC with version 8 or above; can be tested with `g++ --version`) 

- **CMake** (version 3.10 or above; can be tested with `cmake --version`) that is linked to the C++ compiler (e.g. by default or through the `CXX` environment variable)

&rarr; installation via package manager: `sudo apt update && sudo apt install g++ cmake`  
(see also: [gcc install](https://gcc.gnu.org/install/), [gcc binaries](https://gcc.gnu.org/install/binaries.html) or [cmake install](https://cmake.org/install/))

Download the source code of VHoppAS from [Github](https://github.com/phi-hein/VHoppAS) to a folder on the target system (unpack if necessary) and execute the following commands in this folder (which contains the `CMakeLists.txt` file):
- `mkdir build_Release`
- `cd build_Release`
- `cmake -DCMAKE_BUILD_TYPE=Release ..`
- `cmake --build .`

Recommended: Copy the resulting `VHoppAS` executable from the `build_Release` folder to the directory structure where the simulation input files reside (or set up a symbolic link).

## Input files and simulation parameters
_missing_
<!--- 
- layout of the two required input files
- meaning of the different parameters 
- meaning of SimID and RepID
- recommended directory layout (where to place exe/symlink and DOS files)
--->

## Run simulations from the command line
_missing_
<!--- 
- description of the possible command line arguments
- usage example for running a simulation
--->

## Run simulations on a computer cluster with a job scheduler
_missing_
<!--- 
- refer back to the recommended directory layout
- explain the use of array jobs with SimID selection
- usage example for SLURM cluster
--->

### Example job script on a SLURM cluster
File `job_script.sh`:
```
#!<shebang>

# Submit to job queue by command: sbatch job_script.sh
# (submit from the desired working directory)

# specify mail notification
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<email-adress>

# specify runtime limit (d-hh::mm::ss)
#SBATCH --time=1-00:00:00

# ask for memory (M is the default and can therefore be omitted, but could also be K(ilo)|G(iga)|T(era))
#SBATCH --mem-per-cpu=3800M

# name the job
#SBATCH --job-name=hpc_001

# define sub-jobs (selected simulation IDs)
#SBATCH --array=1-5,8

# declare the merged STDOUT/STDERR file
#SBATCH --output=Result(%a)-%A.out

### begin of executable commands
# print job info
date +"%a %d.%m.%Y, %T"
echo "Starting job: ${SLURM_JOB_NAME}"
echo "Job array ID: ${SLURM_ARRAY_JOB_ID}"
echo "Job array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Directory: ${SLURM_SUBMIT_DIR}"

# load runtime libraries of used compiler
module load DEVELOP
module load gcc/11
echo

# execute simulation
cd ${SLURM_SUBMIT_DIR}
../VHoppAS -sim ${SLURM_ARRAY_TASK_ID} -input Input.txt

exitcode=$?
echo
echo "Exit code: ${exitcode}"
date +"%a %d.%m.%Y, %T"
exit ${exitcode}
```

## Output files and simulation results
_missing_
<!--- 
- layout of the output files
- reference to technical notes document that explains the underlying math/definitions/algorithms
- mention the use of "Incomplete SimIDs" for restarting jobs
--->