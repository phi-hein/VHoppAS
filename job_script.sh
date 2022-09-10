#!/usr/local_rwth/bin/zsh

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