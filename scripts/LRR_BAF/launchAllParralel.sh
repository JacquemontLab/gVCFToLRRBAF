#!/bin/bash

# ===============================
# SLURM Job Configuration
# ===============================

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3980M
#SBATCH --account=rrg-jacquese
#SBATCH --time=01:30:00

# ===============================
# Environment Setup
# ===============================

export IPATH_NO_CPUAFFINITY=1          # Prevent Intel MPI from enforcing CPU binding
export OMP_NUM_THREADS=1               # Ensure all code runs single-threaded (no OpenMP over-subscription)

# Load required software modules
module load StdEnv/2020
module load gcc/9.3.0 openmpi/4.0.3 mpi4py/3.1.3

# ===============================
# Start Timer
# ===============================

start=`date +%s`                       # Capture the start time (UNIX timestamp)

# ===============================
# Define Input Directory
# ===============================

# Path to the directory containing the input files to process
task_directory="/home/$USER/projects/rrg-jacquese/All_user_common_folder/ANALYSIS_DIRECTORY/SeynabouWGScnvOnReadCount/SPARK/WGS_pipeline/WGS_pipeline/data"

# ===============================
# Launch Parallel Processing
# ===============================

# Use srun to run the MPI-enabled Python script
# Arguments passed to the script:
#   1. $PWD → the current working directory (used to find other scripts or config)
#   2. $task_directory → the path to the input files
#   3. $SLURM_CPUS_PER_TASK → number of CPU cores per task (passed to control internal script behavior if needed)
srun --cpus-per-task=$SLURM_CPUS_PER_TASK --cpu-bind=rank \
     python3 mpiFinalReport.py "$PWD" "$task_directory" $SLURM_CPUS_PER_TASK

# Wait until all srun-launched tasks finish
wait

# ===============================
# Calculate and Print Runtime
# ===============================

end=`date +%s`                         # Capture the end time
runtime=$((end - start))              # Compute total runtime in seconds

# Print the runtime
echo "Total execution time (in seconds): $runtime"
