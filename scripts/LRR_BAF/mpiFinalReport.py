# ---------- Import required modules ----------
import os                     # For interacting with the OS (env vars, path joining, etc.)
import subprocess             # To execute shell scripts
import time                   # For performance timing
from mpi4py import MPI        # MPI support via mpi4py
import sys                    # For reading command-line arguments
from collections import defaultdict  # To group tasks by MPI rank

# ---------- Read command-line arguments ----------
mypath = sys.argv[1]              # Path to the directory containing the shell script
task_directory = sys.argv[2]      # Path to the folder with all sample files
cpu_per_batch = int(sys.argv[3]) # Number of CPU cores to assign to each process
                                 # (used only if taskset is active)

# ---------------------------------------------------------------------
# Function: Split a list into `x` groups using round-robin distribution
# ---------------------------------------------------------------------
def split_by_round_robin(lst, x):
    """
    Distribute items in `lst` into `x` groups as evenly as possible.
    Each group gets every x-th item (round-robin).
    """
    groups = defaultdict(list)
    for idx, item in enumerate(lst):
        group_id = idx % x
        groups[group_id].append(item)
    return [groups[i] for i in range(x)]  # Return ordered list of groups


# ---------------------------------------------------------------------
# Function: Launch script for a single sample
# ---------------------------------------------------------------------
def launchRead3(sample):

    # Retrieve the MPI local rank (rank ID on the current node)
    local_rank_str = os.environ.get("SLURM_LOCALID") or os.environ.get("OMPI_COMM_WORLD_LOCAL_RANK")
    local_rank = int(local_rank_str)

    # Compute CPU core range based on the local rank and cores assigned per rank
    start_cpu = cpu_per_batch * local_rank
    end_cpu = start_cpu + cpu_per_batch - 1
    cpu_cores = f"{start_cpu},{end_cpu}"

    # Debug output: show rank and cores
    print(f"Local rank: {local_rank}")
    print(f"Assigned CPU cores: {cpu_cores}")

    # Track execution time
    start_time = time.perf_counter()

    # Build the command to run the pipeline
    cmd = [mypath + "/launchFinalReport.sh", sample]

    # Run the command (this blocks until it finishes)
    subprocess.call(cmd)

    # Report elapsed time in minutes
    end_time = time.perf_counter()
    elapsed_minutes = (end_time - start_time) / 60
    print(f"Elapsed time: {elapsed_minutes:.3f} minutes")


# ---------------------------------------------------------------------
# Function: Assign and execute a batch of tasks per MPI rank
# ---------------------------------------------------------------------
def process_tasks(batch_rank):
    """
    Main function called by each MPI process.
    Each process gets its assigned list of tasks and processes them sequentially.
    """
    comm = MPI.COMM_WORLD
    number_of_batches = comm.Get_size()  # Total number of MPI ranks

    # List all task files from the input directory
    all_tasks = os.listdir(task_directory)

    # Distribute the tasks across all MPI ranks
    batched_tasks = split_by_round_robin(all_tasks, number_of_batches)

    # Debugging info
    print("batched_task:", batched_tasks)
    print("lots:", len(batched_tasks))
    print("processus:", number_of_batches)

    # Get the subset of tasks assigned to this process
    tasks_for_rank = batched_tasks[batch_rank]

    # Loop through and process each sample
    for task_file in tasks_for_rank:
        sample = os.path.join(task_directory, task_file)
        print("Executing finalReport on coreID:", batch_rank, "sampleID=", os.path.basename(sample))
        launchRead3(sample)

    # Indicate completion
    print("Process", batch_rank, "has finished all tasks.")


# ---------------------------------------------------------------------
# Entry point for the script
# ---------------------------------------------------------------------
if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    batch_rank = comm.Get_rank()    # Get the rank of this process
    process_tasks(batch_rank)       # Start processing assigned tasks
