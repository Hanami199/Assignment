#!/bin/bash
#SBATCH --nodes=1                   # 1 node
#SBATCH --ntasks=1                  # total MPI tasks across nodes
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8          # OpenMP threads per MPI task
#SBATCH --mem=0                     # use all available memory
#SBATCH --partition=EPYC

#SBATCH -t 00:00:20                 # 5 minutes for profiling and test runs
#SBATCH --job-name=threads_rank_ratio

# Set OpenMP variables
export OMP_NUM_THREADS=8
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_AFFINITY=TRUE

# Load MPI module if needed
module load openMPI/5.0.5

mpicc -O3 -g -fno-omit-frame-pointer -fopenmp -march=native \
      -Iinclude src/stencil_template_parallel.c -o stencil

# --- discover core count on this node ---
CORES=$(lscpu | awk -F: '/^CPU\(s\)/{gsub(/ /,"",$2);print $2}')
echo "Node has $CORES CPUs"

# --- OpenMP placement ---
export OMP_PLACES=cores
export OMP_PROC_BIND=close
export OMP_DISPLAY_AFFINITY=TRUE

# rank counts to try (we'll only keep those that divide CORES)
RANKS_LIST=(1 2 4 8 16 32 64 128)

OUT=results_$(date +%y%m%d_%H%M%S)
mkdir -p "$OUT"

for r in "${RANKS_LIST[@]}"; do
  # skip if more ranks than cores or doesn't divide evenly
  (( r > CORES )) && continue
  (( CORES % r != 0 )) && continue

  t=$(( CORES / r ))                  # threads per rank so r*t == CORES
  export OMP_NUM_THREADS=$t

  echo ">>> r=${r}  t=${t}"
  # time + output per config
  /usr/bin/time -f "ELAPSED=%E MAXRSS_KB=%M" -o "$OUT/time_r${r}_t${t}.txt" \
  srun --ntasks=${r} \
       --cpus-per-task=${t} \
       --cpu-bind=cores \
       ./stencil -o 0 > "$OUT/out_r${r}_t${t}.log" 2>> "$OUT/time_r${r}_t${t}.txt"
done

echo "Sweep complete. See $OUT/ for logs and timings."