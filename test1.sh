#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128          
#SBATCH --hint=nomultithread
#SBATCH --mem=0
#SBATCH --partition=EPYC
#SBATCH -t 00:10:00
#SBATCH --job-name=HPC_Exam
#SBATCH --exclusive                   

module load openMPI/5.0.5 

# compile
mpicc -march=native -O3 -std=c17 -fopenmp -Iinclude src/stencil_template_parallel.c -o stencil_parallel

# affinity 
export OMP_PLACES=cores
export OMP_PROC_BIND=close # try close later?
export OMP_DISPLAY_AFFINITY=TRUE

mkdir -p outputs
RESULTS="outputs/results_${SLURM_JOB_ID}.csv"
echo "timestamp,threads,elapsed_s,maxrss_kb" > "$RESULTS"

for nt in 1 4 16 32 64 128; do
  export OMP_NUM_THREADS=$nt
  echo "Running with $nt threads"
  ts=$(date +"%Y%m%d_%H%M%S")
  log="outputs/output_${ts}_1Task_${nt}Threads.log"
  # time the step and append to CSV:
  /usr/bin/time -f "%e,%M" -o /tmp/time.$$ \
    srun --ntasks=1 --cpus-per-task=$nt --cpu-bind=cores ./stencil_parallel -o 0 -v 1 > "$log"
  read elapsed_kb < /tmp/time.$$
  elapsed=${elapsed_kb%%,*}
  maxrss=${elapsed_kb##*,}
  echo "$ts,$nt,$elapsed,$maxrss" >> "$RESULTS"
  rm -f /tmp/time.$$
done

echo "Done. Logs in outputs/, summary in $RESULTS"
