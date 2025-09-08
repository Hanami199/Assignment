#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks=32
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16          
#SBATCH --mem=0
#SBATCH --partition=EPYC
#SBATCH -t 00:15:00
#SBATCH --job-name=HPC_Exam
#SBATCH --exclusive                   

module load openMPI/5.0.5 

# compile
mpicc -march=native -O3 -std=c17 -fopenmp -Iinclude src/stencil_template_parallel.c -o stencil_parallel

# affinity 
export OMP_PLACES=cores
export OMP_PROC_BIND=close # try close later?
export OMP_DISPLAY_AFFINITY=TRUE

mkdir -p outputs_strong_scale
RESULTS="outputs_strong_scale/results_${SLURM_JOB_ID}.csv"
echo "timestamp,threads,elapsed_s,maxrss_kb" > "$RESULTS"


 export OMP_NUM_THREADS=16
 echo "Running with 16 threads"
 ts=$(date +"%Y%m%d_%H%M%S")
 log="outputs_strong_scale/output_16_8Task_Threads.log"
 # time the step and append to CSV:
 /usr/bin/time -f "%e,%M" -o /tmp/time.$$ \
    srun --ntasks=32 --cpus-per-task=16 --cpu-bind=cores ./stencil_parallel -x 15000 -y 15000 -o 0 -v 1 > "$log"
 read elapsed_kb < /tmp/time.$$
 elapsed=${elapsed_kb%%,*}
 maxrss=${elapsed_kb##*,}
 echo "$ts,$nt,$elapsed,$maxrss" >> "$RESULTS"
 rm -f /tmp/time.$$

echo "Done. Logs in outputs/, summary in $RESULTS"
