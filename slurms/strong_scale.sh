#!/bin/sh
#SBATCH --nodes=16
#SBATCH --tasks=128
#SBATCH --cpus-per-task=14
#SBATCH --ntasks-per-node=8
#SBATCH --mem=0
#SBATCH --partition dcgp_usr_prod
#SBATCH -A uTS25_Tornator_0
#SBATCH -t 00:15:00
#SBATCH --job-name=test
#SBATCH --exclusive                                                                   
# ======================================================= 
module load openmpi/4.1.6--gcc--12.2.0

# compile
mpicc -march=native -O3 -std=c17 -fopenmp -Iinclude src/stencil_template_parallel.c -o stencil_parallel

# affinity 
export OMP_PLACES=cores
export OMP_PROC_BIND=close # try close later?
export OMP_DISPLAY_AFFINITY=TRUE
export OMP_NUM_THREADS=14

mkdir -p outputs_weak_scale_leo
RESULTS="outputs_weak_scale_leo/results_${SLURM_JOB_NUM_NODES}Node_${SLURM_NTASKS}Tasks_${OMP_NUM_THREADS}Threads.csv"
echo "timestamp,threads,elapsed_s,maxrss_kb" > "$RESULTS"

 echo "Running with 14 threads"
 ts=$(date +"%Y%m%d_%H%M%S")
 log="outputs_weak_scale_leo/output_${SLURM_JOB_NUM_NODES}Node_${SLURM_NTASKS}Tasks_${OMP_NUM_THREADS}Threads_weak.log"
 # time the step and append to CSV:
 /usr/bin/time -f "%e,%M" -o /tmp/time.$$ \
    srun --ntasks=${SLURM_NTASKS} --cpus-per-task=${OMP_NUM_THREADS} --cpu-bind=cores ./stencil_parallel -x 80000 -y 80000 -o 0 -v 1 > "$log"
 read elapsed_kb < /tmp/time.$$
 elapsed=${elapsed_kb%%,*}
 maxrss=${elapsed_kb##*,}
 echo "$ts,$nt,$elapsed,$maxrss" >> "$RESULTS"
 rm -f /tmp/time.$$

echo "Done. Logs in outputs/, summary in $RESULTS"
