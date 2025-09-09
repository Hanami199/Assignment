#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=14          
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

mkdir -p outputs_strong_scale_leo
RESULTS="outputs_strong_scale_leo/results_1_node.csv"
echo "timestamp,threads,elapsed_s,maxrss_kb" > "$RESULTS"


 export OMP_NUM_THREADS=14
 echo "Running with 14 threads"
 ts=$(date +"%Y%m%d_%H%M%S")
 log="outputs_strong_scale_leo/output_14_8Task_1Node.log"
 # time the step and append to CSV:
 /usr/bin/time -f "%e,%M" -o /tmp/time.$$ \
    srun --ntasks=32 --cpus-per-task=16 --cpu-bind=cores ./stencil_parallel -x 15000 -y 15000 -o 0 -v 1 > "$log"
 read elapsed_kb < /tmp/time.$$
 elapsed=${elapsed_kb%%,*}
 maxrss=${elapsed_kb##*,}
 echo "$ts,$nt,$elapsed,$maxrss" >> "$RESULTS"
 rm -f /tmp/time.$$

echo "Done. Logs in outputs/, summary in $RESULTS"
