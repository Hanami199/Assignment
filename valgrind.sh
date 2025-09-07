#!/bin/bash
#SBATCH --job-name=valgrind_stencil
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4         # cachegrind: start single-threaded
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --hint=nomultithread
#SBATCH --exclusive

module load openMPI/5.0.5 

BIN=stencil_parallel
SRC=src/stencil_template_parallel.c
mkdir -p outputs_valgrind

# Build with debug symbols
mpicc -O3 -march=native -std=c17 -g -fno-omit-frame-pointer -fopenmp -Iinclude "$SRC" -o "$BIN"

# Keep runs small (Valgrind is 10–100× slower). Adjust your app args:
APP_ARGS="-x 512 -y 512 -o 0 -v 0 -n 100"   # <-- replace with a tiny test case if you have size flags

# Pin to one core/NUMA for stability
export OMP_NUM_THREADS=4
export OMP_PLACES=cores
export OMP_PROC_BIND=close

# ===== Cache models from your node =====
#
# Per-core caches (bytes)
I1="32768,8,64"
D1="32768,8,64"
L2="524288,8,64"
# Shared LLC slice (bytes)
L3="16777216,16,64"

mkdir -p outputs_valgrind

echo "[*] L2-view (LL = L2 = $L2)"
srun -n1 -c1 --cpu-bind=cores \
  valgrind --tool=cachegrind --cache-sim=yes --branch-sim=yes \
  --I1=$I1 --D1=$D1 --LL=$L2 \
  --log-file=outputs_valgrind/cachegrind_L2.log \
  --cachegrind-out-file=outputs_valgrind/cachegrind_L2.out \
  ./stencil_parallel $APP_ARGS

cg_annotate --auto=yes outputs_valgrind/cachegrind_L2.out \
  > outputs_valgrind/cachegrind_L2_annotate.txt

echo "[*] L3-view (LL = L3 slice = $L3)"
srun -n1 -c1 --cpu-bind=cores \
  valgrind --tool=cachegrind --cache-sim=yes --branch-sim=yes \
  --I1=$I1 --D1=$D1 --LL=$L3 \
  --log-file=outputs_valgrind/cachegrind_L3.log \
  --cachegrind-out-file=outputs_valgrind/cachegrind_L3.out \
  ./stencil_parallel $APP_ARGS

cg_annotate --auto=yes outputs_valgrind/cachegrind_L3.out \
 g> outputs_valgrind/cachegrind_L3_annotate.txt
