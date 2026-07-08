#!/bin/bash
# Short N=2 Lohse validation runs.
#
# Submit from the repo root after building:
#   mkdir -p out/n2_validation_slurm_logs
#   sbatch slurm/n2_validation.sbatch
#
# Array task 0: ecg1d_2ninteraction
# Array task 1: ecg1d_2gaussian

#SBATCH --job-name=ecg_n2_val
#SBATCH --array=0-1
#SBATCH --partition=zen5
#SBATCH --exclude=bar
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=out/n2_validation_slurm_logs/n2_%a.out

EXES=(./build/ecg1d_2ninteraction ./build/ecg1d_2gaussian)
EXE=${EXES[$SLURM_ARRAY_TASK_ID]}

cd "$SLURM_SUBMIT_DIR" || exit 1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "[n2-validation] task $SLURM_ARRAY_TASK_ID -> $EXE on $(hostname), $OMP_NUM_THREADS threads"
exec "$EXE"
