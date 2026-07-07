#!/bin/bash
# dt convergence scan for the lohse_10_5 config (see txt_task/rk45_adaptive_dt_proposal.md).
#
# Nine array tasks, one per dt in {0.002 .. 0.010}. Each runs a pre-built,
# self-contained binary (out/dtscan_lohse/bin/ecg1d_dt<tag>, built by
# /tmp/build_dtscan_bins.sh) that simulates the first T = 1.5 code units of
# the Lohse pump and writes out/dtscan_lohse/dt<tag>/.
#
# Submit FROM THE REPO ROOT (binaries resolve schedule/basis paths relative
# to the working directory):
#   mkdir -p out/dtscan_lohse/slurm_logs
#   sbatch slurm/dtscan_lohse.sbatch
#
#SBATCH --job-name=ecg_dtscan
#SBATCH --array=0-8
#SBATCH --partition=zen5
#SBATCH --exclude=bar
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=out/dtscan_lohse/slurm_logs/dtscan_%a.out

DTS=(0p002 0p003 0p004 0p005 0p006 0p007 0p008 0p009 0p010)
TAG=${DTS[$SLURM_ARRAY_TASK_ID]}

cd "$SLURM_SUBMIT_DIR" || exit 1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "[dtscan] task $SLURM_ARRAY_TASK_ID -> dt tag $TAG on $(hostname), $OMP_NUM_THREADS threads"
exec ./out/dtscan_lohse/bin/ecg1d_dt$TAG
