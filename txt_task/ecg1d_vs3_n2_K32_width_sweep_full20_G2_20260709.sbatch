#!/usr/bin/env bash
#
# Submit:
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/txt_task/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full20/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full20/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full20/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709/vs3n2fu2_*.out
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full20/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709/vs3n2fu2_*.err
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full20/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full20/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709.jobid)"

#SBATCH --job-name=vs3n2fu2
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=4-00:00:00
#SBATCH --array=0-1
#SBATCH --exclude=bar
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full20/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709/%x_%A_%a.out
#SBATCH --error=slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full20/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709/%x_%A_%a.err

set -euo pipefail

meta_dir="slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full20/ecg1d_vs3_n2_K32_width_sweep_full20_G2_20260709"
mkdir -p "${meta_dir}" out/vs3_n2_dt0p01_T160pi_K32_sweep_full20

case_sigmas=("0.75" "1.0")
case_labels=("gauss_sigma0p750" "gauss_sigma1p000")
idx="${SLURM_ARRAY_TASK_ID}"
sigma="${case_sigmas[$idx]}"
label="${case_labels[$idx]}"

echo "job_id=${SLURM_JOB_ID}"
echo "array_job_id=${SLURM_ARRAY_JOB_ID:-unset}"
echo "array_task_id=${SLURM_ARRAY_TASK_ID:-unset}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "excluded_nodes=bar"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "binary=./build/ecg1d_vs3_n2_dt0005_full"
echo "case_label=${label}"
echo "mode=gauss"
echo "sigma=${sigma}"
echo "planning_estimate=full-cycle dt/20 fixed dt=0.0005, estimated 1005310 steps; 32-core estimate gauss ~64.2 h; requested walltime 96 h"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false
export MALLOC_ARENA_MAX=2

srun --cpu-bind=cores ./build/ecg1d_vs3_n2_dt0005_full gauss "${sigma}"

echo "end_time=$(date -Is)"
