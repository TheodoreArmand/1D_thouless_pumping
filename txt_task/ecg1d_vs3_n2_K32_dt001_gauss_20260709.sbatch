#!/usr/bin/env bash
#
# Submit:
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/txt_task/ecg1d_vs3_n2_K32_dt001_gauss_20260709.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709/ecg1d_vs3_n2_K32_dt001_gauss_20260709.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709/ecg1d_vs3_n2_K32_dt001_gauss_20260709.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709/vs3n2k32g_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709/ecg1d_vs3_n2_K32_dt001_gauss_20260709.jobid).out
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709/vs3n2k32g_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709/ecg1d_vs3_n2_K32_dt001_gauss_20260709.jobid).err
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709/ecg1d_vs3_n2_K32_dt001_gauss_20260709.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709/ecg1d_vs3_n2_K32_dt001_gauss_20260709.jobid)"

#SBATCH --job-name=vs3n2k32g
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=96
#SBATCH --exclude=bar
#SBATCH --mem=64G
#SBATCH --time=4-00:00:00
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709/%x_%j.out
#SBATCH --error=slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709/%x_%j.err

set -euo pipefail

mkdir -p slurm/vs3_n2_dt0p01_T160pi_K32/gauss/ecg1d_vs3_n2_K32_dt001_gauss_20260709 out/vs3_n2_dt0p01_T160pi/gauss

echo "job_id=${SLURM_JOB_ID}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "submit_dir=${SLURM_SUBMIT_DIR}"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "time_limit=${SLURM_TIMELIMIT:-unset}"
echo "planning_estimate=K32 dt0.01 with dt0.001 in s=[0.2,0.3],[0.7,0.8] gaussian: about 140747 steps, about 0.80 s/step, about 31.3 h; requested 96 h"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun ./build/ecg1d_vs3_n2_dt001_full gauss

echo "end_time=$(date -Is)"
