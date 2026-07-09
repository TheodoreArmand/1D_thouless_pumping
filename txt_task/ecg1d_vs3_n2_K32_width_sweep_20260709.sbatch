#!/usr/bin/env bash
#
# Submit:
#   mkdir -p /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep/ecg1d_vs3_n2_K32_width_sweep_20260709
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/txt_task/ecg1d_vs3_n2_K32_width_sweep_20260709.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep/ecg1d_vs3_n2_K32_width_sweep_20260709/ecg1d_vs3_n2_K32_width_sweep_20260709.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep/ecg1d_vs3_n2_K32_width_sweep_20260709/ecg1d_vs3_n2_K32_width_sweep_20260709.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep/ecg1d_vs3_n2_K32_width_sweep_20260709/vs3n2k32sw_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep/ecg1d_vs3_n2_K32_width_sweep_20260709/ecg1d_vs3_n2_K32_width_sweep_20260709.jobid)_0.out
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep/ecg1d_vs3_n2_K32_width_sweep_20260709/ecg1d_vs3_n2_K32_width_sweep_20260709.jobid)" --format=JobID,JobName%16,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep/ecg1d_vs3_n2_K32_width_sweep_20260709/ecg1d_vs3_n2_K32_width_sweep_20260709.jobid)"

#SBATCH --job-name=vs3n2k32sw
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-7%4
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/vs3_n2_dt0p01_T160pi_K32_sweep/ecg1d_vs3_n2_K32_width_sweep_20260709/%x_%A_%a.out
#SBATCH --error=slurm/vs3_n2_dt0p01_T160pi_K32_sweep/ecg1d_vs3_n2_K32_width_sweep_20260709/%x_%A_%a.err

set -euo pipefail

mkdir -p slurm/vs3_n2_dt0p01_T160pi_K32_sweep/ecg1d_vs3_n2_K32_width_sweep_20260709
mkdir -p out/vs3_n2_dt0p01_T160pi_K32_sweep

labels=(free gauss_sigma0p500 gauss_sigma0p750 gauss_sigma1p000 gauss_sigma1p250 gauss_sigma1p500 gauss_sigma2p000 gauss_sigma3p000)
modes=(free gauss gauss gauss gauss gauss gauss gauss)
sigmas=(none 0.5 0.75 1.0 1.25 1.5 2.0 3.0)

idx="${SLURM_ARRAY_TASK_ID}"
label="${labels[$idx]}"
mode="${modes[$idx]}"
sigma="${sigmas[$idx]}"

echo "job_id=${SLURM_JOB_ID}"
echo "array_job_id=${SLURM_ARRAY_JOB_ID}"
echo "array_task_id=${SLURM_ARRAY_TASK_ID}"
echo "case_label=${label}"
echo "mode=${mode}"
echo "sigma_gauss=${sigma}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "submit_dir=${SLURM_SUBMIT_DIR}"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "time_limit=${SLURM_TIMELIMIT:-unset}"
echo "planning_estimate=K32 local-dt run, 140747 steps; free about 13.7 h, gauss conservatively about 31.3 h per width after pair-cache/eigensolver optimization"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export MALLOC_ARENA_MAX=2

if [[ "${mode}" == "free" ]]; then
  srun ./build/ecg1d_vs3_n2_dt001_full free
else
  srun ./build/ecg1d_vs3_n2_dt001_full gauss "${sigma}"
fi

echo "end_time=$(date -Is)"
