#!/usr/bin/env bash
#
# Submit:
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/txt_task/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em4/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em4/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em4/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710/vs3n2cd4_*.out
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em4/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710/vs3n2cd4_*.err
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em4/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em4/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710.jobid)"

#SBATCH --job-name=vs3n2cd4
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=4-03:00:00
#SBATCH --array=0-1
#SBATCH --exclude=bar
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em4/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710/%x_%A_%a.out
#SBATCH --error=slurm/vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em4/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710/%x_%A_%a.err

set -euo pipefail

meta_dir="slurm/vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em4/ecg1d_vs3_n2_K32_fine20_COMdiag_rcond1em4_free_sigma1_20260710"
mkdir -p "${meta_dir}" out/vs3_n2_dt0p01_T160pi_K32_sweep_fine20_COMdiag_rcond1em4

case_modes=("free" "gauss")
case_sigmas=("" "1.0")
case_labels=("free" "gauss_sigma1p000")
idx="${SLURM_ARRAY_TASK_ID}"
mode="${case_modes[$idx]}"
sigma="${case_sigmas[$idx]}"
label="${case_labels[$idx]}"

echo "job_id=${SLURM_JOB_ID}"
echo "array_job_id=${SLURM_ARRAY_JOB_ID:-unset}"
echo "array_task_id=${SLURM_ARRAY_TASK_ID:-unset}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "excluded_nodes=bar"
echo "submit_dir=${SLURM_SUBMIT_DIR}"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "binary=./build/ecg1d_vs3_n2_dt001_fine20_comdiag_rcond"
echo "basis=initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32_COMdiag.csv"
echo "basis_variant=COMdiag_reserve_no_Aseed"
echo "rcond_case=rcond1em4"
echo "rcond=1e-4"
echo "case_label=${label}"
echo "mode=${mode}"
echo "sigma=${sigma}"
echo "planning_estimate=fine20 local dt in s=0.1-0.4,0.6-0.9; estimated 623294 steps; free 60.6 h, gauss 138.5 h; requested cluster maximum walltime 99 h"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false
export MALLOC_ARENA_MAX=2

if [[ "${mode}" == "free" ]]; then
  srun --cpu-bind=cores ./build/ecg1d_vs3_n2_dt001_fine20_comdiag_rcond rcond1em4 free
else
  srun --cpu-bind=cores ./build/ecg1d_vs3_n2_dt001_fine20_comdiag_rcond rcond1em4 gauss "${sigma}"
fi

echo "end_time=$(date -Is)"
