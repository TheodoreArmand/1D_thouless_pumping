#!/usr/bin/env bash
#
# Submit:
#   mkdir -p slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_k32_sweep_20260709
#   sbatch --parsable txt_task/ecg1d_vs3_n2_grid_k32_sweep_20260709.sbatch | tee slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_k32_sweep_20260709/ecg1d_vs3_n2_grid_k32_sweep_20260709.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_k32_sweep_20260709/ecg1d_vs3_n2_grid_k32_sweep_20260709.jobid)"
#   tail -f slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_k32_sweep_20260709/vs3n2gk32_*.out
#   tail -f slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_k32_sweep_20260709/vs3n2gk32_*.err
#
# After completion:
#   sacct -j "$(cat slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_k32_sweep_20260709/ecg1d_vs3_n2_grid_k32_sweep_20260709.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_k32_sweep_20260709/ecg1d_vs3_n2_grid_k32_sweep_20260709.jobid)"

#SBATCH --job-name=vs3n2gk32
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH --array=0-7%2
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_k32_sweep_20260709/%x_%A_%a.out
#SBATCH --error=slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_k32_sweep_20260709/%x_%A_%a.err

set -euo pipefail

meta_dir="slurm/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_k32_sweep_20260709"
mkdir -p "${meta_dir}"

echo "job_id=${SLURM_JOB_ID}"
echo "array_job_id=${SLURM_ARRAY_JOB_ID:-unset}"
echo "array_task_id=${SLURM_ARRAY_TASK_ID:-unset}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "submit_dir=${SLURM_SUBMIT_DIR}"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export MALLOC_ARENA_MAX=2

labels=(
  "free"
  "gauss_sigma0p500"
  "gauss_sigma0p750"
  "gauss_sigma1p000"
  "gauss_sigma1p250"
  "gauss_sigma1p500"
  "gauss_sigma2p000"
  "gauss_sigma3p000"
)
g_values=("0.0" "0.3" "0.3" "0.3" "0.3" "0.3" "0.3" "0.3")
sigmas=("1.0" "0.5" "0.75" "1.0" "1.25" "1.5" "2.0" "3.0")

idx="${SLURM_ARRAY_TASK_ID}"
label="${labels[$idx]}"
g_over_er="${g_values[$idx]}"
sigma="${sigmas[$idx]}"

out_dir="out/vs3_n2_dt0p01_T160pi_K32_sweep_full_depth_schedule/grid_ref/${label}/g1024_dt0p01_T160pi"
mkdir -p "${out_dir}"

echo "case=${label}"
echo "g_over_er=${g_over_er}"
echo "sigma=${sigma}"
echo "out_dir=${out_dir}"

srun python3 rice_mele_reference/lohes_experience/n2_grid_reference.py \
  --basis initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K32.csv \
  --schedule rice_mele_reference/Vs3Vl3_3_3/gap_adaptive_vs3vl3_full_depth_schedule.csv \
  --grid 1024 \
  --cells 16 \
  --a-lat 8.0 \
  --vs-over-Er 3.0 \
  --vl-over-Er 3.0 \
  --dt 0.01 \
  --tmax 502.6548245743669 \
  --period 502.6548245743669 \
  --sample-every 25 \
  --progress-every 250 \
  --fft-workers "${SLURM_CPUS_PER_TASK:-1}" \
  --g-over-Er "${g_over_er}" \
  --sigma "${sigma}" \
  --out "${out_dir}/n2_grid_reference_vs3_${label}_g1024_dt0p01_T160pi.npz" \
  --csv "${out_dir}/n2_grid_reference_vs3_${label}_g1024_dt0p01_T160pi.csv"

echo "end_time=$(date -Is)"
