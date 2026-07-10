#!/usr/bin/env bash
#
# Submit:
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/txt_task/ecg1d_vs3_n2_grid_ref_20260708.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/ecg1d_vs3_n2_grid_ref_20260708.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/ecg1d_vs3_n2_grid_ref_20260708.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/vs3n2grid_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/ecg1d_vs3_n2_grid_ref_20260708.jobid)_0.out
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/vs3n2grid_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/ecg1d_vs3_n2_grid_ref_20260708.jobid)_1.out
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/ecg1d_vs3_n2_grid_ref_20260708.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/ecg1d_vs3_n2_grid_ref_20260708.jobid)"

#SBATCH --job-name=vs3n2grid
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=0-1
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/%x_%A_%a.out
#SBATCH --error=slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708/%x_%A_%a.err

set -euo pipefail

mkdir -p slurm/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/ecg1d_vs3_n2_grid_ref_20260708

case "${SLURM_ARRAY_TASK_ID}" in
  0)
    mode=free
    g_over_er=0.0
    ;;
  1)
    mode=gauss
    g_over_er=0.3
    ;;
  *)
    echo "unexpected SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
    ;;
esac

out_dir="out/vs3_n2_dt0p01_T160pi_full_depth_schedule/grid_ref/${mode}/g1024_dt0p01_T160pi"
mkdir -p "${out_dir}"

echo "job_id=${SLURM_JOB_ID}"
echo "array_task_id=${SLURM_ARRAY_TASK_ID}"
echo "mode=${mode}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "submit_dir=${SLURM_SUBMIT_DIR}"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "time_limit=${SLURM_TIMELIMIT:-unset}"
echo "grid=1024"
echo "dt=0.01"
echo "tmax=502.6548245743669"
echo "period=502.6548245743669"
echo "vs_over_Er=3.0"
echo "vl_over_Er=3.0"
echo "g_over_Er=${g_over_er}"
echo "fft_workers=${SLURM_CPUS_PER_TASK:-1}"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun python3 rice_mele_reference/lohes_experience/n2_grid_reference.py \
  --basis initial_state/Vs3Er_Vl3Er/initial_pathpad_N2_K24.csv \
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
  --sigma 1.0 \
  --out "${out_dir}/n2_grid_reference_vs3_${mode}_g1024_dt0p01_T160pi.npz" \
  --csv "${out_dir}/n2_grid_reference_vs3_${mode}_g1024_dt0p01_T160pi.csv"

echo "end_time=$(date -Is)"
