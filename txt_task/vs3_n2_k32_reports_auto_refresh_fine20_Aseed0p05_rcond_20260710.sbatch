#!/usr/bin/env bash
#
# Submit:
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/txt_task/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710/vs3n2htmlas_*.out
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710/vs3n2htmlas_*.err
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710.jobid)"

#SBATCH --job-name=vs3n2htmlas
#SBATCH --account=gyqyan
#SBATCH --partition=wild
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=4-00:00:00
#SBATCH --exclude=bar
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710/%x_%j.out
#SBATCH --error=slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710/%x_%j.err

set -euo pipefail

meta_dir="slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_Aseed0p05_rcond_20260710"
mkdir -p "${meta_dir}"

echo "job_id=${SLURM_JOB_ID}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "excluded_nodes=bar"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "reports=fine20_Aseed0p05_rcond1em4 fine20_Aseed0p05_rcond1em5"
echo "interval_seconds=600"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

srun python3 rice_mele_reference/refresh_vs3_n2_k32_reports.py \
  --interval 600 \
  --reports fine20_Aseed0p05_rcond1em4 fine20_Aseed0p05_rcond1em5

echo "end_time=$(date -Is)"
