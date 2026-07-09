#!/usr/bin/env bash
#
# Submit:
#   mkdir -p slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709
#   sbatch --parsable txt_task/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709.sbatch | tee slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709.jobid)"
#   tail -f slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709/%x_%j.out
#   tail -f slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709/%x_%j.err
#
# After completion:
#   sacct -j "$(cat slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709.jobid)"

#SBATCH --job-name=vs3n2htmlrc
#SBATCH --account=gyqyan
#SBATCH --partition=wild
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=4-00:00:00
#SBATCH --exclude=bar
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709/%x_%j.out
#SBATCH --error=slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709/%x_%j.err

set -euo pipefail

meta_dir="slurm/vs3_n2_reports_auto_refresh/vs3_n2_k32_reports_auto_refresh_fine20_rcond_20260709"
mkdir -p "${meta_dir}"

echo "job_id=${SLURM_JOB_ID}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "submit_dir=${SLURM_SUBMIT_DIR}"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "reports=fine20_rcond1em5,fine20_rcond1em6"
echo "refresh_interval_seconds=600"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/matplotlib-cache

srun python3 rice_mele_reference/refresh_vs3_n2_k32_reports.py --interval 600 --reports fine20_rcond1em5 fine20_rcond1em6

echo "end_time=$(date -Is)"
