#!/usr/bin/env bash
#
# Submit:
#   mkdir -p slurm/bench_legacy_prb_3_3_300step_8thread/ecg1d_legacy_prb_3_3_300step_8thread_20260709
#   sbatch --parsable txt_task/ecg1d_legacy_prb_3_3_300step_8thread_20260709.sbatch | tee slurm/bench_legacy_prb_3_3_300step_8thread/ecg1d_legacy_prb_3_3_300step_8thread_20260709/ecg1d_legacy_prb_3_3_300step_8thread_20260709.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat slurm/bench_legacy_prb_3_3_300step_8thread/ecg1d_legacy_prb_3_3_300step_8thread_20260709/ecg1d_legacy_prb_3_3_300step_8thread_20260709.jobid)"
#   tail -f slurm/bench_legacy_prb_3_3_300step_8thread/ecg1d_legacy_prb_3_3_300step_8thread_20260709/n1vs3b8_*.out
#   tail -f slurm/bench_legacy_prb_3_3_300step_8thread/ecg1d_legacy_prb_3_3_300step_8thread_20260709/n1vs3b8_*.err
#
# After completion:
#   sacct -j "$(cat slurm/bench_legacy_prb_3_3_300step_8thread/ecg1d_legacy_prb_3_3_300step_8thread_20260709/ecg1d_legacy_prb_3_3_300step_8thread_20260709.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat slurm/bench_legacy_prb_3_3_300step_8thread/ecg1d_legacy_prb_3_3_300step_8thread_20260709/ecg1d_legacy_prb_3_3_300step_8thread_20260709.jobid)"

#SBATCH --job-name=n1vs3b8
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/bench_legacy_prb_3_3_300step_8thread/ecg1d_legacy_prb_3_3_300step_8thread_20260709/%x_%j.out
#SBATCH --error=slurm/bench_legacy_prb_3_3_300step_8thread/ecg1d_legacy_prb_3_3_300step_8thread_20260709/%x_%j.err

set -euo pipefail

meta_dir="slurm/bench_legacy_prb_3_3_300step_8thread/ecg1d_legacy_prb_3_3_300step_8thread_20260709"
mkdir -p "${meta_dir}" out/bench_legacy_prb_3_3_300step_8thread

echo "job_id=${SLURM_JOB_ID}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "submit_dir=${SLURM_SUBMIT_DIR}"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "config=legacy_prb_3_3 current N=1 K=16 dt=0.002 total_time=0.6 steps=300"
echo "output_root=out/bench_legacy_prb_3_3_300step_8thread"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export MALLOC_ARENA_MAX=2

srun ./build/ecg1d_legacy_prb_3_3_300step

echo "end_time=$(date -Is)"
