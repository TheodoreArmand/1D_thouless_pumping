#!/usr/bin/env bash
#
# Submit:
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/txt_task/ecg1d_vs3_n2_64core_bench_gradcache_20260708.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708/ecg1d_vs3_n2_64core_bench_gradcache_20260708.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708/ecg1d_vs3_n2_64core_bench_gradcache_20260708.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708/vs3n2_gc64_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708/ecg1d_vs3_n2_64core_bench_gradcache_20260708.jobid).out
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708/vs3n2_gc64_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708/ecg1d_vs3_n2_64core_bench_gradcache_20260708.jobid).err
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708/ecg1d_vs3_n2_64core_bench_gradcache_20260708.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708/ecg1d_vs3_n2_64core_bench_gradcache_20260708.jobid)"

#SBATCH --job-name=vs3n2_gc64
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --exclude=bar
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708/%x_%j.out
#SBATCH --error=slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708/%x_%j.err

set -euo pipefail

mkdir -p slurm/bench_vs3_n2_64core_10step_gradcache/ecg1d_vs3_n2_64core_bench_gradcache_20260708 out/bench_vs3_n2_64core_10step

echo "job_id=${SLURM_JOB_ID}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "submit_dir=${SLURM_SUBMIT_DIR}"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun ./build/ecg1d_vs3_n2_64core_bench

echo "end_time=$(date -Is)"
