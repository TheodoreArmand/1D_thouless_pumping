#!/usr/bin/env bash
#
# Submit:
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/tmp/n1_vs3_100step_20260708/n1_vs3_100step_64core.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/tmp/n1_vs3_100step_20260708/slurm/n1_vs3_100step_64core.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/tmp/n1_vs3_100step_20260708/slurm/n1_vs3_100step_64core.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/tmp/n1_vs3_100step_20260708/slurm/n1vs3_100_64_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/tmp/n1_vs3_100step_20260708/slurm/n1_vs3_100step_64core.jobid).out
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/tmp/n1_vs3_100step_20260708/slurm/n1vs3_100_64_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/tmp/n1_vs3_100step_20260708/slurm/n1_vs3_100step_64core.jobid).err
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/tmp/n1_vs3_100step_20260708/slurm/n1_vs3_100step_64core.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/tmp/n1_vs3_100step_20260708/slurm/n1_vs3_100step_64core.jobid)"

#SBATCH --job-name=n1vs3_100_64
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --exclude=bar
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=tmp/n1_vs3_100step_20260708/slurm/%x_%j.out
#SBATCH --error=tmp/n1_vs3_100step_20260708/slurm/%x_%j.err

set -euo pipefail

mkdir -p tmp/n1_vs3_100step_20260708/slurm tmp/n1_vs3_100step_20260708/out

echo "job_id=${SLURM_JOB_ID}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "submit_dir=${SLURM_SUBMIT_DIR}"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "time_limit=${SLURM_TIMELIMIT:-unset}"
echo "config=legacy_prb_3_3 N=1 K=16 dt=0.002 total_time=0.2 steps=100"
echo "output_root=tmp/n1_vs3_100step_20260708/out"
echo "start_time=$(date -Is)"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
export OMP_PROC_BIND=close
export OMP_PLACES=cores

srun ./tmp/n1_vs3_100step_20260708/bin/ecg1d_n1_vs3_100step

echo "end_time=$(date -Is)"
