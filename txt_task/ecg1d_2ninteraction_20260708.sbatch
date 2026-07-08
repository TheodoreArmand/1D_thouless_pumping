#!/usr/bin/env bash
#
# Submit:
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/txt_task/ecg1d_2ninteraction_20260708.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708/ecg1d_2ninteraction_20260708.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708/ecg1d_2ninteraction_20260708.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708/ecg_2nfree_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708/ecg1d_2ninteraction_20260708.jobid).out
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708/ecg_2nfree_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708/ecg1d_2ninteraction_20260708.jobid).err
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708/ecg1d_2ninteraction_20260708.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708/ecg1d_2ninteraction_20260708.jobid)"

#SBATCH --job-name=ecg_2nfree
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708/%x_%j.out
#SBATCH --error=slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708/%x_%j.err

set -euo pipefail

mkdir -p slurm/n2_free_validation/a8p000_K24_tmax12p000_VsER5p000_VlER2p500/ecg1d_2ninteraction_20260708 out/n2_free_validation

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

srun ./build/ecg1d_2ninteraction

echo "end_time=$(date -Is)"
