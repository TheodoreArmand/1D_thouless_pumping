#!/usr/bin/env bash
#
# Submit:
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/txt_task/ecg1d_vs3_n2_k32_fullstep_threads_20260709.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_k32_fullstep_threads/ecg1d_vs3_n2_k32_fullstep_threads_20260709/ecg1d_vs3_n2_k32_fullstep_threads_20260709.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_k32_fullstep_threads/ecg1d_vs3_n2_k32_fullstep_threads_20260709/ecg1d_vs3_n2_k32_fullstep_threads_20260709.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_k32_fullstep_threads/ecg1d_vs3_n2_k32_fullstep_threads_20260709/vs3k32fs_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_k32_fullstep_threads/ecg1d_vs3_n2_k32_fullstep_threads_20260709/ecg1d_vs3_n2_k32_fullstep_threads_20260709.jobid).out
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_k32_fullstep_threads/ecg1d_vs3_n2_k32_fullstep_threads_20260709/vs3k32fs_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_k32_fullstep_threads/ecg1d_vs3_n2_k32_fullstep_threads_20260709/ecg1d_vs3_n2_k32_fullstep_threads_20260709.jobid).err
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_k32_fullstep_threads/ecg1d_vs3_n2_k32_fullstep_threads_20260709/ecg1d_vs3_n2_k32_fullstep_threads_20260709.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_vs3_n2_k32_fullstep_threads/ecg1d_vs3_n2_k32_fullstep_threads_20260709/ecg1d_vs3_n2_k32_fullstep_threads_20260709.jobid)"

#SBATCH --job-name=vs3k32fs
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/bench_vs3_n2_k32_fullstep_threads/ecg1d_vs3_n2_k32_fullstep_threads_20260709/%x_%j.out
#SBATCH --error=slurm/bench_vs3_n2_k32_fullstep_threads/ecg1d_vs3_n2_k32_fullstep_threads_20260709/%x_%j.err

set -euo pipefail

tag="ecg1d_vs3_n2_k32_fullstep_threads_20260709"
slurm_dir="slurm/bench_vs3_n2_k32_fullstep_threads/${tag}"
out_dir="out/bench_vs3_n2_k32_fullstep_threads/${tag}"
mkdir -p "${slurm_dir}" "${out_dir}"

echo "job_id=${SLURM_JOB_ID}"
echo "job_name=${SLURM_JOB_NAME}"
echo "partition=${SLURM_JOB_PARTITION}"
echo "node_list=${SLURM_JOB_NODELIST}"
echo "submit_dir=${SLURM_SUBMIT_DIR}"
echo "workdir=$(pwd)"
echo "cpus_per_task=${SLURM_CPUS_PER_TASK}"
echo "mem_per_node=${SLURM_MEM_PER_NODE:-unset}"
echo "start_time=$(date -Is)"
echo "binary=./build/ecg1d_vs3_n2_k32_fullstep_threads_bench"
echo "benchmark=K32 full RK4 TDVP step, free/gauss, steps=20, s0=0.25, dt=5e-4"

test -x ./build/ecg1d_vs3_n2_k32_fullstep_threads_bench

export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

summary="${out_dir}/summary.csv"
echo "mode,threads,steps,s0,seconds_per_full_rk4_step,wall_seconds,final_norm,log" > "${summary}"

run_case() {
  local mode="$1"
  local threads="$2"
  local steps=20
  local s0=0.25
  local log="${out_dir}/${mode}_threads${threads}.log"

  echo "run_start=$(date -Is) mode=${mode} threads=${threads}"
  srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task="${threads}" --cpu-bind=cores \
    env OMP_NUM_THREADS="${threads}" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false \
      ./build/ecg1d_vs3_n2_k32_fullstep_threads_bench "${mode}" "${steps}" "${s0}" \
    | tee "${log}"

  local seconds_per_step wall_seconds final_norm
  seconds_per_step="$(awk -F= '/^seconds_per_full_rk4_step=/ {print $2}' "${log}")"
  wall_seconds="$(awk -F= '/^wall_seconds=/ {print $2}' "${log}")"
  final_norm="$(awk -F= '/^final_norm=/ {print $2}' "${log}")"
  printf "%s,%s,%s,%s,%s,%s,%s,%s\n" \
    "${mode}" "${threads}" "${steps}" "${s0}" \
    "${seconds_per_step}" "${wall_seconds}" "${final_norm}" "${log}" >> "${summary}"
  echo "run_end=$(date -Is) mode=${mode} threads=${threads}"
}

for threads in 8 16 32; do
  run_case free "${threads}"
done

for threads in 8 16 32; do
  run_case gauss "${threads}"
done

echo "summary=${summary}"
cat "${summary}"
echo "end_time=$(date -Is)"
