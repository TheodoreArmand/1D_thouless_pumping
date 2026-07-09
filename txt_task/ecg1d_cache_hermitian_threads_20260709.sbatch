#!/usr/bin/env bash
#
# Submit:
#   sbatch --parsable /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/txt_task/ecg1d_cache_hermitian_threads_20260709.sbatch | tee /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_cache_hermitian_threads/ecg1d_cache_hermitian_threads_20260709/ecg1d_cache_hermitian_threads_20260709.jobid
#
# Monitor:
#   squeue -u zexuan
#   squeue -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_cache_hermitian_threads/ecg1d_cache_hermitian_threads_20260709/ecg1d_cache_hermitian_threads_20260709.jobid)"
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_cache_hermitian_threads/ecg1d_cache_hermitian_threads_20260709/ecgchthr_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_cache_hermitian_threads/ecg1d_cache_hermitian_threads_20260709/ecg1d_cache_hermitian_threads_20260709.jobid).out
#   tail -f /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_cache_hermitian_threads/ecg1d_cache_hermitian_threads_20260709/ecgchthr_$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_cache_hermitian_threads/ecg1d_cache_hermitian_threads_20260709/ecg1d_cache_hermitian_threads_20260709.jobid).err
#
# After completion:
#   sacct -j "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_cache_hermitian_threads/ecg1d_cache_hermitian_threads_20260709/ecg1d_cache_hermitian_threads_20260709.jobid)" --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
#
# Cancel:
#   scancel "$(cat /home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping/slurm/bench_cache_hermitian_threads/ecg1d_cache_hermitian_threads_20260709/ecg1d_cache_hermitian_threads_20260709.jobid)"

#SBATCH --job-name=ecgchthr
#SBATCH --account=gyqyan
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --chdir=/home/gyqyan/zexuan/ecg1d_code/1D_thouless_pumping
#SBATCH --output=slurm/bench_cache_hermitian_threads/ecg1d_cache_hermitian_threads_20260709/%x_%j.out
#SBATCH --error=slurm/bench_cache_hermitian_threads/ecg1d_cache_hermitian_threads_20260709/%x_%j.err

set -euo pipefail

tag="ecg1d_cache_hermitian_threads_20260709"
slurm_dir="slurm/bench_cache_hermitian_threads/${tag}"
out_dir="out/bench_cache_hermitian_threads/${tag}"
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
echo "binary=./build/ecg1d_cache_hermitian_validate"
echo "benchmark=K32 free/gauss, repeats=20, s=0.25"

test -x ./build/ecg1d_cache_hermitian_validate

export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=false
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

summary="${out_dir}/summary.csv"
echo "mode,threads,repeats,s,seconds_uncached,seconds_cached,speedup_cached,seconds_half_C,seconds_full_C,speedup_half_C,validation,log" > "${summary}"

run_case() {
  local mode="$1"
  local threads="$2"
  local repeats=20
  local s=0.25
  local log="${out_dir}/${mode}_threads${threads}.log"

  echo "run_start=$(date -Is) mode=${mode} threads=${threads}"
  srun --exclusive --nodes=1 --ntasks=1 --cpus-per-task="${threads}" --cpu-bind=cores \
    env OMP_NUM_THREADS="${threads}" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
      OMP_PROC_BIND=close OMP_PLACES=cores OMP_DYNAMIC=false \
      ./build/ecg1d_cache_hermitian_validate "${mode}" 32 "${repeats}" "${s}" \
    | tee "${log}"

  local seconds_uncached seconds_cached speedup_cached seconds_half seconds_full speedup_half validation
  seconds_uncached="$(awk -F= '/^seconds_per_C_plus_grad_uncached=/ {print $2}' "${log}")"
  seconds_cached="$(awk -F= '/^seconds_per_C_plus_grad_pair_cached=/ {print $2}' "${log}")"
  speedup_cached="$(awk -F= '/^speedup_cached_vs_uncached=/ {print $2}' "${log}")"
  seconds_half="$(awk -F= '/^seconds_half_C_pair_cached=/ {print $2}' "${log}")"
  seconds_full="$(awk -F= '/^seconds_full_C_reference_once=/ {print $2}' "${log}")"
  speedup_half="$(awk -F= '/^speedup_half_C_vs_full_C=/ {print $2}' "${log}")"
  validation="$(awk -F= '/^validation=/ {print $2}' "${log}")"
  printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" \
    "${mode}" "${threads}" "${repeats}" "${s}" \
    "${seconds_uncached}" "${seconds_cached}" "${speedup_cached}" \
    "${seconds_half}" "${seconds_full}" "${speedup_half}" \
    "${validation}" "${log}" >> "${summary}"
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
