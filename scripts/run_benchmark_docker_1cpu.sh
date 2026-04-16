#!/usr/bin/env bash
# Build and run the kernel benchmark inside Docker (default 4-CPU cap).
# Usage:
#   ./scripts/run_benchmark_docker_1cpu.sh
#   ./scripts/run_benchmark_docker_1cpu.sh --sub-frames 250
# Optional env:
#   PHENOMS_BENCH_IMAGE   Docker image tag (default: phenoms-bench:1cpu)
#   PHENOMS_BENCH_CPUS    Docker CPU quota (default: 4)
#   PHENOMS_BENCH_OUT_SUBDIR  Under ./phenom_outputs/benchmarks/ (default: docker_4cpu)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TAG="${PHENOMS_BENCH_IMAGE:-phenoms-bench:1cpu}"
BENCH_CPUS="${PHENOMS_BENCH_CPUS:-4}"
OUT_SUBDIR="${PHENOMS_BENCH_OUT_SUBDIR:-docker_4cpu}"
OUT_HOST="${ROOT}/phenom_outputs/benchmarks/${OUT_SUBDIR}"
mkdir -p "${OUT_HOST}"

docker build -f "${ROOT}/docker/benchmark-one-cpu.Dockerfile" -t "${TAG}" "${ROOT}"

ARGS=("$@")
if [[ ${#ARGS[@]} -eq 0 ]]; then
  ARGS=(--sub-frames 15)
fi

exec docker run --cpus="${BENCH_CPUS}" --rm \
  -v "${OUT_HOST}:/out" \
  "${TAG}" \
  python scripts/benchmark_kernel.py \
  --out-dir /out \
  "${ARGS[@]}"
