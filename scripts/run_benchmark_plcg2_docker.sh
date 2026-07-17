#!/usr/bin/env bash
# Build and run PLCg2 truncation benchmark inside Docker (default 4-CPU cap).
# Usage:
#   ./scripts/run_benchmark_plcg2_docker.sh
#   ./scripts/run_benchmark_plcg2_docker.sh --sub-frames 500 --n-windows 10 --truncation-step 100
# Optional env:
#   PHENOMS_BENCH_IMAGE        Docker image tag (default: phenoms-bench:1cpu)
#   PHENOMS_BENCH_CPUS         Docker CPU quota (default: 4)
#   PHENOMS_PLCG2_OUT_SUBDIR   Under ./phenom_outputs/benchmarks/ (default: plcg2_docker_4cpu_500f)
#   PHENOMS_PLCG2_PREP_PDB     Host path to preprocessed multi-frame PDB (optional)
#   PHENOMS_PLCG2_PREP_TRAJ    Host path to preprocessed trajectory (e.g., .xtc) (optional)
#   PHENOMS_PLCG2_PREP_TOP     Host path to topology for PREP_TRAJ (e.g., .pdb) (optional)
#   PHENOMS_SKIP_DOCKER_BUILD  1 to skip docker build if image already exists
#   PHENOMS_MOUNT_SOURCE       1 to mount local repo at /app (run latest code without rebuild)
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TAG="${PHENOMS_BENCH_IMAGE:-phenoms-bench:1cpu}"
BENCH_CPUS="${PHENOMS_BENCH_CPUS:-4}"
OUT_SUBDIR="${PHENOMS_PLCG2_OUT_SUBDIR:-plcg2_docker_4cpu_500f}"
OUT_HOST="${ROOT}/phenom_outputs/benchmarks/${OUT_SUBDIR}"
PREP_PDB="${PHENOMS_PLCG2_PREP_PDB:-}"
PREP_TRAJ="${PHENOMS_PLCG2_PREP_TRAJ:-}"
PREP_TOP="${PHENOMS_PLCG2_PREP_TOP:-}"
SKIP_BUILD="${PHENOMS_SKIP_DOCKER_BUILD:-0}"
MOUNT_SOURCE="${PHENOMS_MOUNT_SOURCE:-0}"
mkdir -p "${OUT_HOST}"

if [[ "${SKIP_BUILD}" == "1" ]]; then
  if ! docker image inspect "${TAG}" >/dev/null 2>&1; then
    echo "PHENOMS_SKIP_DOCKER_BUILD=1 but image not found: ${TAG}" >&2
    exit 1
  fi
else
  docker build -f "${ROOT}/docker/benchmark-one-cpu.Dockerfile" -t "${TAG}" "${ROOT}"
fi

ARGS=("$@")
if [[ ${#ARGS[@]} -eq 0 ]]; then
  ARGS=(--sub-frames 500 --n-windows 10 --truncation-step 100 --rust-threads 4 --mda-workers 4)
fi

MOUNTS=(-v "${OUT_HOST}:/out")
if [[ "${MOUNT_SOURCE}" == "1" ]]; then
  MOUNTS+=(-v "${ROOT}:/app")
fi
if [[ -n "${PREP_PDB}" ]]; then
  PREP_ABS="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${PREP_PDB}")"
  PREP_DIR="$(dirname "${PREP_ABS}")"
  PREP_BASE="$(basename "${PREP_ABS}")"
  MOUNTS+=(-v "${PREP_DIR}:/prep:ro")
  ARGS+=(--preprocessed-pdb "/prep/${PREP_BASE}")
elif [[ -n "${PREP_TRAJ}" || -n "${PREP_TOP}" ]]; then
  if [[ -z "${PREP_TRAJ}" || -z "${PREP_TOP}" ]]; then
    echo "Set both PHENOMS_PLCG2_PREP_TRAJ and PHENOMS_PLCG2_PREP_TOP." >&2
    exit 1
  fi
  PREP_TRAJ_ABS="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${PREP_TRAJ}")"
  PREP_TOP_ABS="$(python3 -c 'import os,sys; print(os.path.abspath(sys.argv[1]))' "${PREP_TOP}")"
  PREP_TRAJ_DIR="$(dirname "${PREP_TRAJ_ABS}")"
  PREP_TOP_DIR="$(dirname "${PREP_TOP_ABS}")"
  PREP_TRAJ_BASE="$(basename "${PREP_TRAJ_ABS}")"
  PREP_TOP_BASE="$(basename "${PREP_TOP_ABS}")"
  MOUNTS+=(-v "${PREP_TRAJ_DIR}:/prep_traj:ro")
  if [[ "${PREP_TOP_DIR}" == "${PREP_TRAJ_DIR}" ]]; then
    ARGS+=(--preprocessed-traj "/prep_traj/${PREP_TRAJ_BASE}" --preprocessed-top "/prep_traj/${PREP_TOP_BASE}")
  else
    MOUNTS+=(-v "${PREP_TOP_DIR}:/prep_top:ro")
    ARGS+=(--preprocessed-traj "/prep_traj/${PREP_TRAJ_BASE}" --preprocessed-top "/prep_top/${PREP_TOP_BASE}")
  fi
fi

exec docker run --cpus="${BENCH_CPUS}" --rm \
  -e OMP_NUM_THREADS=4 \
  -e OPENBLAS_NUM_THREADS=4 \
  -e MKL_NUM_THREADS=4 \
  -e NUMEXPR_NUM_THREADS=4 \
  -e VECLIB_MAXIMUM_THREADS=4 \
  "${MOUNTS[@]}" \
  "${TAG}" \
  python scripts/benchmark_plcg2_truncations.py \
  --out-dir /out \
  "${ARGS[@]}"
