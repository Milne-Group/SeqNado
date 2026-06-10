#!/usr/bin/env bash
# Build SeqNado containers for linux/amd64 and linux/arm64 using docker buildx.
#
# Usage:
#   ./containers/build.sh [--push] [--tag TAG]
#
# Options:
#   --push       Push images to registry after build (requires --tag)
#   --tag TAG    Image name prefix, e.g. ghcr.io/asmith/seqnado (default: seqnado)
#
# Examples:
#   ./containers/build.sh
#   ./containers/build.sh --tag ghcr.io/asmith/seqnado --push

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PLATFORMS="linux/amd64,linux/arm64"
TAG_PREFIX="seqnado"
PUSH=0

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    --push) PUSH=1; shift ;;
    --tag)  TAG_PREFIX="$2"; shift 2 ;;
    *) echo "Unknown option: $1" >&2; exit 1 ;;
  esac
done

# Ensure buildx is available
if ! docker buildx version &>/dev/null; then
  echo "Error: docker buildx is required. Update Docker Desktop or install the buildx plugin." >&2
  exit 1
fi

BUILD_ARGS=(buildx build --platform "${PLATFORMS}")
[[ "${PUSH}" -eq 1 ]] && BUILD_ARGS+=(--push) || BUILD_ARGS+=(--load)

build() {
  local name="$1"
  local context="$2"
  local tag="${TAG_PREFIX}_${name}"
  echo "==> Building ${tag} for ${PLATFORMS}"
  docker "${BUILD_ARGS[@]}" \
    --build-arg "BUILD_DATE=$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    -t "${tag}" \
    "${context}"
}

build "pipeline" "${SCRIPT_DIR}/pipeline"
build "seqnado"  "${SCRIPT_DIR}/seqnado"
build "ml_cpu"   "${SCRIPT_DIR}/ml_cpu"

# GPU container: CUDA only supports linux/amd64
echo "==> Building ${TAG_PREFIX}_ml_gpu for linux/amd64 (CUDA, no ARM support)"
docker buildx build --platform linux/amd64 \
  $([[ "${PUSH}" -eq 1 ]] && echo "--push" || echo "--load") \
  -t "${TAG_PREFIX}_ml_gpu" \
  "${SCRIPT_DIR}/ml_gpu"

echo ""
echo "Done. Images built for ${PLATFORMS}."
if [[ "${PUSH}" -eq 0 ]]; then
  echo "Note: --load only works for single-platform builds locally."
  echo "To push multi-arch manifests, re-run with --push."
fi
