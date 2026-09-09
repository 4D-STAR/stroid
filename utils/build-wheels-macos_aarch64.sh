#!/usr/bin/env bash
set -euo pipefail

if [[ $(uname -m) != "arm64" ]]; then
  echo "Error: This script is intended to run on an Apple Silicon (arm64) Mac."
  exit 1
fi

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <git-repo-url> [fourdst-wheels-dir]"
  echo "  fourdst-wheels-dir: optional local directory of fourdst wheels to"
  echo "  install from instead of PyPI (for bootstrapping a new stroid wheel)"
  exit 1
fi

for TOOL in pyenv cmake git; do
  if ! command -v "$TOOL" &> /dev/null; then
    echo "Error: ${TOOL} not found."
    exit 1
  fi
done

REPO_URL="$1"
LOCAL_FOURDST_WHEELS="${2:-}"
WORK_DIR="$(pwd)"
WHEEL_DIR="${WORK_DIR}/wheels_macos_aarch64_tmp"
FINAL_WHEEL_DIR="${WORK_DIR}/wheels_macos_aarch64"

echo "➤ Creating wheel output directories"
mkdir -p "${WHEEL_DIR}"
mkdir -p "${FINAL_WHEEL_DIR}"

export MACOSX_DEPLOYMENT_TARGET=15.0


TMPDIR="$(mktemp -d)"
echo "➤ Cloning ${REPO_URL} → ${TMPDIR}/project"
git clone --depth 1 "${REPO_URL}" "${TMPDIR}/project"
cd "${TMPDIR}/project"

FOURDST_PIN="$(grep -oE 'fourdst==[0-9][0-9a-zA-Z.]*' pyproject.toml | head -n1 || true)"
if [[ -n "${FOURDST_PIN}" ]]; then
  echo "➤ Project depends on ${FOURDST_PIN}; wheel repair will exclude fourdst libraries"
fi

PYTHON_VERSIONS=("3.9.23" "3.10.18" "3.11.13" "3.12.11" "3.13.5" "3.14.0rc1" "3.14.0rc1t")

eval "$(pyenv init -)"

for PY_VERSION in "${PYTHON_VERSIONS[@]}"; do
  (
    set -e

    pyenv shell "${PY_VERSION}"
    PY="$(pyenv which python)"

    echo "----------------------------------------------------------------"
    echo "➤ Building for $($PY --version) on macOS arm64"
    echo "----------------------------------------------------------------"

    "$PY" -m pip install --upgrade pip setuptools wheel meson-python delocate
    "$PY" -m pip install meson==1.9.1

    if [[ -n "${FOURDST_PIN}" ]]; then
      if [[ -n "${LOCAL_FOURDST_WHEELS}" ]]; then
        "$PY" -m pip install --force-reinstall \
          --find-links "${LOCAL_FOURDST_WHEELS}" "${FOURDST_PIN}"
      else
        "$PY" -m pip install --force-reinstall "${FOURDST_PIN}"
      fi
    fi

    echo "➤ Building wheel with ccache enabled"
    echo "➤ Found meson version $(meson --version)"

    CC="ccache clang" CXX="ccache clang++" \
      "$PY" -m pip wheel . --no-deps --no-build-isolation \
        -w "${WHEEL_DIR}" -v
    CURRENT_WHEEL=$(find "${WHEEL_DIR}" -name "*.whl" | head -n 1)

    echo "➤ Repairing wheel with delocate"
    if [[ -n "${FOURDST_PIN}" ]]; then
        delocate-wheel --require-archs arm64 \
          -e composition -e logging -e const -e reflect_cpp \
          -w "${FINAL_WHEEL_DIR}" -v "$CURRENT_WHEEL"
    else
        delocate-wheel --require-archs arm64 \
          -w "${FINAL_WHEEL_DIR}" -v "$CURRENT_WHEEL"
    fi

    REPAIRED_WHEEL="${FINAL_WHEEL_DIR}/$(basename "$CURRENT_WHEEL")"

    if [[ -n "${FOURDST_PIN}" ]]; then
      if unzip -l "${REPAIRED_WHEEL}" | grep -E 'libcomposition|liblogging|libconst|libreflect_cpp'; then
        echo "ERROR: repaired wheel contains vendored fourdst libraries"
        exit 1
      fi
    fi

    rm "$CURRENT_WHEEL"
  )
done

# Cleanup
rm -rf "${TMPDIR}"
rm -rf "${WHEEL_DIR}"

echo "All builds complete. Artifacts in ${FINAL_WHEEL_DIR}"