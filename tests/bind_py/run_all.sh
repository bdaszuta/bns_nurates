#!/bin/bash
# Run all binding A/B tests sequentially.
#
# Environment:
#   activate_conda  - python3 + nanobind (sets $CONDA_PREFIX)
#   activate_spack  - gcc runtime (sets $SPACK_ENV)
#   Neither required if python3 and libgsl are on system paths.
set -e
cd "$(dirname "$0")"

# Python: prefer conda environment, fall back to PATH
if [ -n "$CONDA_PREFIX" ]; then
    PY="$CONDA_PREFIX/bin/python3"
else
    PY="$(which python3)"
fi

# libstdc++: if spack is active, preload its libstdc++ to satisfy
# CXXABI_1.3.15 which conda's bundled libstdc++.so.6 often lacks.
# NOTE: this is ONLY for CXXABI compatibility, NOT for GSL (GSL is
# resolved via RPATH baked into the .so at build time).
if [ -n "$SPACK_ENV" ]; then
    SPACK_LIBSTDCXX=""
    if [ -d "$SPACK_ENV/.spack-env/view/lib" ]; then
        SPACK_LIBSTDCXX="$SPACK_ENV/.spack-env/view/lib/libstdc++.so.6"
    fi
    if [ -z "$SPACK_LIBSTDCXX" ] || [ ! -f "$SPACK_LIBSTDCXX" ]; then
        SPACK_LIBSTDCXX=$(find /mnt/nimbus/_Installed/spack -name "libstdc++.so.6" \
            -path "*/gcc-*/lib*" 2>/dev/null | head -1)
    fi
    if [ -n "$SPACK_LIBSTDCXX" ] && [ -f "$SPACK_LIBSTDCXX" ]; then
        export LD_PRELOAD="$SPACK_LIBSTDCXX"
    fi
fi

tests=(
    test_bindings_nanobind.py
)

passed=0
failed=0

for t in "${tests[@]}"; do
    echo "=== $t ==="
    if $PY "$t"; then
        passed=$((passed + 1))
    else
        failed=$((failed + 1))
    fi
    echo
done

echo "=========================="
echo "$passed passed, $failed failed"
if [ $failed -gt 0 ]; then
    exit 1
fi
