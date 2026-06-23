# bns_nurates nanobind Python Bindings

Python bindings for the bns_nurates neutrino reaction library using nanobind.

## Build

```bash
activate_conda   # conda activate dev
activate_spack   # source spack setup-env.sh
make bind
```

Output: `bind_py/build/bns_core.cpython-311-x86_64-linux-gnu.so`

## Import

```python
import sys
sys.path.insert(0, "bind_py")
import bnsnurates  # loads the nanobind wrapper
```

Or use directly (without dict wrappers):

```python
import sys
sys.path.insert(0, "bind_py/build")
import bns_core
```

## API Differences from SWIG

| Feature | SWIG | nanobind |
|---------|------|----------|
| Struct scalars | `eos.nb = value` | `eos.nb = value` (same) |
| C arrays (n[4], etc.) | `m1.n = [a,b,c,d]` (list) | `arr = m1.n()` returns numpy view; write via `arr[0] = val` |
| OpacityParams bools | Dict workaround: `{"use_dU": 1}` | Direct: `op.use_dU = True` |
| OpacityFlags ints | Dict workaround | Direct: `of.use_abs_em = 1` |
| Function outputs | Dict: `result["eta_0"]` | Dict (wrapper): same format |
| MyQuadrature setup | `q.type = kGauleg; q.nx = 6; q.x1=0; q.x2=1` | Same, PLUS `q.dim = 1` |
| Quadrature points/w | `q.points` (list) | `q.points()` (numpy view) |
| H[4][4] access | Broken in SWIG | `m1.H()` returns (4,4) numpy view |
| print_* helpers | Available | Available (copied from SWIG) |

## Test

```bash
make test-bindings
```

Or directly:

```bash
PYTHONPATH=bind_py/build:bind_py:bindings python3 tests/bind_py/test_bindings_nanobind.py
```

The test runs the full MWE (Chiesa+25 point A) with BOTH bindings and verifies
all output values match to 1e-12 relative tolerance.

## Known Issues

- `pair`, `brem`, and `inelastic_scatt` reaction channels segfault in the SWIG
  bindings. The A/B test uses only `use_abs_em` + `use_iso` which are safe.
- `BS_DEBUG` is enabled; uninitialized struct fields trigger assertions.
  Always set `quad.dim = 1` before calling `GaussLegendre`.
- CXXABI_1.3.15 incompatibility between conda and spack libstdc++ is handled
  by `tests/bind_py/run_all.sh` via LD_PRELOAD.
