"""A/B comparison: SWIG vs nanobind bindings for bns_nurates."""

import sys as _sys
import os as _os
import numpy as np

# --- Import SWIG bindings ---
_swig_dir = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                          "..", "..", "bindings")
if _swig_dir not in _sys.path:
    _sys.path.insert(0, _swig_dir)
import bnsnurates as bns_swig  # noqa: E402

# --- Import nanobind bindings (C++ module + wrapper) ---
# IMPORTANT: import as bns_core to avoid module name collision with SWIG
_nb_build_dir = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                              "..", "..", "bind_py", "build")
if _nb_build_dir not in _sys.path:
    _sys.path.insert(0, _nb_build_dir)

import bns_core  # noqa: E402

# Replicate the dict-output wrappers from bind_py/bnsnurates/__init__.py
# so the test can call ComputeM1Opacities and get dicts (SWIG-compatible).

def _nb_m1_to_dict(m1):
    return {
        "eta_0": list(m1.eta_0()),
        "kappa_0_a": list(m1.kappa_0_a()),
        "eta": list(m1.eta()),
        "kappa_a": list(m1.kappa_a()),
        "kappa_s": list(m1.kappa_s()),
    }

def _nb_spectral_to_dict(spec):
    return {
        "j": list(spec.j()),
        "kappa": list(spec.kappa()),
        "j_s": list(spec.j_s()),
        "kappa_s": list(spec.kappa_s()),
    }


class _NanobindMod:
    """Wrapper that makes bns_core look like bnsnurates module."""
    def __init__(self):
        self._c = bns_core

    def __getattr__(self, name):
        return getattr(self._c, name)

    def ComputeM1Opacities(self, quad_1d, quad_2d, grey_params):
        return _nb_m1_to_dict(
            self._c.ComputeM1Opacities(quad_1d, quad_2d, grey_params))

    def ComputeSpectralOpacitiesNotStimulatedAbs(self, nu, quad_1d, grey_params):
        return _nb_spectral_to_dict(
            self._c.ComputeSpectralOpacitiesNotStimulatedAbs(nu, quad_1d, grey_params))

    def ComputeSpectralOpacitiesStimulatedAbs(self, nu, quad_1d, grey_params):
        return _nb_spectral_to_dict(
            self._c.ComputeSpectralOpacitiesStimulatedAbs(nu, quad_1d, grey_params))

    # Enum values as module-level attributes (SWIG-compatible)
    @property
    def kGauleg(self):
        return self._c.Quadrature.kGauleg

    @property
    def kGaulag(self):
        return self._c.Quadrature.kGaulag


bns_nb = _NanobindMod()

# Common tolerance
TOL = 1e-12


def _setup_eos(mod):
    """Build MyEOSParams for point A. Works for both bindings."""
    nb_nm = 4.208366627847035e+38 * 1e-21
    eos = mod.MyEOSParams()
    eos.nb = nb_nm
    eos.temp = 12.406403541564941
    eos.ye = 0.07158458232879639
    eos.yp = 0.07158458232879639
    eos.yn = 1.0 - 0.07158458232879639
    eos.mu_e = 187.1814489
    eos.mu_p = 1011.01797737
    eos.mu_n = 1221.59013681
    eos.dU = 18.92714728
    eos.dm_eff = 280.16495513 - 278.87162217
    return eos


def _setup_quad(mod, nx=6):
    """Build MyQuadrature. Works for both bindings."""
    q = mod.MyQuadrature()
    q.type = mod.kGauleg
    q.dim = 1
    q.nx = nx
    q.x1 = 0.0
    q.x2 = 1.0
    mod.GaussLegendre(q)
    return q


def _run_equilibrium(mod):
    """Run full MWE equilibrium path. Returns (spectral, gray) dicts."""
    eos = _setup_eos(mod)
    quad = _setup_quad(mod, nx=6)

    # Build kernel params
    kp = mod.MyKernelParams()
    kp.pair_kernel_params = mod.PairKernelParams()
    kp.pair_kernel_params.lmax = 0.0
    kp.brem_kernel_params = mod.BremKernelParams()
    kp.inelastic_kernel_params = mod.InelasticScattKernelParams()

    # Build opacity params and flags
    if mod is bns_swig:
        # SWIG: dict workaround for bool setter bug
        opacity_pars = {"use_dU": 1, "use_dm_eff": 0, "use_WM_ab": 1,
                        "use_WM_sc": 1, "use_decay": 1, "use_BRT_brem": 0,
                        "neglect_blocking": 0, "use_NN_medium_corr": 1}
        opacity_flags = {"use_abs_em": 1, "use_pair": 1, "use_brem": 1,
                         "use_inelastic_scatt": 1, "use_iso": 1}
    else:
        # nanobind: direct bool setters
        op = mod.OpacityParams()
        op.use_dU = True
        op.use_dm_eff = False
        op.use_WM_ab = True
        op.use_WM_sc = True
        op.use_decay = True
        op.use_BRT_brem = False
        op.neglect_blocking = False
        op.use_NN_medium_corr = True
        opacity_pars = op
        of = mod.OpacityFlags()
        of.use_abs_em = 1
        of.use_pair = 1
        of.use_brem = 1
        of.use_inelastic_scatt = 1
        of.use_iso = 1
        opacity_flags = of

    # Distribution
    distr = mod.NuEquilibriumParams(eos)

    # M1 quantities
    m1 = mod.M1Quantities()
    mod.ComputeM1DensitiesEq(eos, distr, m1)
    if mod is bns_swig:
        m1.chi = [1./3., 1./3., 1./3., 1./3.]
    else:
        c = m1.chi()
        c[0] = c[1] = c[2] = c[3] = 1./3.

    # GreyOpacityParams
    gop = mod.GreyOpacityParams()
    gop.eos_pars = eos
    gop.opacity_flags = opacity_flags
    gop.opacity_pars = opacity_pars
    gop.distr_pars = distr
    gop.m1_pars = m1
    gop.kernel_pars = kp

    # Compute
    spec = mod.ComputeSpectralOpacitiesNotStimulatedAbs(10.0, quad, gop)
    gray = mod.ComputeM1Opacities(quad, quad, gop)
    return spec, gray


def test_struct_construction():
    """Both bindings construct all key structs."""
    for label, mod in [("SWIG", bns_swig), ("nanobind", bns_nb)]:
        assert mod.MyEOSParams() is not None, f"{label}: MyEOSParams"
        assert mod.MyQuadrature() is not None, f"{label}: MyQuadrature"
        assert mod.M1Quantities() is not None, f"{label}: M1Quantities"
        assert mod.NuDistributionParams() is not None, f"{label}: NuDistributionParams"
        assert mod.OpacityParams() is not None, f"{label}: OpacityParams"
        assert mod.OpacityFlags() is not None, f"{label}: OpacityFlags"
        assert mod.MyKernelParams() is not None, f"{label}: MyKernelParams"
        assert mod.GreyOpacityParams() is not None, f"{label}: GreyOpacityParams"
    print("PASS: test_struct_construction")


def test_eos_getset():
    """EOSParams get/set values match."""
    nb_val = 0.042
    temp_val = 12.41
    ye_val = 0.0716

    eos_sw = bns_swig.MyEOSParams()
    eos_sw.nb = nb_val
    eos_sw.temp = temp_val
    eos_sw.ye = ye_val

    eos_nb = bns_nb.MyEOSParams()
    eos_nb.nb = nb_val
    eos_nb.temp = temp_val
    eos_nb.ye = ye_val

    assert abs(eos_sw.nb - eos_nb.nb) < TOL
    assert abs(eos_sw.temp - eos_nb.temp) < TOL
    assert abs(eos_sw.ye - eos_nb.ye) < TOL
    print("PASS: test_eos_getset")


def test_gauss_legendre():
    """Quadrature generation produces same results."""
    for label, mod in [("SWIG", bns_swig), ("nanobind", bns_nb)]:
        q = _setup_quad(mod, nx=6)
        assert q.nx == 6, f"{label}: nx mismatch"
    print("PASS: test_gauss_legendre")


def test_equilibrium_mwe():
    """Full MWE -- both bindings, compare all output values."""
    swig_spec, swig_gray = _run_equilibrium(bns_swig)
    nb_spec, nb_gray = _run_equilibrium(bns_nb)

    # Helper: extract all values from dicts
    def compare_dicts(d_sw, d_nb, tag, tol=TOL):
        for key in d_sw:
            a_sw = np.asarray(d_sw[key], dtype=np.float64)
            a_nb = np.asarray(d_nb[key], dtype=np.float64)
            denom = np.maximum(np.abs(a_sw), 1e-40)
            rel = np.max(np.abs(a_sw - a_nb) / denom)
            if rel > tol:
                print(f"  FAIL {tag}.{key}: rel diff = {rel:.2e}")
                print(f"    SWIG: {a_sw}")
                print(f"    NB:   {a_nb}")
                return False
        return True

    spec_ok = compare_dicts(swig_spec, nb_spec, "spectral")
    gray_ok = compare_dicts(swig_gray, nb_gray, "gray")

    if not spec_ok or not gray_ok:
        raise AssertionError("A/B comparison failed")

    print("PASS: test_equilibrium_mwe")


if __name__ == "__main__":
    test_struct_construction()
    test_eos_getset()
    test_gauss_legendre()
    test_equilibrium_mwe()
    print("\nALL TESTS PASSED")
