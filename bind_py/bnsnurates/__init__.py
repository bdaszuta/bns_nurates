import sys as _sys
import os as _os

_build_dir = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                            "..", "build")
if _build_dir not in _sys.path:
    _sys.path.insert(0, _build_dir)

from bns_core import *  # noqa: E402, F403

# Enum aliases at module level (SWIG-compatible)
kGauleg = Quadrature.kGauleg  # noqa: F405
kGaulag = Quadrature.kGaulag  # noqa: F405

# --- Dict-output wrappers (SWIG-compatible) ---

def _m1_to_dict(m1):
    return {
        "eta_0": list(m1.eta_0()),
        "kappa_0_a": list(m1.kappa_0_a()),
        "eta": list(m1.eta()),
        "kappa_a": list(m1.kappa_a()),
        "kappa_s": list(m1.kappa_s()),
    }

def _m1nt_to_dict(m1nt):
    d = _m1_to_dict(m1nt)
    d["eta_th"] = list(m1nt.eta_th())
    d["eta_non_th"] = list(m1nt.eta_non_th())
    d["kappa_a_th"] = list(m1nt.kappa_a_th())
    d["kappa_a_non_th"] = list(m1nt.kappa_a_non_th())
    return d

def _spectral_to_dict(spec):
    return {
        "j": list(spec.j()),
        "kappa": list(spec.kappa()),
        "j_s": list(spec.j_s()),
        "kappa_s": list(spec.kappa_s()),
    }

_ComputeM1Opacities = ComputeM1Opacities  # noqa: F405
_ComputeM1OpacitiesNonThermalSeparated = ComputeM1OpacitiesNonThermalSeparated  # noqa: F405
_ComputeSpectralOpacitiesNotStimulatedAbs = ComputeSpectralOpacitiesNotStimulatedAbs  # noqa: F405
_ComputeSpectralOpacitiesStimulatedAbs = ComputeSpectralOpacitiesStimulatedAbs  # noqa: F405

def ComputeM1Opacities(quad_1d, quad_2d, grey_params):
    return _m1_to_dict(_ComputeM1Opacities(quad_1d, quad_2d, grey_params))

def ComputeM1OpacitiesNonThermalSeparated(quad_1d, quad_2d, grey_params):
    return _m1nt_to_dict(_ComputeM1OpacitiesNonThermalSeparated(quad_1d, quad_2d, grey_params))

def ComputeSpectralOpacitiesNotStimulatedAbs(nu, quad_1d, grey_params):
    return _spectral_to_dict(_ComputeSpectralOpacitiesNotStimulatedAbs(nu, quad_1d, grey_params))

def ComputeSpectralOpacitiesStimulatedAbs(nu, quad_1d, grey_params):
    return _spectral_to_dict(_ComputeSpectralOpacitiesStimulatedAbs(nu, quad_1d, grey_params))


# --- Printing helpers (from SWIG %pythoncode) ---

def print_reactions(opacity_flags):
    # Handle both dict (SWIG compat) and OpacityFlags object
    if isinstance(opacity_flags, dict):
        f = opacity_flags
    else:
        f = {"use_abs_em": opacity_flags.use_abs_em,
             "use_pair": opacity_flags.use_pair,
             "use_brem": opacity_flags.use_brem,
             "use_inelastic_scatt": opacity_flags.use_inelastic_scatt,
             "use_iso": opacity_flags.use_iso}
    print('# Included reactions:')
    print('# Neutrino abs on nucleons - e+/e- capture : %d' % f["use_abs_em"])
    print('# Elastic scattering on nucleons           : %d' % f["use_iso"])
    print('# Nucleon-nucleon bremmstrahlung           : %d' % f["use_brem"])
    print('# e+ e- annihilation (and inverse)         : %d' % f["use_pair"])
    print('# Inelastic scattering on e+/e-            : %d' % f["use_inelastic_scatt"])
    print("")

def print_corrections(opacity_pars):
    if isinstance(opacity_pars, dict):
        p = opacity_pars
    else:
        p = {"use_dU": opacity_pars.use_dU,
             "use_dm_eff": opacity_pars.use_dm_eff,
             "use_WM_ab": opacity_pars.use_WM_ab,
             "use_WM_sc": opacity_pars.use_WM_sc,
             "use_decay": opacity_pars.use_decay,
             "use_BRT_brem": opacity_pars.use_BRT_brem,
             "neglect_blocking": opacity_pars.neglect_blocking,
             "use_NN_medium_corr": opacity_pars.use_NN_medium_corr}
    print('# Included corrections:')
    print('# dU correction                        : %d' % p["use_dU"])
    print('# Effective mass correction            : %d' % p["use_dm_eff"])
    print('# Weak magnetism on charged reactions  : %d' % p["use_WM_ab"])
    print('# Weak magnetism on neutral reactions  : %d' % p["use_WM_sc"])
    print('# Include (inverse) nucleon decays     : %d' % p["use_decay"])
    print('# Neglect blocking factors             : %d' % p["neglect_blocking"])
    print('# Bremsstrahlung as in Burrows+2006    : %d' % p["use_BRT_brem"])
    print('# Medium Fischer+16 correction to brem : %d' % p["use_NN_medium_corr"])
    print("")

def print_neutrino_quantities(m1_pars):
    # Handle both M1Quantities objects and dicts
    if hasattr(m1_pars, 'n'):
        n = m1_pars.n()
        J = m1_pars.J()
        chi = m1_pars.chi()
    else:
        n = m1_pars["n"]
        J = m1_pars["J"]
        chi = m1_pars["chi"]
    print("    nue            anue            nux            anux")
    print("n  %13.6e  %13.6e   %13.6e  %13.6e     (cm^-3)" %(n[0], n[1], n[2], n[3]))
    print("J  %13.6e  %13.6e   %13.6e  %13.6e (MeV cm^-3)" %(J[0], J[1], J[2], J[3]))
    print("chi %12.10f   %12.10f    %12.10f   %12.10f"   %(chi[0], chi[1], chi[2], chi[3]))
    print("")

def print_spectral_rates(spec_rates):
    print("      j             j_s           kappa         kappa_s")
    print(" nue %13.6e %13.6e %13.6e %13.6e" %(spec_rates['j'][0], spec_rates['j_s'][0], spec_rates['kappa'][0], spec_rates['kappa_s'][0]))
    print("anue %13.6e %13.6e %13.6e %13.6e" %(spec_rates['j'][1], spec_rates['j_s'][1], spec_rates['kappa'][1], spec_rates['kappa_s'][1]))
    print(" nux %13.6e %13.6e %13.6e %13.6e" %(spec_rates['j'][2], spec_rates['j_s'][2], spec_rates['kappa'][2], spec_rates['kappa_s'][2]))
    print("anux %13.6e %13.6e %13.6e %13.6e" %(spec_rates['j'][3], spec_rates['j_s'][3], spec_rates['kappa'][3], spec_rates['kappa_s'][3]))
    print("")

def print_integrated_rates(gray_rates):
    print("      eta0          eta1          kappa0        kappa1        scat1")
    print(" nue %13.6e %13.6e %13.6e %13.6e %13.6e" %(gray_rates['eta_0'][0], gray_rates['eta'][0], gray_rates['kappa_0_a'][0], gray_rates['kappa_a'][0], gray_rates['kappa_s'][0]))
    print("anue %13.6e %13.6e %13.6e %13.6e %13.6e" %(gray_rates['eta_0'][1], gray_rates['eta'][1], gray_rates['kappa_0_a'][1], gray_rates['kappa_a'][1], gray_rates['kappa_s'][1]))
    print(" nux %13.6e %13.6e %13.6e %13.6e %13.6e" %(gray_rates['eta_0'][2], gray_rates['eta'][2], gray_rates['kappa_0_a'][2], gray_rates['kappa_a'][2], gray_rates['kappa_s'][2]))
    print("anux %13.6e %13.6e %13.6e %13.6e %13.6e" %(gray_rates['eta_0'][3], gray_rates['eta'][3], gray_rates['kappa_0_a'][3], gray_rates['kappa_a'][3], gray_rates['kappa_s'][3]))
    print("")


__all__ = [
    "ComputeM1Opacities",
    "ComputeM1OpacitiesNonThermalSeparated",
    "ComputeSpectralOpacitiesNotStimulatedAbs",
    "ComputeSpectralOpacitiesStimulatedAbs",
    "print_reactions",
    "print_corrections",
    "print_neutrino_quantities",
    "print_spectral_rates",
    "print_integrated_rates",
]
