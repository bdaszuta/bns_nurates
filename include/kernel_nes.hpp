//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_nes.hpp
//  \brief contains kernels for inelastic neutrino scattering
//         on electrons and positrons
//
// Computation of inelastic neutrino-electron and neutrino-positron
// scattering using Eq. (43) from Mezzacappa & Bruenn, ApJ v.410, p.740 (1993)
// https://ui.adsabs.harvard.edu/abs/1993ApJ...410..740M/abstract

#ifndef BNS_NURATES_INCLUDE_KERNEL_NEPS_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_NEPS_HPP_

#include "bns_nurates.hpp"
#include "functions.hpp"
#include "constants.hpp"

// Numerical constants
//---------------------------------------------------------------------------------------------------------------------
constexpr BS_REAL kTaylorSeriesEpsilon = 1e-3;


// Physical constats
//---------------------------------------------------------------------------------------------------------------------

// Saves the expression that are needed for the unapproximated integral
KOKKOS_INLINE_FUNCTION
void ComputeFDIForInelastic(BS_REAL w, BS_REAL wp, BS_REAL eta,
                            BS_REAL* fdi_diff_w, BS_REAL* fdi_diff_abs)
{
    BS_REAL abs_val = fabs(w - wp);

    fdi_diff_w[0] = FDI_p1(eta - wp) - FDI_p1(eta - w);
    fdi_diff_w[1] = FDI_p2(eta - wp) - FDI_p2(eta - w);
    fdi_diff_w[2] = FDI_p3(eta - wp) - FDI_p3(eta - w);
    fdi_diff_w[3] = FDI_p4(eta - wp) - FDI_p4(eta - w);
    fdi_diff_w[4] = FDI_p5(eta - wp) - FDI_p5(eta - w);

    fdi_diff_abs[0] = FDI_p3(eta) - FDI_p3(eta - abs_val);
    fdi_diff_abs[1] = FDI_p4(eta) - FDI_p4(eta - abs_val);
    fdi_diff_abs[2] = FDI_p5(eta) - FDI_p5(eta - abs_val);
}

// Functions that calculates the kernel integral
//=========================================================================================================================================

// Not approximated out kernel integral
KOKKOS_INLINE_FUNCTION
BS_REAL MezzacappaIntOut(BS_REAL w, BS_REAL wp, BS_REAL x, BS_REAL y, int sign,
                         BS_REAL b1, BS_REAL b2, BS_REAL* fdi_diff_w,
                         BS_REAL* fdi_diff_abs)
{
    constexpr BS_REAL one_fifth = 0.2;
    constexpr BS_REAL two       = 2;
    constexpr BS_REAL three     = 3;
    constexpr BS_REAL six       = 6;

    return (((b1 + b2) * (sign * fdi_diff_abs[2] - fdi_diff_w[4]) *
                 one_fifth // All G5 terms

             - b1 * (w + wp) *
                   (

                       fdi_diff_w[3]

                       + two * ((w + wp) * fdi_diff_w[2] +
                                three * w * wp * fdi_diff_w[1])

                           ) // G4(eta - wp) + G3(eta - wp) term

             + sign * ((b1 * x - b2 * y) * fdi_diff_abs[1] +
                       two * (b1 * x * x + b2 * y * y) * fdi_diff_abs[0])) /
                (w * w * wp * wp)

            - six * b1 * fdi_diff_w[0]);
}

// Taylor expansion in the lowest energy of the function MezzacappaIntOut and
// MezzacappaIntIn
KOKKOS_INLINE_FUNCTION
BS_REAL MezzacappaIntOneEnergy(BS_REAL x, BS_REAL y, int sign, BS_REAL b1,
                               BS_REAL b2, const BS_REAL* fdis)
{
    constexpr BS_REAL one_fifth = 0.2;
    constexpr BS_REAL two       = 2;
    constexpr BS_REAL three     = 3;
    constexpr BS_REAL four      = 4;
    constexpr BS_REAL six       = 6;

    return -sign * y *
           (two * (b1 + b2) * fdis[0]

            + (b1 * (y + four * x) + three * b2 * y) * fdis[1]

            + (b1 * (three * y - four * x) + b2 * y) * fdis[2]

            + ((b1 + six * b2) * y * y * one_fifth + b1 * x * y +
               two * b1 * x * x) *
                  fdis[3]

            - ((six * b1 + b2) * y * y * one_fifth - three * b1 * x * y +
               two * b1 * x * x) *
                  fdis[4]) /
           (x * x);
}

// Taylor expansion in both energies of the function MezzacappaIntOut and
// MezzacappaIntIn
KOKKOS_INLINE_FUNCTION
BS_REAL MezzacappaIntTwoEnergies(BS_REAL w, BS_REAL wp, BS_REAL x, BS_REAL y,
                                 BS_REAL b1, BS_REAL b2, const BS_REAL* fdis)
{
    constexpr BS_REAL two        = 2;
    constexpr BS_REAL four       = 4;
    constexpr BS_REAL eleven     = 11;
    constexpr BS_REAL twenty     = 20;
    constexpr BS_REAL twentyfive = 25;
    constexpr BS_REAL thirthy    = 30;

    return y * (w - wp) *
           ((b1 + b2) *
                (four * fdis[2]

                 + fdis[0] *
                       (eleven * y * y - twentyfive * x * y + twenty * x * x) /
                       thirthy)

            + (b1 - b2) * (two * x - y) * fdis[1]) /
           (x * x);
}

// Functions for kernel calculation
//===================================================================================================

// Calculates and saves the neutrino electron scattering in and out kernel for
// every neutrino species
KOKKOS_INLINE_FUNCTION
MyKernelOutput NESKernels(InelasticScattKernelParams* kernel_params,
                          MyEOSParams* eos_params)
{
    const BS_REAL T                    = eos_params->temp;
    const BS_REAL w                    = kernel_params->omega / T;
    const BS_REAL wp                   = kernel_params->omega_prime / T;
    const BS_REAL x                    = fmax(w, wp);
    const BS_REAL y                    = fmin(w, wp);
    const BS_REAL eta_e                = eos_params->mu_e / T;
    const BS_REAL exp_factor           = NEPSExpFunc(wp - w);
    const BS_REAL exp_factor_exchanged = NEPSExpFunc(w - wp);

    MyKernelOutput output;

    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;
    constexpr BS_REAL six  = 6;

    if (y > eta_e * kTaylorSeriesEpsilon)
    {
        BS_REAL fdi_diff_abs[3], fdi_diff_w[5];

        const int sign = 2 * signbit(wp - w) - 1;

        ComputeFDIForInelastic(w, wp, eta_e, fdi_diff_w, fdi_diff_abs);

        output.abs[id_nue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BPlus, kBS_NEPS_BZero,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_anue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BZero, kBS_NEPS_BPlus,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_nux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BMinus, kBS_NEPS_BZero,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_anux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BZero, kBS_NEPS_BMinus,
                             fdi_diff_w, fdi_diff_abs);
    }
    else if (x > eta_e * kTaylorSeriesEpsilon)
    {
        const int sign = 2 * signbit(wp - w) - 1;

        const BS_REAL fdis[5] = {
            FDI_p2(eta_e - x) - FDI_p2(eta_e), FDI_p1(eta_e - x), FDI_p1(eta_e),
            FDI_0(eta_e - x) + y * FermiDistr(zero, one, eta_e - x) / six,
            FDI_0(eta_e) - y * FermiDistr(zero, one, eta_e) / six};

        output.abs[id_nue] = kBS_NEPS_Const * POW2(T) *
                             MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BPlus,
                                                    kBS_NEPS_BZero, fdis);
        output.abs[id_anue] = kBS_NEPS_Const * POW2(T) *
                              MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BZero,
                                                     kBS_NEPS_BPlus, fdis);
        output.abs[id_nux] = kBS_NEPS_Const * POW2(T) *
                             MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BMinus,
                                                    kBS_NEPS_BZero, fdis);
        output.abs[id_anux] = kBS_NEPS_Const * POW2(T) *
                              MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BZero,
                                                     kBS_NEPS_BMinus, fdis);
    }
    else
    {
        const BS_REAL fdis[3] = {FermiDistr(0., 1., eta_e), FDI_0(eta_e),
                                 FDI_p1(eta_e)};

        output.abs[id_nue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BPlus,
                                     kBS_NEPS_BZero, fdis);
        output.abs[id_anue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BZero,
                                     kBS_NEPS_BPlus, fdis);
        output.abs[id_nux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BMinus,
                                     kBS_NEPS_BZero, fdis);
        output.abs[id_anux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BZero,
                                     kBS_NEPS_BMinus, fdis);
    }

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        output.em[idx]  = output.abs[idx];
        output.abs[idx] = output.abs[idx] * exp_factor;
        output.em[idx]  = -output.em[idx] * exp_factor_exchanged;
    }

    return output;
}

// Calculates and saves the neutrino positron scattering in and out kernel for
// every neutrino species
KOKKOS_INLINE_FUNCTION
MyKernelOutput NPSKernels(InelasticScattKernelParams* kernel_params,
                          MyEOSParams* eos_params)
{
    const BS_REAL T                    = eos_params->temp;
    const BS_REAL w                    = kernel_params->omega / T;
    const BS_REAL wp                   = kernel_params->omega_prime / T;
    const BS_REAL x                    = fmax(w, wp);
    const BS_REAL y                    = fmin(w, wp);
    const BS_REAL eta_p                = -eos_params->mu_e / T;
    const BS_REAL exp_factor           = NEPSExpFunc(wp - w);
    const BS_REAL exp_factor_exchanged = NEPSExpFunc(w - wp);

    MyKernelOutput output;

    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;
    constexpr BS_REAL six  = 6;

    if (y > eta_p * kTaylorSeriesEpsilon)
    {
        BS_REAL fdi_diff_abs[3], fdi_diff_w[5];

        const int sign = 2 * signbit(wp - w) - 1;

        ComputeFDIForInelastic(w, wp, eta_p, fdi_diff_w, fdi_diff_abs);

        output.abs[id_nue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BZero, kBS_NEPS_BPlus,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_anue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BPlus, kBS_NEPS_BZero,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_nux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BZero, kBS_NEPS_BMinus,
                             fdi_diff_w, fdi_diff_abs);
        output.abs[id_anux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOut(w, wp, x, y, sign, kBS_NEPS_BMinus, kBS_NEPS_BZero,
                             fdi_diff_w, fdi_diff_abs);
    }
    else if (x > eta_p * kTaylorSeriesEpsilon)
    {
        const int sign = 2 * signbit(wp - w) - 1;

        const BS_REAL fdis[5] = {
            FDI_p2(eta_p - x) - FDI_p2(eta_p), FDI_p1(eta_p - x), FDI_p1(eta_p),
            FDI_0(eta_p - x) + y * FermiDistr(zero, one, eta_p - x) / six,
            FDI_0(eta_p) - y * FermiDistr(zero, one, eta_p) / six};

        output.abs[id_nue] = kBS_NEPS_Const * POW2(T) *
                             MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BZero,
                                                    kBS_NEPS_BPlus, fdis);
        output.abs[id_anue] = kBS_NEPS_Const * POW2(T) *
                              MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BPlus,
                                                     kBS_NEPS_BZero, fdis);
        output.abs[id_nux] = kBS_NEPS_Const * POW2(T) *
                             MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BZero,
                                                    kBS_NEPS_BMinus, fdis);
        output.abs[id_anux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntOneEnergy(x, y, sign, kBS_NEPS_BMinus, kBS_NEPS_BZero,
                                   fdis);
    }
    else
    {
        const BS_REAL fdis[3] = {FermiDistr(0., 1., eta_p), FDI_0(eta_p),
                                 FDI_p1(eta_p)};

        output.abs[id_nue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BZero,
                                     kBS_NEPS_BPlus, fdis);
        output.abs[id_anue] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BPlus,
                                     kBS_NEPS_BZero, fdis);
        output.abs[id_nux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BZero,
                                     kBS_NEPS_BMinus, fdis);
        output.abs[id_anux] =
            kBS_NEPS_Const * POW2(T) *
            MezzacappaIntTwoEnergies(w, wp, x, y, kBS_NEPS_BMinus,
                                     kBS_NEPS_BZero, fdis);
    }

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        output.em[idx]  = output.abs[idx];
        output.abs[idx] = output.abs[idx] * exp_factor;
        output.em[idx]  = -output.em[idx] * exp_factor_exchanged;
    }

    return output;
}

// Calculates the full in and out kernels
KOKKOS_INLINE_FUNCTION
MyKernelOutput InelasticScattKernels(InelasticScattKernelParams* kernel_params,
                                     MyEOSParams* eos_params)
{
    MyKernelOutput nes_kernel = NESKernels(kernel_params, eos_params);
    MyKernelOutput nps_kernel = NPSKernels(kernel_params, eos_params);

    MyKernelOutput tot_kernel = {0};

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        tot_kernel.em[idx]  = nes_kernel.em[idx] + nps_kernel.em[idx];
        tot_kernel.abs[idx] = nes_kernel.abs[idx] + nps_kernel.abs[idx];
    }

    return tot_kernel;
}

KOKKOS_INLINE_FUNCTION
void InelasticKernelsTable(const int n, BS_REAL* nu_array,
                           GreyOpacityParams* grey_pars, M1MatrixKokkos2D* out)
{
    MyKernelOutput inel_1, inel_2;

    InelasticScattKernelParams inelastic_pars =
        grey_pars->kernel_pars.inelastic_kernel_params;
    for (int i = 0; i < n; ++i)
    {

        for (int j = i; j < n; ++j)
        {

            // compute the pair kernels
            inelastic_pars.omega       = nu_array[i];
            inelastic_pars.omega_prime = nu_array[j];
            inel_1 =
                InelasticScattKernels(&inelastic_pars, &grey_pars->eos_pars);

            inelastic_pars.omega       = nu_array[j];
            inelastic_pars.omega_prime = nu_array[i];
            inel_2 =
                InelasticScattKernels(&inelastic_pars, &grey_pars->eos_pars);


            for (int idx = 0; idx < total_num_species; ++idx)
            {
                out->m1_mat_em[idx][i][j] = inel_1.em[idx];
                out->m1_mat_em[idx][j][i] = inel_2.em[idx];

                out->m1_mat_ab[idx][i][j] = inel_1.abs[idx];
                out->m1_mat_ab[idx][j][i] = inel_2.abs[idx];
            }
        }
    }

    return;
}

// ============================================================================
// Precomputed NEPS constants for FDI reuse across kernel evaluations
// ============================================================================

// Precomputed FDI values at a single eta, constant across all (w, wp) pairs
struct NEPSPrecomputedEta
{
    BS_REAL eta;             // The eta value itself (eta_e or eta_p)
    BS_REAL fdi_p3_eta;      // FDI_p3(eta)
    BS_REAL fdi_p4_eta;      // FDI_p4(eta)
    BS_REAL fdi_p5_eta;      // FDI_p5(eta)
    BS_REAL fdi_p2_eta;      // FDI_p2(eta)
    BS_REAL fdi_p1_eta;      // FDI_p1(eta)
    BS_REAL fdi_0_eta;       // FDI_0(eta)
    BS_REAL fermi_distr_eta; // FermiDistr(0, 1, eta)
};

struct NEPSPrecomputed
{
    NEPSPrecomputedEta nes; // For NES: eta = eta_e = mu_e / T
    NEPSPrecomputedEta nps; // For NPS: eta = eta_p = -mu_e / T
    BS_REAL T;              // temperature (NOT T_inv, to preserve omega/T)
};

KOKKOS_INLINE_FUNCTION
NEPSPrecomputedEta PrecomputeNEPSEta(BS_REAL eta)
{
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;

    NEPSPrecomputedEta pre;
    pre.eta             = eta;
    pre.fdi_p3_eta      = FDI_p3(eta);
    pre.fdi_p4_eta      = FDI_p4(eta);
    pre.fdi_p5_eta      = FDI_p5(eta);
    pre.fdi_p2_eta      = FDI_p2(eta);
    pre.fdi_p1_eta      = FDI_p1(eta);
    pre.fdi_0_eta       = FDI_0(eta);
    pre.fermi_distr_eta = FermiDistr(zero, one, eta);
    return pre;
}

KOKKOS_INLINE_FUNCTION
NEPSPrecomputed PrecomputeNEPSParams(MyEOSParams* eos_params)
{
    NEPSPrecomputed pre;
    pre.T = eos_params->temp;

    const BS_REAL eta_e = eos_params->mu_e / pre.T;
    const BS_REAL eta_p = -eta_e;

    pre.nes = PrecomputeNEPSEta(eta_e);
    pre.nps = PrecomputeNEPSEta(eta_p);

    return pre;
}

// Version of ComputeFDIForInelastic that reuses precomputed FDI(eta) values
// and uses combined FDI evaluation to share exp() across orders.
// Expression trees for the differences are identical to the original.
KOKKOS_INLINE_FUNCTION
void ComputeFDIForInelasticFast(BS_REAL w, BS_REAL wp, BS_REAL eta,
                                const NEPSPrecomputedEta* pre,
                                BS_REAL* fdi_diff_w, BS_REAL* fdi_diff_abs)
{
    BS_REAL abs_val = fabs(w - wp);

    // Evaluate FDI orders 1-5 at (eta-w) and (eta-wp) with shared exp()
    BS_REAL p1_w, p2_w, p3_w, p4_w, p5_w;
    BS_REAL p1_wp, p2_wp, p3_wp, p4_wp, p5_wp;
    FDI_p12345(eta - w, &p1_w, &p2_w, &p3_w, &p4_w, &p5_w);
    FDI_p12345(eta - wp, &p1_wp, &p2_wp, &p3_wp, &p4_wp, &p5_wp);

    fdi_diff_w[0] = p1_wp - p1_w;
    fdi_diff_w[1] = p2_wp - p2_w;
    fdi_diff_w[2] = p3_wp - p3_w;
    fdi_diff_w[3] = p4_wp - p4_w;
    fdi_diff_w[4] = p5_wp - p5_w;

    // Evaluate FDI orders 3-5 at (eta-|w-wp|) with shared exp()
    BS_REAL p3_abs, p4_abs, p5_abs;
    FDI_p345(eta - abs_val, &p3_abs, &p4_abs, &p5_abs);

    // Reuse precomputed FDI_pN(eta) instead of recomputing
    fdi_diff_abs[0] = pre->fdi_p3_eta - p3_abs;
    fdi_diff_abs[1] = pre->fdi_p4_eta - p4_abs;
    fdi_diff_abs[2] = pre->fdi_p5_eta - p5_abs;
}

// Generic fast kernel for a single scattering type (NES or NPS).
// Uses omega/T division (not omega*T_inv) to preserve bitwise exactness.
// b1_arr/b2_arr hold the (b1,b2) pair for each species.
KOKKOS_INLINE_FUNCTION
MyKernelOutput NEPSSingleTypeFast(BS_REAL omega, BS_REAL omega_prime,
                                  const NEPSPrecomputedEta* pre, BS_REAL T,
                                  const BS_REAL* b1_arr, const BS_REAL* b2_arr)
{
    const BS_REAL w                    = omega / T;
    const BS_REAL wp                   = omega_prime / T;
    const BS_REAL x                    = fmax(w, wp);
    const BS_REAL y                    = fmin(w, wp);
    const BS_REAL eta                  = pre->eta;
    const BS_REAL exp_factor           = NEPSExpFunc(wp - w);
    const BS_REAL exp_factor_exchanged = NEPSExpFunc(w - wp);

    MyKernelOutput output;

    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;
    constexpr BS_REAL six  = 6;

    if (y > eta * kTaylorSeriesEpsilon)
    {
        BS_REAL fdi_diff_abs[3], fdi_diff_w[5];

        const int sign = 2 * signbit(wp - w) - 1;

        ComputeFDIForInelasticFast(w, wp, eta, pre, fdi_diff_w, fdi_diff_abs);

        for (int idx = 0; idx < total_num_species; ++idx)
        {
            output.abs[idx] =
                kBS_NEPS_Const * POW2(T) *
                MezzacappaIntOut(w, wp, x, y, sign, b1_arr[idx], b2_arr[idx],
                                 fdi_diff_w, fdi_diff_abs);
        }
    }
    else if (x > eta * kTaylorSeriesEpsilon)
    {
        const int sign = 2 * signbit(wp - w) - 1;

        // Reuse precomputed FDI values at eta
        const BS_REAL fdis[5] = {
            FDI_p2(eta - x) - pre->fdi_p2_eta, FDI_p1(eta - x), pre->fdi_p1_eta,
            FDI_0(eta - x) + y * FermiDistr(zero, one, eta - x) / six,
            pre->fdi_0_eta - y * pre->fermi_distr_eta / six};

        for (int idx = 0; idx < total_num_species; ++idx)
        {
            output.abs[idx] = kBS_NEPS_Const * POW2(T) *
                              MezzacappaIntOneEnergy(x, y, sign, b1_arr[idx],
                                                     b2_arr[idx], fdis);
        }
    }
    else
    {
        // Reuse precomputed FDI values at eta
        const BS_REAL fdis[3] = {pre->fermi_distr_eta, pre->fdi_0_eta,
                                 pre->fdi_p1_eta};

        for (int idx = 0; idx < total_num_species; ++idx)
        {
            output.abs[idx] = kBS_NEPS_Const * POW2(T) *
                              MezzacappaIntTwoEnergies(w, wp, x, y, b1_arr[idx],
                                                       b2_arr[idx], fdis);
        }
    }

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        output.em[idx]  = output.abs[idx];
        output.abs[idx] = output.abs[idx] * exp_factor;
        output.em[idx]  = -output.em[idx] * exp_factor_exchanged;
    }

    return output;
}

// Fast combined NES+NPS kernel using precomputed constants
KOKKOS_INLINE_FUNCTION
MyKernelOutput InelasticScattKernelsFast(BS_REAL omega, BS_REAL omega_prime,
                                         const NEPSPrecomputed* pre)
{
    // NES b1/b2 per species: nue=(BPlus,BZero), anue=(BZero,BPlus),
    //                        nux=(BMinus,BZero), anux=(BZero,BMinus)
    const BS_REAL nes_b1[total_num_species] = {kBS_NEPS_BPlus, kBS_NEPS_BZero,
                                               kBS_NEPS_BMinus, kBS_NEPS_BZero};
    const BS_REAL nes_b2[total_num_species] = {kBS_NEPS_BZero, kBS_NEPS_BPlus,
                                               kBS_NEPS_BZero, kBS_NEPS_BMinus};

    // NPS b1/b2 per species (swapped relative to NES)
    const BS_REAL nps_b1[total_num_species] = {kBS_NEPS_BZero, kBS_NEPS_BPlus,
                                               kBS_NEPS_BZero, kBS_NEPS_BMinus};
    const BS_REAL nps_b2[total_num_species] = {kBS_NEPS_BPlus, kBS_NEPS_BZero,
                                               kBS_NEPS_BMinus, kBS_NEPS_BZero};

    MyKernelOutput nes_kernel = NEPSSingleTypeFast(
        omega, omega_prime, &pre->nes, pre->T, nes_b1, nes_b2);
    MyKernelOutput nps_kernel = NEPSSingleTypeFast(
        omega, omega_prime, &pre->nps, pre->T, nps_b1, nps_b2);

    MyKernelOutput tot_kernel = {0};

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        tot_kernel.em[idx]  = nes_kernel.em[idx] + nps_kernel.em[idx];
        tot_kernel.abs[idx] = nes_kernel.abs[idx] + nps_kernel.abs[idx];
    }

    return tot_kernel;
}

#endif // BNS_NURATES_INCLUDE_KERNEL_NEPS_HPP_
