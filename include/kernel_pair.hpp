//=================================================
// bns-nurates neutrino opacities code
// Copyright(C) XXX, licensed under the YYY License
// ================================================
//! \file  kernel_pair.hpp
//  \brief contains pair kernels and associated helper functions

#ifndef BNS_NURATES_INCLUDE_KERNEL_PAIR_HPP_
#define BNS_NURATES_INCLUDE_KERNEL_PAIR_HPP_

#include "bns_nurates.hpp"
#include "constants.hpp"
#include "functions.hpp"


KOKKOS_INLINE_FUNCTION
void PairPsi(const int l, const BS_REAL y, const BS_REAL z, const BS_REAL eta,
             BS_REAL* psi_out)
{
    constexpr BS_REAL zero             = 0;
    constexpr BS_REAL four             = 4;
    constexpr BS_REAL fifteen          = 15;
    constexpr BS_REAL twenty           = 20;
    constexpr BS_REAL forty            = 40;
    constexpr BS_REAL eighty           = 80;
    constexpr BS_REAL onehundredtwenty = 120;
    constexpr BS_REAL twohundred       = 200;

    BS_ASSERT(l == 0);

    BS_ASSERT(isfinite(y) and y >= zero);
    BS_ASSERT(isfinite(y) and z >= zero);
    BS_ASSERT(isfinite(eta));

    /* The check on eta is necessary because, if eta is very large, the electron
    phase space is full and the reaction is suppressed. The kernel should in
    principle go to 0 automatically, however in some cases numerical
    cancellation results in small but negative (and so unphysical) rates. This
    check fixes that. */
    /* TODO: The threshold of 200 is somewhat arbitrary, can we find some better
    motivated number? */
    if (eta - y - z > twohundred)
    {
        psi_out[0] = zero;
        psi_out[1] = zero;
    }
    else
    {
        const BS_REAL FDI_p1_emy  = FDI_p1(eta - y);
        const BS_REAL FDI_p1_emz  = FDI_p1(eta - z);
        const BS_REAL FDI_p1_epy  = FDI_p1(eta + y);
        const BS_REAL FDI_p1_epz  = FDI_p1(eta + z);
        const BS_REAL FDI_p2_emy  = FDI_p2(eta - y);
        const BS_REAL FDI_p2_emz  = FDI_p2(eta - z);
        const BS_REAL FDI_p2_epy  = FDI_p2(eta + y);
        const BS_REAL FDI_p2_epz  = FDI_p2(eta + z);
        const BS_REAL FDI_p3_e    = FDI_p3(eta);
        const BS_REAL FDI_p3_emy  = FDI_p3(eta - y);
        const BS_REAL FDI_p3_emz  = FDI_p3(eta - z);
        const BS_REAL FDI_p3_epy  = FDI_p3(eta + y);
        const BS_REAL FDI_p3_epz  = FDI_p3(eta + z);
        const BS_REAL FDI_p3_emyz = FDI_p3(eta - y - z);
        const BS_REAL FDI_p3_epyz = FDI_p3(eta + y + z);
        const BS_REAL FDI_p4_e    = FDI_p4(eta);
        const BS_REAL FDI_p4_emy  = FDI_p4(eta - y);
        const BS_REAL FDI_p4_emz  = FDI_p4(eta - z);
        const BS_REAL FDI_p4_epy  = FDI_p4(eta + y);
        const BS_REAL FDI_p4_epz  = FDI_p4(eta + z);
        const BS_REAL FDI_p4_emyz = FDI_p4(eta - y - z);
        const BS_REAL FDI_p4_epyz = FDI_p4(eta + y + z);
        const BS_REAL FDI_p5_emy  = FDI_p5(eta - y);
        const BS_REAL FDI_p5_emz  = FDI_p5(eta - z);
        const BS_REAL FDI_p5_epy  = FDI_p5(eta + y);
        const BS_REAL FDI_p5_epz  = FDI_p5(eta + z);
        const BS_REAL FDI_p5_emyz = FDI_p5(eta - y - z);
        const BS_REAL FDI_p5_epyz = FDI_p5(eta + y + z);

        const BS_REAL x0 = twenty * FDI_p4_emy;
        const BS_REAL x1 = twenty * FDI_p4_epz;
        const BS_REAL x2 =
            onehundredtwenty * FDI_p2_emy - onehundredtwenty * FDI_p2_epz;
        const BS_REAL x3  = -forty * FDI_p3_emy + forty * FDI_p3_epz;
        const BS_REAL x4  = forty * FDI_p3_e;
        const BS_REAL x5  = forty * FDI_p3_emyz - x4;
        const BS_REAL x6  = -twenty * FDI_p4_e;
        const BS_REAL x7  = twenty * FDI_p4_emyz + x6;
        const BS_REAL x8  = -forty * FDI_p3_epyz + x4;
        const BS_REAL x9  = twenty * FDI_p4_epyz + x6;
        const BS_REAL x10 = -four * FDI_p5_emy + four * FDI_p5_emyz -
                            four * FDI_p5_emz + four * FDI_p5_epy -
                            four * FDI_p5_epyz + four * FDI_p5_epz;
        const BS_REAL x11 = twenty * FDI_p4_epy;
        const BS_REAL x12 = twenty * FDI_p4_emz;
        const BS_REAL x13 = -forty * FDI_p3_emz + forty * FDI_p3_epy;
        const BS_REAL x14 =
            onehundredtwenty * FDI_p2_emz - onehundredtwenty * FDI_p2_epy;

        const BS_REAL aux = fifteen * POW2(y * z);

        psi_out[0] =
            x10 +
            y * (-x0 + x1 + x7 +
                 y * (x3 + x5 +
                      z * (x2 + z * (-onehundredtwenty * FDI_p1_emy +
                                     onehundredtwenty * FDI_p1_epz))) +
                 z * (eighty * FDI_p3_emy - eighty * FDI_p3_epz - x2 * z)) +
            z * (x0 - x1 + x9 + z * (x3 + x8));

        psi_out[1] =
            x10 + y * (-x11 + x12 + x9 + y * (x13 + x8)) +
            z * (x11 - x12 + x7 +
                 y * (eighty * FDI_p3_emz - eighty * FDI_p3_epy - x14 * y) +
                 z * (x13 + x5 +
                      y * (x14 + y * (-onehundredtwenty * FDI_p1_emz +
                                      onehundredtwenty * FDI_p1_epy))));

        psi_out[0] /= aux;
        psi_out[1] /= aux;
    }
}

/* Calculate Phi_l(y,z) from Eqn. (10) of Pons et. al. (1998)
 *
 * Inputs:
 *      l:            mode number
 *      omega:        neutrino energy [MeV]
 *      omega_prime:  anti-neutrino energy [MeV]
 *      temp:         temperature [MeV]
 *      e_x:          neutrino species type (0: elentron, 1: mu/tau)
 *
 * Output:
 *      Phi_l(y,z) = (G^2 temp^2)/(pi (1 - e^{y+z})) [alpha1 Psi_l(y,z) + alpha2
 * Psi_l(z,y)]
 */
KOKKOS_INLINE_FUNCTION
void PairPhi(const BS_REAL omega, const BS_REAL omega_prime, const int l,
             const BS_REAL eta, const BS_REAL T, BS_REAL* phi_out)
{
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;

    const BS_REAL y = omega / T;
    const BS_REAL z = omega_prime / T;

    const BS_REAL aux = kBS_Pair_Phi * POW2(T) / (one - SafeExp(y + z));

    BS_REAL pair_psi[2] = {zero};

    PairPsi(l, y, z, eta, pair_psi);

    phi_out[0] = (POW2(kBS_Pair_Alpha1_0) * pair_psi[0] +
                  POW2(kBS_Pair_Alpha2_0) * pair_psi[1]) *
                 aux;
    phi_out[1] = (POW2(kBS_Pair_Alpha1_0) * pair_psi[1] +
                  POW2(kBS_Pair_Alpha2_0) * pair_psi[0]) *
                 aux;
    phi_out[2] = (POW2(kBS_Pair_Alpha1_1) * pair_psi[0] +
                  POW2(kBS_Pair_Alpha2_1) * pair_psi[1]) *
                 aux;
    phi_out[3] = (POW2(kBS_Pair_Alpha1_1) * pair_psi[1] +
                  POW2(kBS_Pair_Alpha2_1) * pair_psi[0]) *
                 aux;
}

KOKKOS_INLINE_FUNCTION
MyKernelOutput PairKernels(const MyEOSParams* eos_pars,
                           const PairKernelParams* kernel_pars)
{
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL half = 0.5;

    // EOS specific parameters
    const BS_REAL T   = eos_pars->temp;
    const BS_REAL eta = eos_pars->mu_e / T;

    // kernel specific parameters
    const BS_REAL omega       = kernel_pars->omega;
    const BS_REAL omega_prime = kernel_pars->omega_prime;

    BS_REAL pair_phi[4] = {zero};

    PairPhi(omega, omega_prime, 0, eta, T, pair_phi);

    MyKernelOutput pair_kernel;

    pair_kernel.em[id_nue]  = half * pair_phi[0];
    pair_kernel.em[id_anue] = half * pair_phi[1];
    pair_kernel.em[id_nux]  = half * pair_phi[2];
    pair_kernel.em[id_anux] = half * pair_phi[3];

    for (int idx = 0; idx < total_num_species; ++idx)
    {
        pair_kernel.abs[idx] =
            SafeExp((omega + omega_prime) / T) * pair_kernel.em[idx];
    }

    return pair_kernel;
}

KOKKOS_INLINE_FUNCTION
void PairKernels(const MyEOSParams* eos_pars,
                 const PairKernelParams* kernel_pars, MyKernelOutput* out_for,
                 MyKernelOutput* out_inv)
{
    *out_for = PairKernels(eos_pars, kernel_pars);

    out_inv->em[id_nue]  = out_for->em[id_anue];
    out_inv->em[id_anue] = out_for->em[id_nue];
    out_inv->em[id_nux]  = out_for->em[id_anux];
    out_inv->em[id_anux] = out_for->em[id_nux];

    out_inv->abs[id_nue]  = out_for->abs[id_anue];
    out_inv->abs[id_anue] = out_for->abs[id_nue];
    out_inv->abs[id_nux]  = out_for->abs[id_anux];
    out_inv->abs[id_anux] = out_for->abs[id_nux];
}

// ============================================================
// Precomputed FDI values for fast Pair kernel evaluation.
// PairPsi evaluates 28 FDI calls at arguments involving
//   eta, eta+-y, eta+-z, eta+-(y+z)
// where y = omega/T, z = omega_prime/T, eta = mu_e/T.
//
// For a given energy point, we can precompute all 10 FDI values
// at eta+x and eta-x (for FDI_p1..p5). The eta-only values
// (FDI_p3(eta), FDI_p4(eta)) are constant per thermodynamic point.
// This leaves only 6 FDI calls per (y,z) pair for eta+-(y+z).
// ============================================================

// FDI values at eta+x and eta-x for a single dimensionless energy x
struct PairFDIAtPoint
{
    BS_REAL fdi_p1_minus; // FDI_p1(eta - x)
    BS_REAL fdi_p1_plus;  // FDI_p1(eta + x)
    BS_REAL fdi_p2_minus; // FDI_p2(eta - x)
    BS_REAL fdi_p2_plus;  // FDI_p2(eta + x)
    BS_REAL fdi_p3_minus; // FDI_p3(eta - x)
    BS_REAL fdi_p3_plus;  // FDI_p3(eta + x)
    BS_REAL fdi_p4_minus; // FDI_p4(eta - x)
    BS_REAL fdi_p4_plus;  // FDI_p4(eta + x)
    BS_REAL fdi_p5_minus; // FDI_p5(eta - x)
    BS_REAL fdi_p5_plus;  // FDI_p5(eta + x)
};

// Precomputed data for Pair kernels: constant per thermodynamic point
struct PairPrecomputed
{
    BS_REAL eta;        // mu_e / T
    BS_REAL T;          // temperature
    BS_REAL FDI_p3_eta; // FDI_p3(eta)
    BS_REAL FDI_p4_eta; // FDI_p4(eta)
};

// Compute FDI values at eta+-x for a single dimensionless energy x
// Uses standalone FDI calls to preserve bitwise exactness
KOKKOS_INLINE_FUNCTION
PairFDIAtPoint PrecomputePairFDIAtPoint(const BS_REAL eta, const BS_REAL x)
{
    PairFDIAtPoint pt;
    const BS_REAL em = eta - x;
    const BS_REAL ep = eta + x;
    pt.fdi_p1_minus  = FDI_p1(em);
    pt.fdi_p2_minus  = FDI_p2(em);
    pt.fdi_p3_minus  = FDI_p3(em);
    pt.fdi_p4_minus  = FDI_p4(em);
    pt.fdi_p5_minus  = FDI_p5(em);
    pt.fdi_p1_plus   = FDI_p1(ep);
    pt.fdi_p2_plus   = FDI_p2(ep);
    pt.fdi_p3_plus   = FDI_p3(ep);
    pt.fdi_p4_plus   = FDI_p4(ep);
    pt.fdi_p5_plus   = FDI_p5(ep);
    return pt;
}

// Precompute the eta-only FDI values
KOKKOS_INLINE_FUNCTION
PairPrecomputed PrecomputePairParams(const BS_REAL eta, const BS_REAL T)
{
    PairPrecomputed pre;
    pre.eta        = eta;
    pre.T          = T;
    pre.FDI_p3_eta = FDI_p3(eta);
    pre.FDI_p4_eta = FDI_p4(eta);
    return pre;
}

// Fast PairPsi using precomputed FDI values for y and z.
// Only evaluates 6 FDI calls for eta+-(y+z).
KOKKOS_INLINE_FUNCTION
void PairPsiFast(const BS_REAL y, const BS_REAL z, const PairPrecomputed& pre,
                 const PairFDIAtPoint& fdi_y, const PairFDIAtPoint& fdi_z,
                 BS_REAL* psi_out)
{
    constexpr BS_REAL zero             = 0;
    constexpr BS_REAL four             = 4;
    constexpr BS_REAL fifteen          = 15;
    constexpr BS_REAL twenty           = 20;
    constexpr BS_REAL forty            = 40;
    constexpr BS_REAL eighty           = 80;
    constexpr BS_REAL onehundredtwenty = 120;
    constexpr BS_REAL twohundred       = 200;

    if (pre.eta - y - z > twohundred)
    {
        psi_out[0] = zero;
        psi_out[1] = zero;
    }
    else
    {
        // y-point FDI values
        const BS_REAL FDI_p1_emy = fdi_y.fdi_p1_minus;
        const BS_REAL FDI_p1_epy = fdi_y.fdi_p1_plus;
        const BS_REAL FDI_p2_emy = fdi_y.fdi_p2_minus;
        const BS_REAL FDI_p2_epy = fdi_y.fdi_p2_plus;
        const BS_REAL FDI_p3_emy = fdi_y.fdi_p3_minus;
        const BS_REAL FDI_p3_epy = fdi_y.fdi_p3_plus;
        const BS_REAL FDI_p4_emy = fdi_y.fdi_p4_minus;
        const BS_REAL FDI_p4_epy = fdi_y.fdi_p4_plus;
        const BS_REAL FDI_p5_emy = fdi_y.fdi_p5_minus;
        const BS_REAL FDI_p5_epy = fdi_y.fdi_p5_plus;

        // z-point FDI values
        const BS_REAL FDI_p1_emz = fdi_z.fdi_p1_minus;
        const BS_REAL FDI_p1_epz = fdi_z.fdi_p1_plus;
        const BS_REAL FDI_p2_emz = fdi_z.fdi_p2_minus;
        const BS_REAL FDI_p2_epz = fdi_z.fdi_p2_plus;
        const BS_REAL FDI_p3_emz = fdi_z.fdi_p3_minus;
        const BS_REAL FDI_p3_epz = fdi_z.fdi_p3_plus;
        const BS_REAL FDI_p4_emz = fdi_z.fdi_p4_minus;
        const BS_REAL FDI_p4_epz = fdi_z.fdi_p4_plus;
        const BS_REAL FDI_p5_emz = fdi_z.fdi_p5_minus;
        const BS_REAL FDI_p5_epz = fdi_z.fdi_p5_plus;

        // eta-only values (precomputed)
        const BS_REAL FDI_p3_e = pre.FDI_p3_eta;
        const BS_REAL FDI_p4_e = pre.FDI_p4_eta;

        // Only these 6 FDI calls remain per (y,z) pair
        const BS_REAL yz_sum      = y + z;
        const BS_REAL FDI_p3_emyz = FDI_p3(pre.eta - yz_sum);
        const BS_REAL FDI_p4_emyz = FDI_p4(pre.eta - yz_sum);
        const BS_REAL FDI_p5_emyz = FDI_p5(pre.eta - yz_sum);
        const BS_REAL FDI_p3_epyz = FDI_p3(pre.eta + yz_sum);
        const BS_REAL FDI_p4_epyz = FDI_p4(pre.eta + yz_sum);
        const BS_REAL FDI_p5_epyz = FDI_p5(pre.eta + yz_sum);

        // Exact same arithmetic as original PairPsi
        const BS_REAL x0 = twenty * FDI_p4_emy;
        const BS_REAL x1 = twenty * FDI_p4_epz;
        const BS_REAL x2 =
            onehundredtwenty * FDI_p2_emy - onehundredtwenty * FDI_p2_epz;
        const BS_REAL x3  = -forty * FDI_p3_emy + forty * FDI_p3_epz;
        const BS_REAL x4  = forty * FDI_p3_e;
        const BS_REAL x5  = forty * FDI_p3_emyz - x4;
        const BS_REAL x6  = -twenty * FDI_p4_e;
        const BS_REAL x7  = twenty * FDI_p4_emyz + x6;
        const BS_REAL x8  = -forty * FDI_p3_epyz + x4;
        const BS_REAL x9  = twenty * FDI_p4_epyz + x6;
        const BS_REAL x10 = -four * FDI_p5_emy + four * FDI_p5_emyz -
                            four * FDI_p5_emz + four * FDI_p5_epy -
                            four * FDI_p5_epyz + four * FDI_p5_epz;
        const BS_REAL x11 = twenty * FDI_p4_epy;
        const BS_REAL x12 = twenty * FDI_p4_emz;
        const BS_REAL x13 = -forty * FDI_p3_emz + forty * FDI_p3_epy;
        const BS_REAL x14 =
            onehundredtwenty * FDI_p2_emz - onehundredtwenty * FDI_p2_epy;

        const BS_REAL aux = fifteen * POW2(y * z);

        psi_out[0] =
            x10 +
            y * (-x0 + x1 + x7 +
                 y * (x3 + x5 +
                      z * (x2 + z * (-onehundredtwenty * FDI_p1_emy +
                                     onehundredtwenty * FDI_p1_epz))) +
                 z * (eighty * FDI_p3_emy - eighty * FDI_p3_epz - x2 * z)) +
            z * (x0 - x1 + x9 + z * (x3 + x8));

        psi_out[1] =
            x10 + y * (-x11 + x12 + x9 + y * (x13 + x8)) +
            z * (x11 - x12 + x7 +
                 y * (eighty * FDI_p3_emz - eighty * FDI_p3_epy - x14 * y) +
                 z * (x13 + x5 +
                      y * (x14 + y * (-onehundredtwenty * FDI_p1_emz +
                                      onehundredtwenty * FDI_p1_epy))));

        psi_out[0] /= aux;
        psi_out[1] /= aux;
    }
}

// Fast PairKernels using precomputed FDI values.
KOKKOS_INLINE_FUNCTION
MyKernelOutput PairKernelsFast(const BS_REAL omega, const BS_REAL omega_prime,
                               const PairPrecomputed& pre,
                               const PairFDIAtPoint& fdi_y,
                               const PairFDIAtPoint& fdi_z)
{
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL half = 0.5;
    constexpr BS_REAL one  = 1;

    const BS_REAL T = pre.T;
    const BS_REAL y = omega / T;
    const BS_REAL z = omega_prime / T;

    const BS_REAL aux = kBS_Pair_Phi * POW2(T) / (one - SafeExp(y + z));

    BS_REAL pair_psi[2] = {zero};
    PairPsiFast(y, z, pre, fdi_y, fdi_z, pair_psi);

    BS_REAL pair_phi[4];
    pair_phi[0] = (POW2(kBS_Pair_Alpha1_0) * pair_psi[0] +
                   POW2(kBS_Pair_Alpha2_0) * pair_psi[1]) *
                  aux;
    pair_phi[1] = (POW2(kBS_Pair_Alpha1_0) * pair_psi[1] +
                   POW2(kBS_Pair_Alpha2_0) * pair_psi[0]) *
                  aux;
    pair_phi[2] = (POW2(kBS_Pair_Alpha1_1) * pair_psi[0] +
                   POW2(kBS_Pair_Alpha2_1) * pair_psi[1]) *
                  aux;
    pair_phi[3] = (POW2(kBS_Pair_Alpha1_1) * pair_psi[1] +
                   POW2(kBS_Pair_Alpha2_1) * pair_psi[0]) *
                  aux;

    MyKernelOutput pair_kernel;

    pair_kernel.em[id_nue]  = half * pair_phi[0];
    pair_kernel.em[id_anue] = half * pair_phi[1];
    pair_kernel.em[id_nux]  = half * pair_phi[2];
    pair_kernel.em[id_anux] = half * pair_phi[3];

    const BS_REAL exp_factor = SafeExp((omega + omega_prime) / T);
    for (int idx = 0; idx < total_num_species; ++idx)
    {
        pair_kernel.abs[idx] = exp_factor * pair_kernel.em[idx];
    }

    return pair_kernel;
}

// Fast PairKernels with both forward and inverse output.
KOKKOS_INLINE_FUNCTION
void PairKernelsFast(const BS_REAL omega, const BS_REAL omega_prime,
                     const PairPrecomputed& pre, const PairFDIAtPoint& fdi_y,
                     const PairFDIAtPoint& fdi_z, MyKernelOutput* out_for,
                     MyKernelOutput* out_inv)
{
    *out_for = PairKernelsFast(omega, omega_prime, pre, fdi_y, fdi_z);

    out_inv->em[id_nue]  = out_for->em[id_anue];
    out_inv->em[id_anue] = out_for->em[id_nue];
    out_inv->em[id_nux]  = out_for->em[id_anux];
    out_inv->em[id_anux] = out_for->em[id_nux];

    out_inv->abs[id_nue]  = out_for->abs[id_anue];
    out_inv->abs[id_anue] = out_for->abs[id_nue];
    out_inv->abs[id_nux]  = out_for->abs[id_anux];
    out_inv->abs[id_anux] = out_for->abs[id_nux];
}

KOKKOS_INLINE_FUNCTION
void PairKernelsTable(const int n, const BS_REAL* nu_array,
                      GreyOpacityParams* grey_pars, M1MatrixKokkos2D* out)
{
    constexpr BS_REAL zero = 0;
    constexpr BS_REAL one  = 1;

    MyKernelOutput pair_1, pair_2;

    grey_pars->kernel_pars.pair_kernel_params.cos_theta = one;
    grey_pars->kernel_pars.pair_kernel_params.filter    = zero;
    grey_pars->kernel_pars.pair_kernel_params.lmax      = zero;
    grey_pars->kernel_pars.pair_kernel_params.mu        = one;
    grey_pars->kernel_pars.pair_kernel_params.mu_prime  = one;

    for (int i = 0; i < n; ++i)
    {
        grey_pars->kernel_pars.pair_kernel_params.omega = nu_array[i];

        for (int j = i; j < n; ++j)
        {
            grey_pars->kernel_pars.pair_kernel_params.omega_prime = nu_array[j];

            PairKernels(&grey_pars->eos_pars,
                        &grey_pars->kernel_pars.pair_kernel_params, &pair_1,
                        &pair_2);

            for (int idx = 0; idx < total_num_species; ++idx)
            {
                out->m1_mat_em[idx][i][j] = pair_1.em[idx];
                out->m1_mat_em[idx][j][i] = pair_2.em[idx];

                out->m1_mat_ab[idx][i][j] = pair_1.abs[idx];
                out->m1_mat_ab[idx][j][i] = pair_2.abs[idx];
            }
        }
    }

    return;
}

#endif // BNS_NURATES_INCLUDE_KERNEL_PAIR_HPP_
