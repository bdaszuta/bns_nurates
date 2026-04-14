#ifndef BNS_NURATES_PROFILE_H_
#define BNS_NURATES_PROFILE_H_

// Profiling infrastructure for bns_nurates opacity library.
// Include this header only from profiling drivers.

#include <chrono>
#include <cstdio>

struct Timer
{
    std::chrono::high_resolution_clock::time_point t_start;

    void tic()
    {
        t_start = std::chrono::high_resolution_clock::now();
    }

    double toc() const
    {
        auto t_end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double, std::micro>(t_end - t_start)
            .count();
    }
};

struct ProfileEntry
{
    const char* name;
    double total_us; // accumulated microseconds
    int indent;      // 0 = top-level, 1 = sub-reaction
};

constexpr int N_PROFILE_ENTRIES = 26;

enum ProfileIdx
{
    IDX_GAUSS_LEGENDRE = 0,
    IDX_NU_EQ_PARAMS,
    IDX_M1_DENSITIES_EQ,
    IDX_SPECTRAL_EQ,
    IDX_SPECTRAL_EQ_PAIR_BREM,
    IDX_SPECTRAL_EQ_NEPS,
    IDX_SPECTRAL_EQ_BETA,
    IDX_SPECTRAL_EQ_ISO,
    IDX_M1_EQ,
    IDX_M1_EQ_ISO,
    IDX_M1_EQ_BETA,
    IDX_M1_EQ_PAIR_BREM,
    IDX_M1_EQ_NEPS,
    IDX_M1_EQ_ASSEMBLY,
    IDX_DISTR_FROM_M1,
    IDX_SPECTRAL_M1,
    IDX_SPECTRAL_M1_PAIR_BREM,
    IDX_SPECTRAL_M1_NEPS,
    IDX_SPECTRAL_M1_BETA,
    IDX_SPECTRAL_M1_ISO,
    IDX_M1_M1,
    IDX_M1_M1_ISO,
    IDX_M1_M1_BETA,
    IDX_M1_M1_PAIR_BREM,
    IDX_M1_M1_NEPS,
    IDX_M1_M1_ASSEMBLY,
};

static ProfileEntry g_profile_entries[N_PROFILE_ENTRIES] = {
    {"GaussLegendre", 0.0, 0},
    {"NuEquilibriumParams", 0.0, 0},
    {"ComputeM1DensitiesEq", 0.0, 0},
    {"ComputeSpectralOpacities [eq]", 0.0, 0},
    {"  Pair+Brem integration", 0.0, 1},
    {"  NEPS integration", 0.0, 1},
    {"  Beta processes", 0.0, 1},
    {"  Iso scattering", 0.0, 1},
    {"ComputeM1Opacities [eq]", 0.0, 0},
    {"  Iso scattering (1D)", 0.0, 1},
    {"  Beta processes (1D)", 0.0, 1},
    {"  Pair+Brem (2D)", 0.0, 1},
    {"  NEPS / Inelastic scatt (2D)", 0.0, 1},
    {"  Assembly", 0.0, 1},
    {"CalculateDistrParamsFromM1", 0.0, 0},
    {"ComputeSpectralOpacities [m1]", 0.0, 0},
    {"  Pair+Brem integration", 0.0, 1},
    {"  NEPS integration", 0.0, 1},
    {"  Beta processes", 0.0, 1},
    {"  Iso scattering", 0.0, 1},
    {"ComputeM1Opacities [m1]", 0.0, 0},
    {"  Iso scattering (1D)", 0.0, 1},
    {"  Beta processes (1D)", 0.0, 1},
    {"  Pair+Brem (2D)", 0.0, 1},
    {"  NEPS / Inelastic scatt (2D)", 0.0, 1},
    {"  Assembly", 0.0, 1},
};

static void PrintProfilingSummary(int n_iterations)
{
    double total_us = 0.0;
    for (int i = 0; i < N_PROFILE_ENTRIES; ++i)
    {
        if (g_profile_entries[i].indent == 0)
            total_us += g_profile_entries[i].total_us;
    }

    printf("\n");
    printf("==================================================================="
           "==========\n");
    printf("  PROFILING SUMMARY  (N = %d iterations)\n", n_iterations);
    printf("==================================================================="
           "==========\n");
    printf("%-42s %12s %12s %8s\n", "Function", "Total (ms)", "Avg (us)",
           "% total");
    printf("-------------------------------------------------------------------"
           "----------\n");

    for (int i = 0; i < N_PROFILE_ENTRIES; ++i)
    {
        double total_ms = g_profile_entries[i].total_us / 1000.0;
        double avg_us   = g_profile_entries[i].total_us / n_iterations;
        double pct      = (total_us > 0) ?
                              (g_profile_entries[i].total_us / total_us * 100.0) :
                              0.0;

        if (g_profile_entries[i].indent == 0)
        {
            printf("%-42s %12.3f %12.3f %7.2f%%\n", g_profile_entries[i].name,
                   total_ms, avg_us, pct);
        }
        else
        {
            printf("  %-40s %12.3f %12.3f %7.2f%%\n", g_profile_entries[i].name,
                   total_ms, avg_us, pct);
        }
    }

    printf("-------------------------------------------------------------------"
           "----------\n");
    printf("%-42s %12.3f %12.3f %7.2f%%\n", "TOTAL", total_us / 1000.0,
           total_us / n_iterations, 100.0);
    printf("==================================================================="
           "==========\n");
}

#endif // BNS_NURATES_PROFILE_H_
