// ================================================
// nanobind bindings for bns_nurates
// ================================================

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include <memory>

#include "bns_nurates.hpp"
#include "m1_opacities.hpp"
#include "kernel_pair.hpp"
#include "kernel_brem.hpp"
#include "kernel_nes.hpp"

namespace nb = nanobind;
using namespace nb::literals;

using Arr4   = nb::ndarray<nb::numpy, BS_REAL, nb::shape<total_num_species>>;
using Arr20  = nb::ndarray<nb::numpy, BS_REAL, nb::shape<BS_N_MAX>>;
using Arr16  = nb::ndarray<nb::numpy, BS_REAL, nb::shape<16>>;
using ArrM1H = nb::ndarray<nb::numpy, BS_REAL, nb::shape<total_num_species, 4>>;

// ------------------------------------------------------------------
// M1OpacitiesNonThermalSeparated (defined here, not in headers)
// ------------------------------------------------------------------
struct M1OpacitiesNonThermalSeparated
{
    BS_REAL eta_0[total_num_species];
    BS_REAL kappa_0_a[total_num_species];
    BS_REAL eta[total_num_species];
    BS_REAL kappa_a[total_num_species];
    BS_REAL kappa_s[total_num_species];
    BS_REAL eta_th[total_num_species];
    BS_REAL eta_non_th[total_num_species];
    BS_REAL kappa_a_th[total_num_species];
    BS_REAL kappa_a_non_th[total_num_species];
};

inline M1OpacitiesNonThermalSeparated
ComputeM1OpacitiesNonThermalSeparated(MyQuadrature* quad_1d,
                                      MyQuadrature* quad_2d,
                                      GreyOpacityParams* grey_opacity_params)
{
    M1OpacitiesNonThermalSeparated out;
    M1Opacities base = ComputeM1OpacitiesGenericFormalism(
        quad_1d, quad_2d, grey_opacity_params, 1);
    for (int i = 0; i < total_num_species; ++i)
    {
        out.eta_0[i]          = base.eta_0[i];
        out.kappa_0_a[i]      = base.kappa_0_a[i];
        out.eta[i]            = base.eta[i];
        out.kappa_a[i]        = base.kappa_a[i];
        out.kappa_s[i]        = base.kappa_s[i];
        out.eta_th[i]         = base.eta[i];
        out.eta_non_th[i]     = 0.0;
        out.kappa_a_th[i]     = base.kappa_a[i];
        out.kappa_a_non_th[i] = 0.0;
    }
    return out;
}

// ------------------------------------------------------------------
// Shared reference helpers for array-view lifetime management
// ------------------------------------------------------------------

template <typename T>
struct shared_ref
{
    std::shared_ptr<T> data;
    int refs;
};

// ==================================================================
// Module definition
// ==================================================================

NB_MODULE(bns_core, m)
{

    // --------------------------------------------------------------
    // Enum
    // --------------------------------------------------------------
    nb::enum_<Quadrature>(m, "Quadrature")
        .value("kGauleg", kGauleg)
        .value("kGaulag", kGaulag);

    // --------------------------------------------------------------
    // Constants
    // --------------------------------------------------------------
    m.attr("total_num_species") = total_num_species;
    m.attr("BS_N_MAX")          = BS_N_MAX;
    m.attr("id_nue")            = id_nue;
    m.attr("id_anue")           = id_anue;
    m.attr("id_nux")            = id_nux;
    m.attr("id_anux")           = id_anux;

    // ==============================================================
    // STRUCT BINDINGS (bottom-up order)
    // ==============================================================

    // --- 1. BremKernelParams ---
    nb::class_<BremKernelParams>(m, "BremKernelParams")
        .def(nb::init<>())
        .def_rw("omega", &BremKernelParams::omega)
        .def_rw("omega_prime", &BremKernelParams::omega_prime)
        .def_rw("l", &BremKernelParams::l)
        .def_rw("use_NN_medium_corr", &BremKernelParams::use_NN_medium_corr);

    // --- 2. PairKernelParams ---
    nb::class_<PairKernelParams>(m, "PairKernelParams")
        .def(nb::init<>())
        .def_rw("omega", &PairKernelParams::omega)
        .def_rw("omega_prime", &PairKernelParams::omega_prime)
        .def_rw("cos_theta", &PairKernelParams::cos_theta)
        .def_rw("mu", &PairKernelParams::mu)
        .def_rw("mu_prime", &PairKernelParams::mu_prime)
        .def_rw("lmax", &PairKernelParams::lmax)
        .def_rw("filter", &PairKernelParams::filter);

    // --- 3. InelasticScattKernelParams ---
    nb::class_<InelasticScattKernelParams>(m, "InelasticScattKernelParams")
        .def(nb::init<>())
        .def_rw("omega", &InelasticScattKernelParams::omega)
        .def_rw("omega_prime", &InelasticScattKernelParams::omega_prime);

    // --- 4. MyKernelParams ---
    nb::class_<MyKernelParams>(m, "MyKernelParams")
        .def(nb::init<>())
        .def_rw("pair_kernel_params", &MyKernelParams::pair_kernel_params)
        .def_rw("brem_kernel_params", &MyKernelParams::brem_kernel_params)
        .def_rw("inelastic_kernel_params",
                &MyKernelParams::inelastic_kernel_params);

    // --- 5. MyKernelOutput: em[4] (def), abs[4] (def) ---
    nb::class_<MyKernelOutput>(m, "MyKernelOutput")
        .def(nb::init<>())
        .def("em",
             [](MyKernelOutput& obj) -> Arr4
             {
                 auto* ref = new shared_ref<MyKernelOutput>{
                     std::make_shared<MyKernelOutput>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<MyKernelOutput>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->em, {total_num_species},
                             std::move(owner));
             })
        .def("abs",
             [](MyKernelOutput& obj) -> Arr4
             {
                 auto* ref = new shared_ref<MyKernelOutput>{
                     std::make_shared<MyKernelOutput>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<MyKernelOutput>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->abs, {total_num_species},
                             std::move(owner));
             });

    // --- 6. MyEOSParams: all scalars ---
    nb::class_<MyEOSParams>(m, "MyEOSParams")
        .def(nb::init<>())
        .def_rw("nb", &MyEOSParams::nb)
        .def_rw("temp", &MyEOSParams::temp)
        .def_rw("ye", &MyEOSParams::ye)
        .def_rw("yp", &MyEOSParams::yp)
        .def_rw("yn", &MyEOSParams::yn)
        .def_rw("mu_p", &MyEOSParams::mu_p)
        .def_rw("mu_n", &MyEOSParams::mu_n)
        .def_rw("mu_e", &MyEOSParams::mu_e)
        .def_rw("mu_mu", &MyEOSParams::mu_mu)
        .def_rw("dU", &MyEOSParams::dU)
        .def_rw("dm_eff", &MyEOSParams::dm_eff);

    // --- 7. MyOpacity: abs[4] (def), em[4] (def) ---
    nb::class_<MyOpacity>(m, "MyOpacity")
        .def(nb::init<>())
        .def("abs",
             [](MyOpacity& obj) -> Arr4
             {
                 auto* ref = new shared_ref<MyOpacity>{
                     std::make_shared<MyOpacity>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<MyOpacity>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->abs, {total_num_species},
                             std::move(owner));
             })
        .def("em",
             [](MyOpacity& obj) -> Arr4
             {
                 auto* ref = new shared_ref<MyOpacity>{
                     std::make_shared<MyOpacity>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<MyOpacity>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->em, {total_num_species},
                             std::move(owner));
             });

    // --- 8. OpacityParams: all bool def_rw ---
    nb::class_<OpacityParams>(m, "OpacityParams")
        .def(nb::init<>())
        .def_rw("use_dU", &OpacityParams::use_dU)
        .def_rw("use_dm_eff", &OpacityParams::use_dm_eff)
        .def_rw("use_WM_ab", &OpacityParams::use_WM_ab)
        .def_rw("use_WM_sc", &OpacityParams::use_WM_sc)
        .def_rw("use_decay", &OpacityParams::use_decay)
        .def_rw("use_BRT_brem", &OpacityParams::use_BRT_brem)
        .def_rw("use_NN_medium_corr", &OpacityParams::use_NN_medium_corr)
        .def_rw("neglect_blocking", &OpacityParams::neglect_blocking);

    // --- 9. OpacityFlags: all int def_rw ---
    nb::class_<OpacityFlags>(m, "OpacityFlags")
        .def(nb::init<>())
        .def_rw("use_abs_em", &OpacityFlags::use_abs_em)
        .def_rw("use_pair", &OpacityFlags::use_pair)
        .def_rw("use_brem", &OpacityFlags::use_brem)
        .def_rw("use_inelastic_scatt", &OpacityFlags::use_inelastic_scatt)
        .def_rw("use_iso", &OpacityFlags::use_iso);

    // --- 10. NuDistributionParams: 7 arrays of size 4 (def) ---
    nb::class_<NuDistributionParams>(m, "NuDistributionParams")
        .def(nb::init<>())
        .def("w_t",
             [](NuDistributionParams& obj) -> Arr4
             {
                 auto* ref = new shared_ref<NuDistributionParams>{
                     std::make_shared<NuDistributionParams>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<NuDistributionParams>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->w_t, {total_num_species},
                             std::move(owner));
             })
        .def("temp_t",
             [](NuDistributionParams& obj) -> Arr4
             {
                 auto* ref = new shared_ref<NuDistributionParams>{
                     std::make_shared<NuDistributionParams>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<NuDistributionParams>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->temp_t, {total_num_species},
                             std::move(owner));
             })
        .def("eta_t",
             [](NuDistributionParams& obj) -> Arr4
             {
                 auto* ref = new shared_ref<NuDistributionParams>{
                     std::make_shared<NuDistributionParams>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<NuDistributionParams>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->eta_t, {total_num_species},
                             std::move(owner));
             })
        .def("w_f",
             [](NuDistributionParams& obj) -> Arr4
             {
                 auto* ref = new shared_ref<NuDistributionParams>{
                     std::make_shared<NuDistributionParams>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<NuDistributionParams>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->w_f, {total_num_species},
                             std::move(owner));
             })
        .def("temp_f",
             [](NuDistributionParams& obj) -> Arr4
             {
                 auto* ref = new shared_ref<NuDistributionParams>{
                     std::make_shared<NuDistributionParams>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<NuDistributionParams>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->temp_f, {total_num_species},
                             std::move(owner));
             })
        .def("c_f",
             [](NuDistributionParams& obj) -> Arr4
             {
                 auto* ref = new shared_ref<NuDistributionParams>{
                     std::make_shared<NuDistributionParams>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<NuDistributionParams>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->c_f, {total_num_species},
                             std::move(owner));
             })
        .def("beta_f",
             [](NuDistributionParams& obj) -> Arr4
             {
                 auto* ref = new shared_ref<NuDistributionParams>{
                     std::make_shared<NuDistributionParams>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<NuDistributionParams>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->beta_f, {total_num_species},
                             std::move(owner));
             });

    // --- 11. M1Quantities: n[4], J[4], H[4][4], chi[4] ---
    nb::class_<M1Quantities>(m, "M1Quantities")
        .def(nb::init<>())
        .def("n",
             [](M1Quantities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1Quantities>{
                     std::make_shared<M1Quantities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<M1Quantities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->n, {total_num_species},
                             std::move(owner));
             })
        .def("J",
             [](M1Quantities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1Quantities>{
                     std::make_shared<M1Quantities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<M1Quantities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->J, {total_num_species},
                             std::move(owner));
             })
        .def("H",
             [](M1Quantities& obj) -> ArrM1H
             {
                 auto* ref = new shared_ref<M1Quantities>{
                     std::make_shared<M1Quantities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<M1Quantities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return ArrM1H(&ref->data->H[0][0], {total_num_species, 4},
                               std::move(owner));
             })
        .def("chi",
             [](M1Quantities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1Quantities>{
                     std::make_shared<M1Quantities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<M1Quantities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->chi, {total_num_species},
                             std::move(owner));
                 })
                 .def("set_n",
                 [](M1Quantities& obj, Arr4 arr)
                 {
                  for (int i = 0; i < total_num_species; ++i)
                      obj.n[i] = arr(i);
                 })
                 .def("set_J",
                 [](M1Quantities& obj, Arr4 arr)
                 {
                  for (int i = 0; i < total_num_species; ++i)
                      obj.J[i] = arr(i);
                 })
                 .def("set_chi",
                 [](M1Quantities& obj, Arr4 arr)
                 {
                  for (int i = 0; i < total_num_species; ++i)
                      obj.chi[i] = arr(i);
                 });

    // --- 12. M1Opacities: eta_0[4], kappa_0_a[4], eta[4], kappa_a[4],
    // kappa_s[4] ---
    nb::class_<M1Opacities>(m, "M1Opacities")
        .def(nb::init<>())
        .def("eta_0",
             [](M1Opacities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1Opacities>{
                     std::make_shared<M1Opacities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<M1Opacities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->eta_0, {total_num_species},
                             std::move(owner));
             })
        .def("kappa_0_a",
             [](M1Opacities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1Opacities>{
                     std::make_shared<M1Opacities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<M1Opacities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->kappa_0_a, {total_num_species},
                             std::move(owner));
             })
        .def("eta",
             [](M1Opacities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1Opacities>{
                     std::make_shared<M1Opacities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<M1Opacities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->eta, {total_num_species},
                             std::move(owner));
             })
        .def("kappa_a",
             [](M1Opacities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1Opacities>{
                     std::make_shared<M1Opacities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<M1Opacities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->kappa_a, {total_num_species},
                             std::move(owner));
             })
        .def("kappa_s",
             [](M1Opacities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1Opacities>{
                     std::make_shared<M1Opacities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<M1Opacities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->kappa_s, {total_num_species},
                             std::move(owner));
             });

    // --- 13. SpectralOpacities: j[4], kappa[4], j_s[4], kappa_s[4] ---
    nb::class_<SpectralOpacities>(m, "SpectralOpacities")
        .def(nb::init<>())
        .def("j",
             [](SpectralOpacities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<SpectralOpacities>{
                     std::make_shared<SpectralOpacities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<SpectralOpacities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->j, {total_num_species},
                             std::move(owner));
             })
        .def("kappa",
             [](SpectralOpacities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<SpectralOpacities>{
                     std::make_shared<SpectralOpacities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<SpectralOpacities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->kappa, {total_num_species},
                             std::move(owner));
             })
        .def("j_s",
             [](SpectralOpacities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<SpectralOpacities>{
                     std::make_shared<SpectralOpacities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<SpectralOpacities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->j_s, {total_num_species},
                             std::move(owner));
             })
        .def("kappa_s",
             [](SpectralOpacities& obj) -> Arr4
             {
                 auto* ref = new shared_ref<SpectralOpacities>{
                     std::make_shared<SpectralOpacities>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<SpectralOpacities>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->kappa_s, {total_num_species},
                             std::move(owner));
             });

    // --- 14. GreyOpacityParams: all value members def_rw ---
    nb::class_<GreyOpacityParams>(m, "GreyOpacityParams")
        .def(nb::init<>())
        .def_rw("opacity_pars", &GreyOpacityParams::opacity_pars)
        .def_rw("kernel_pars", &GreyOpacityParams::kernel_pars)
        .def_rw("eos_pars", &GreyOpacityParams::eos_pars)
        .def_rw("distr_pars", &GreyOpacityParams::distr_pars)
        .def_rw("m1_pars", &GreyOpacityParams::m1_pars)
        .def_rw("opacity_flags", &GreyOpacityParams::opacity_flags);

    // --- 15. MyQuadrature: scalar def_rw, points/w def ---
    nb::class_<MyQuadrature>(m, "MyQuadrature")
        .def(nb::init<>())
        .def_rw("type", &MyQuadrature::type)
        .def_rw("alpha", &MyQuadrature::alpha)
        .def_rw("dim", &MyQuadrature::dim)
        .def_rw("nx", &MyQuadrature::nx)
        .def_rw("ny", &MyQuadrature::ny)
        .def_rw("nz", &MyQuadrature::nz)
        .def_rw("x1", &MyQuadrature::x1)
        .def_rw("x2", &MyQuadrature::x2)
        .def_rw("y1", &MyQuadrature::y1)
        .def_rw("y2", &MyQuadrature::y2)
        .def_rw("z1", &MyQuadrature::z1)
        .def_rw("z2", &MyQuadrature::z2)
        .def("points",
             [](MyQuadrature& obj) -> Arr20
             {
                 auto* ref = new shared_ref<MyQuadrature>{
                     std::make_shared<MyQuadrature>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<MyQuadrature>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr20(ref->data->points, {BS_N_MAX}, std::move(owner));
             })
        .def("w",
             [](MyQuadrature& obj) -> Arr20
             {
                 auto* ref = new shared_ref<MyQuadrature>{
                     std::make_shared<MyQuadrature>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<shared_ref<MyQuadrature>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr20(ref->data->w, {BS_N_MAX}, std::move(owner));
             });

    // --- 16. MyFunction: dim int def_rw (skip function pointers) ---
    nb::class_<MyFunction>(m, "MyFunction")
        .def(nb::init<>())
        .def_rw("dim", &MyFunction::dim);

    // --- 17. MyQuadratureIntegrand: n int def_rw, integrand[16] def ---
    nb::class_<MyQuadratureIntegrand>(m, "MyQuadratureIntegrand")
        .def(nb::init<>())
        .def_rw("n", &MyQuadratureIntegrand::n)
        .def("integrand",
             [](MyQuadratureIntegrand& obj) -> Arr16
             {
                 auto* ref = new shared_ref<MyQuadratureIntegrand>{
                     std::make_shared<MyQuadratureIntegrand>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<MyQuadratureIntegrand>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr16(ref->data->integrand, {16}, std::move(owner));
             });

    // --- 18. MyFunctionMultiD: dim int def_rw (skip function pointers) ---
    nb::class_<MyFunctionMultiD>(m, "MyFunctionMultiD")
        .def(nb::init<>())
        .def_rw("dim", &MyFunctionMultiD::dim);

    // ==============================================================
    // Internal structs needed by SpectralIntegrandParams
    // ==============================================================

    // --- PairFDIAtPoint ---
    nb::class_<PairFDIAtPoint>(m, "PairFDIAtPoint")
        .def(nb::init<>())
        .def_rw("fdi_p1_minus", &PairFDIAtPoint::fdi_p1_minus)
        .def_rw("fdi_p1_plus", &PairFDIAtPoint::fdi_p1_plus)
        .def_rw("fdi_p2_minus", &PairFDIAtPoint::fdi_p2_minus)
        .def_rw("fdi_p2_plus", &PairFDIAtPoint::fdi_p2_plus)
        .def_rw("fdi_p3_minus", &PairFDIAtPoint::fdi_p3_minus)
        .def_rw("fdi_p3_plus", &PairFDIAtPoint::fdi_p3_plus)
        .def_rw("fdi_p4_minus", &PairFDIAtPoint::fdi_p4_minus)
        .def_rw("fdi_p4_plus", &PairFDIAtPoint::fdi_p4_plus)
        .def_rw("fdi_p5_minus", &PairFDIAtPoint::fdi_p5_minus)
        .def_rw("fdi_p5_plus", &PairFDIAtPoint::fdi_p5_plus);

    // --- PairPrecomputed ---
    nb::class_<PairPrecomputed>(m, "PairPrecomputed")
        .def(nb::init<>())
        .def_rw("eta", &PairPrecomputed::eta)
        .def_rw("T", &PairPrecomputed::T)
        .def_rw("FDI_p3_eta", &PairPrecomputed::FDI_p3_eta)
        .def_rw("FDI_p4_eta", &PairPrecomputed::FDI_p4_eta);

    // --- NEPSPrecomputedEta ---
    nb::class_<NEPSPrecomputedEta>(m, "NEPSPrecomputedEta")
        .def(nb::init<>())
        .def_rw("eta", &NEPSPrecomputedEta::eta)
        .def_rw("fdi_p3_eta", &NEPSPrecomputedEta::fdi_p3_eta)
        .def_rw("fdi_p4_eta", &NEPSPrecomputedEta::fdi_p4_eta)
        .def_rw("fdi_p5_eta", &NEPSPrecomputedEta::fdi_p5_eta)
        .def_rw("fdi_p2_eta", &NEPSPrecomputedEta::fdi_p2_eta)
        .def_rw("fdi_p1_eta", &NEPSPrecomputedEta::fdi_p1_eta)
        .def_rw("fdi_0_eta", &NEPSPrecomputedEta::fdi_0_eta)
        .def_rw("fermi_distr_eta", &NEPSPrecomputedEta::fermi_distr_eta);

    // --- NEPSPrecomputed ---
    nb::class_<NEPSPrecomputed>(m, "NEPSPrecomputed")
        .def(nb::init<>())
        .def_rw("nes", &NEPSPrecomputed::nes)
        .def_rw("nps", &NEPSPrecomputed::nps)
        .def_rw("T", &NEPSPrecomputed::T);

    // --- BremPrecomputedChannel ---
    nb::class_<BremPrecomputedChannel>(m, "BremPrecomputedChannel")
        .def(nb::init<>())
        .def_rw("eta_star", &BremPrecomputedChannel::eta_star)
        .def_rw("y", &BremPrecomputedChannel::y)
        .def_rw("gamma", &BremPrecomputedChannel::gamma)
        .def_rw("half_gamma_gb", &BremPrecomputedChannel::half_gamma_gb)
        .def_rw("f_u", &BremPrecomputedChannel::f_u)
        .def_rw("pow_eta_neg5half", &BremPrecomputedChannel::pow_eta_neg5half)
        .def_rw("exp_neg_y_12", &BremPrecomputedChannel::exp_neg_y_12)
        .def_rw("pow_y_c2", &BremPrecomputedChannel::pow_y_c2)
        .def_rw("sqrt_eta", &BremPrecomputedChannel::sqrt_eta)
        .def_rw("eta_sq", &BremPrecomputedChannel::eta_sq)
        .def_rw("y4", &BremPrecomputedChannel::y4)
        .def_rw("y_sq", &BremPrecomputedChannel::y_sq)
        .def_rw("h_brem", &BremPrecomputedChannel::h_brem)
        .def_rw("p_brem", &BremPrecomputedChannel::p_brem);

    // --- BremPrecomputed ---
    nb::class_<BremPrecomputed>(m, "BremPrecomputed")
        .def(nb::init<>())
        .def_rw("ch_nn", &BremPrecomputed::ch_nn)
        .def_rw("ch_pp", &BremPrecomputed::ch_pp)
        .def_rw("ch_np", &BremPrecomputed::ch_np)
        .def_rw("nn", &BremPrecomputed::nn)
        .def_rw("np", &BremPrecomputed::np)
        .def_rw("n_mean", &BremPrecomputed::n_mean)
        .def_rw("T", &BremPrecomputed::T)
        .def_rw("medium_corr", &BremPrecomputed::medium_corr);

    // --- SpectralIntegrandParams: has_* bools def_rw, g_nu_fixed[4] def,
    // grey_pars def_rw ---
    nb::class_<SpectralIntegrandParams>(m, "SpectralIntegrandParams")
        .def(nb::init<>())
        .def_rw("grey_pars", &SpectralIntegrandParams::grey_pars)
        .def_rw("pair_pre", &SpectralIntegrandParams::pair_pre)
        .def_rw("pair_fdi_omega", &SpectralIntegrandParams::pair_fdi_omega)
        .def_rw("neps_pre", &SpectralIntegrandParams::neps_pre)
        .def_rw("brem_pre", &SpectralIntegrandParams::brem_pre)
        .def("g_nu_fixed",
             [](SpectralIntegrandParams& obj) -> Arr4
             {
                 auto* ref = new shared_ref<SpectralIntegrandParams>{
                     std::make_shared<SpectralIntegrandParams>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r =
                             static_cast<shared_ref<SpectralIntegrandParams>*>(
                                 p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->g_nu_fixed, {total_num_species},
                             std::move(owner));
             })
        .def_rw("has_pair_pre", &SpectralIntegrandParams::has_pair_pre)
        .def_rw("has_neps_pre", &SpectralIntegrandParams::has_neps_pre)
        .def_rw("has_brem_pre", &SpectralIntegrandParams::has_brem_pre)
        .def_rw("has_g_nu_fixed", &SpectralIntegrandParams::has_g_nu_fixed);

    // --- M1OpacitiesNonThermalSeparated: 9 arrays ---
    nb::class_<M1OpacitiesNonThermalSeparated>(m,
                                               "M1OpacitiesNonThermalSeparated")
        .def(nb::init<>())
        .def("eta_0",
             [](M1OpacitiesNonThermalSeparated& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1OpacitiesNonThermalSeparated>{
                     std::make_shared<M1OpacitiesNonThermalSeparated>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<
                             shared_ref<M1OpacitiesNonThermalSeparated>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->eta_0, {total_num_species},
                             std::move(owner));
             })
        .def("kappa_0_a",
             [](M1OpacitiesNonThermalSeparated& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1OpacitiesNonThermalSeparated>{
                     std::make_shared<M1OpacitiesNonThermalSeparated>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<
                             shared_ref<M1OpacitiesNonThermalSeparated>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->kappa_0_a, {total_num_species},
                             std::move(owner));
             })
        .def("eta",
             [](M1OpacitiesNonThermalSeparated& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1OpacitiesNonThermalSeparated>{
                     std::make_shared<M1OpacitiesNonThermalSeparated>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<
                             shared_ref<M1OpacitiesNonThermalSeparated>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->eta, {total_num_species},
                             std::move(owner));
             })
        .def("kappa_a",
             [](M1OpacitiesNonThermalSeparated& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1OpacitiesNonThermalSeparated>{
                     std::make_shared<M1OpacitiesNonThermalSeparated>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<
                             shared_ref<M1OpacitiesNonThermalSeparated>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->kappa_a, {total_num_species},
                             std::move(owner));
             })
        .def("kappa_s",
             [](M1OpacitiesNonThermalSeparated& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1OpacitiesNonThermalSeparated>{
                     std::make_shared<M1OpacitiesNonThermalSeparated>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<
                             shared_ref<M1OpacitiesNonThermalSeparated>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->kappa_s, {total_num_species},
                             std::move(owner));
             })
        .def("eta_th",
             [](M1OpacitiesNonThermalSeparated& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1OpacitiesNonThermalSeparated>{
                     std::make_shared<M1OpacitiesNonThermalSeparated>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<
                             shared_ref<M1OpacitiesNonThermalSeparated>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->eta_th, {total_num_species},
                             std::move(owner));
             })
        .def("eta_non_th",
             [](M1OpacitiesNonThermalSeparated& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1OpacitiesNonThermalSeparated>{
                     std::make_shared<M1OpacitiesNonThermalSeparated>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<
                             shared_ref<M1OpacitiesNonThermalSeparated>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->eta_non_th, {total_num_species},
                             std::move(owner));
             })
        .def("kappa_a_th",
             [](M1OpacitiesNonThermalSeparated& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1OpacitiesNonThermalSeparated>{
                     std::make_shared<M1OpacitiesNonThermalSeparated>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<
                             shared_ref<M1OpacitiesNonThermalSeparated>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->kappa_a_th, {total_num_species},
                             std::move(owner));
             })
        .def("kappa_a_non_th",
             [](M1OpacitiesNonThermalSeparated& obj) -> Arr4
             {
                 auto* ref = new shared_ref<M1OpacitiesNonThermalSeparated>{
                     std::make_shared<M1OpacitiesNonThermalSeparated>(obj), 1};
                 nb::capsule owner(
                     ref,
                     [](void* p) noexcept
                     {
                         auto* r = static_cast<
                             shared_ref<M1OpacitiesNonThermalSeparated>*>(p);
                         if (--r->refs == 0)
                             delete r;
                     });
                 return Arr4(ref->data->kappa_a_non_th, {total_num_species},
                             std::move(owner));
             });

    // ==============================================================
    // FUNCTIONS
    // ==============================================================

    m.def(
        "GaussLegendre", [](MyQuadrature& quad) { GaussLegendre(&quad); },
        "quad"_a);

    m.def(
        "NuEquilibriumParams",
        [](const MyEOSParams& eos) -> NuDistributionParams
        { return NuEquilibriumParams(&eos); }, "eos_pars"_a);

    m.def(
        "CalculateDistrParamsFromM1",
        [](const M1Quantities& m1,
           const MyEOSParams& eos) -> NuDistributionParams
        { return CalculateDistrParamsFromM1(&m1, &eos); },
        "m1_pars"_a, "eos_pars"_a);

    m.def(
        "ComputeM1DensitiesEq",
        [](const MyEOSParams& eos, const NuDistributionParams& distr,
           M1Quantities& m1) { ComputeM1DensitiesEq(&eos, &distr, &m1); },
        "eos_pars"_a, "distr_pars"_a, "m1_pars"_a);

    m.def(
        "TotalNuF",
        [](BS_REAL omega, const NuDistributionParams& distr,
           int nuid) -> BS_REAL { return TotalNuF(omega, &distr, nuid); },
        "omega"_a, "distr_pars"_a, "nuid"_a);

    m.def(
        "ComputeSpectralOpacitiesNotStimulatedAbs",
        [](BS_REAL nu, MyQuadrature& quad_1d,
           GreyOpacityParams& params) -> SpectralOpacities
        {
            return ComputeSpectralOpacitiesNotStimulatedAbs(nu, &quad_1d,
                                                            &params);
        },
        "nu"_a, "quad_1d"_a, "grey_opacity_params"_a);

    m.def(
        "ComputeSpectralOpacitiesStimulatedAbs",
        [](BS_REAL nu, MyQuadrature& quad_1d,
           GreyOpacityParams& params) -> SpectralOpacities
        {
            return ComputeSpectralOpacitiesStimulatedAbs(nu, &quad_1d, &params);
        },
        "nu"_a, "quad_1d"_a, "grey_opacity_params"_a);

    m.def(
        "ComputeM1Opacities",
        [](const MyQuadrature& quad_1d, const MyQuadrature& quad_2d,
           GreyOpacityParams& params) -> M1Opacities
        { return ComputeM1Opacities(&quad_1d, &quad_2d, &params); },
        "quad_1d"_a, "quad_2d"_a, "grey_opacity_params"_a);

    m.def(
        "ComputeM1OpacitiesNotStimulated",
        [](MyQuadrature& quad_1d, MyQuadrature& quad_2d,
           GreyOpacityParams& params) -> M1Opacities
        {
            return ComputeM1OpacitiesNotStimulated(&quad_1d, &quad_2d, &params);
        },
        "quad_1d"_a, "quad_2d"_a, "grey_opacity_params"_a);

    m.def(
        "ComputeM1OpacitiesNonThermalSeparated",
        [](MyQuadrature& quad_1d, MyQuadrature& quad_2d,
           GreyOpacityParams& params) -> M1OpacitiesNonThermalSeparated
        {
            return ComputeM1OpacitiesNonThermalSeparated(&quad_1d, &quad_2d,
                                                         &params);
        },
        "quad_1d"_a, "quad_2d"_a, "grey_opacity_params"_a);

    // ==============================================================
    // DEFAULT GLOBALS
    // ==============================================================
    m.attr("quadrature_default")         = nb::cast(quadrature_default);
    m.attr("opacity_params_default_all") = nb::cast(opacity_params_default_all);
    m.attr("opacity_params_default_none") =
        nb::cast(opacity_params_default_none);
    m.attr("opacity_flags_default_all")  = nb::cast(opacity_flags_default_all);
    m.attr("opacity_flags_default_none") = nb::cast(opacity_flags_default_none);
}
