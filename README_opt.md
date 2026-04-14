# opt
Attempt at some minor optimizations.

```bash
# compile n' run cut'n'pasta

# env
activate_spack
spack env activate num_openmpi_gcc

# clean compile
rm -rf build
mkdir -p build
cd build

cmake -DENABLE_KOKKOS=OFF ..
cmake --build . --target MWE_profile

# run
./MWE_profile 100

cd ..
```

## Pre:
```bash
Input parameters
----------------
Neutrino energy (for spectral rates)               :  1.000000e+01 (MeV)
Baryon number density                              :  4.208367e+38 (cm^-3)
Temperature                                        :  1.240640e+01 (MeV)
Electron fraction                                  :         0.072
Relativistic electron chemical potential           :  1.871814e+02 (MeV)
Relativistic proton chemical potential             :  1.011018e+03 (MeV)
Relativistic neutron chemical potential            :  1.221590e+03 (MeV)
Effective nucleon mass difference (n - p)          :  1.293333e+00 (MeV)
Effective nucleon potential difference (n - p)     :  1.892715e+01 (MeV)
Electron neutrinos number density 'n'              :  3.739749e+33 (cm^-3)
Electron antineutrinos number density 'n'          :  1.217496e+35 (cm^-3)
Heavy-type neutrinos number density 'n'            :  2.243850e+34 (cm^-3)
Heavy-type antineutrinos number density 'n'        :  2.243850e+34 (cm^-3)
Electron neutrinos energy density 'J'              :  1.246583e+35 (MeV cm^-3)
Electron antineutrinos energy density 'J'          :  5.360307e+36 (MeV cm^-3)
Heavy-type neutrinos energy density 'J'            :  8.726082e+35 (MeV cm^-3)
Heavy-type antineutrinos energy density 'J'        :  8.726082e+35 (MeV cm^-3)
Electron neutrinos Eddington parameter 'chi'       : 0.33333333333
Electron antineutrinos Eddington parameter 'chi'   : 0.33333333333
Heavy-type neutrinos Eddington parameter 'chi'     : 0.33333333333
Heavy-type antineutrinos Eddington parameter 'chi' : 0.33333333333

Reconstructed neutrino densities assuming equilibrium
-----------------------------------------------------
    nue            anue            nux            anux
n   3.752388e+33   1.099744e+35    2.270204e+34   2.270204e+34     (cm^-3)
J   1.409347e+35   4.867722e+36    8.875866e+35   8.875866e+35 (MeV cm^-3)
chi 0.3333333333   0.3333333333    0.3333333333   0.3333333333

Spectral rates assuming equilibrium
------------------------------
     j             j_s           kappa         kappa_s
 nue 2.281657e+05  3.761940e+05  1.122778e-04  1.851209e-04 
anue 1.681473e+06  4.422965e+06  1.905980e-05  5.013511e-05 
 nux 2.912999e+05  1.829562e+06  2.175586e-05  1.366417e-04 
anux 2.882248e+05  1.829562e+06  2.152620e-05  1.366417e-04 

Gray rates assuming equilibrium
------------------------------
     eta0          eta1          kappa0        kappa1        scat1
 nue 9.274316e+40  1.002239e+43  8.244294e-04  2.372098e-03  6.313839e-03 
anue 8.238882e+41  4.928994e+43  2.498941e-04  3.377628e-04  7.355463e-03 
 nux 4.864064e+40  2.659807e+42  7.146835e-05  9.995824e-05  6.533912e-03 
anux 4.633231e+40  2.498959e+42  6.807669e-05  9.391341e-05  6.533912e-03 

Spectral rates reconstructing distribution function
------------------------------
     j             j_s           kappa         kappa_s
 nue 2.277815e+05  4.826419e+05  1.200658e-04  1.815701e-04 
anue 1.825383e+06  4.669667e+06  1.777519e-05  4.190601e-05 
 nux 2.907223e+05  1.830833e+06  2.158623e-05  1.365993e-04 
anux 2.876876e+05  1.830833e+06  2.135649e-05  1.365993e-04 

Gray rates reconstructing distribution function
----------------------------------------------
     eta0          eta1          kappa0        kappa1        scat1
 nue 8.591723e+40  9.484546e+42  4.602474e-04  1.058498e-03  4.980347e-03 
anue 8.981812e+41  5.310832e+43  2.464567e-04  3.276931e-04  7.125131e-03 
 nux 4.783883e+40  2.602177e+42  7.081153e-05  9.892904e-05  6.464914e-03 
anux 4.558913e+40  2.445997e+42  6.746718e-05  9.297810e-05  6.464914e-03 


Units
-----
Spectral emissivity 'j'/'j_s'   :           s^-1
Spectral imfp 'kappa'/'kappa_s' :          cm^-1
Gray number emissivity 'eta0'   :     cm^-3 s^-1
Gray energy emissivity 'eta1'   : MeV cm^-3 s^-1
Gray number opacity 'kappa0'    :          cm^-1
Gray energy opacity 'kappa1'    :          cm^-1
Gray scattering opacity 'scat1' :          cm^-1

=============================================================================
  PROFILING SUMMARY  (N = 100 iterations)
=============================================================================
Function                                     Total (ms)     Avg (us)  % total
-----------------------------------------------------------------------------
GaussLegendre                                     0.275        2.752    0.12%
NuEquilibriumParams                               0.153        1.530    0.06%
ComputeM1DensitiesEq                              0.172        1.717    0.07%
ComputeSpectralOpacities [eq]                    57.969      579.688   24.47%
    Pair+Brem integration                        37.575      375.751   15.86%
    NEPS integration                             19.181      191.814    8.10%
    Beta processes                                0.194        1.938    0.08%
    Iso scattering                                0.179        1.785    0.08%
ComputeM1Opacities [eq]                          59.561      595.607   25.14%
    Iso scattering (1D)                           2.197       21.973    0.93%
    Beta processes (1D)                           1.260       12.595    0.53%
    Pair+Brem (2D)                               31.996      319.962   13.51%
    NEPS / Inelastic scatt (2D)                  22.988      229.880    9.70%
    Assembly                                      0.159        1.588    0.07%
CalculateDistrParamsFromM1                        0.227        2.267    0.10%
ComputeSpectralOpacities [m1]                    58.536      585.364   24.71%
    Pair+Brem integration                        37.904      379.039   16.00%
    NEPS integration                             19.404      194.039    8.19%
    Beta processes                                0.192        1.916    0.08%
    Iso scattering                                0.183        1.827    0.08%
ComputeM1Opacities [m1]                          59.977      599.773   25.32%
    Iso scattering (1D)                           2.226       22.257    0.94%
    Beta processes (1D)                           1.277       12.774    0.54%
    Pair+Brem (2D)                               32.125      321.248   13.56%
    NEPS / Inelastic scatt (2D)                  23.234      232.338    9.81%
    Assembly                                      0.164        1.638    0.07%
-----------------------------------------------------------------------------
TOTAL                                           236.870     2368.698  100.00%
=============================================================================
```

## Pos:
```bash
Reconstructed neutrino densities assuming equilibrium
-----------------------------------------------------
    nue            anue            nux            anux
n   3.752388e+33   1.099744e+35    2.270204e+34   2.270204e+34     (cm^-3)
J   1.409347e+35   4.867722e+36    8.875866e+35   8.875866e+35 (MeV cm^-3)
chi 0.3333333333   0.3333333333    0.3333333333   0.3333333333

Spectral rates assuming equilibrium
------------------------------
     j             j_s           kappa         kappa_s
 nue 2.281657e+05  3.761940e+05  1.122778e-04  1.851209e-04 
anue 1.681473e+06  4.422965e+06  1.905980e-05  5.013511e-05 
 nux 2.912999e+05  1.829562e+06  2.175586e-05  1.366417e-04 
anux 2.882248e+05  1.829562e+06  2.152620e-05  1.366417e-04 

Gray rates assuming equilibrium
------------------------------
     eta0          eta1          kappa0        kappa1        scat1
 nue 9.274316e+40  1.002239e+43  8.244294e-04  2.372098e-03  6.313839e-03 
anue 8.238882e+41  4.928994e+43  2.498941e-04  3.377628e-04  7.355463e-03 
 nux 4.864064e+40  2.659807e+42  7.146835e-05  9.995824e-05  6.533912e-03 
anux 4.633231e+40  2.498959e+42  6.807669e-05  9.391341e-05  6.533912e-03 

Spectral rates reconstructing distribution function
------------------------------
     j             j_s           kappa         kappa_s
 nue 2.277815e+05  4.826419e+05  1.200658e-04  1.815701e-04 
anue 1.825383e+06  4.669667e+06  1.777519e-05  4.190601e-05 
 nux 2.907223e+05  1.830833e+06  2.158623e-05  1.365993e-04 
anux 2.876876e+05  1.830833e+06  2.135649e-05  1.365993e-04 

Gray rates reconstructing distribution function
----------------------------------------------
     eta0          eta1          kappa0        kappa1        scat1
 nue 8.591723e+40  9.484546e+42  4.602474e-04  1.058498e-03  4.980347e-03 
anue 8.981812e+41  5.310832e+43  2.464567e-04  3.276931e-04  7.125131e-03 
 nux 4.783883e+40  2.602177e+42  7.081153e-05  9.892904e-05  6.464914e-03 
anux 4.558913e+40  2.445997e+42  6.746718e-05  9.297810e-05  6.464914e-03 


Units
-----
Spectral emissivity 'j'/'j_s'   :           s^-1
Spectral imfp 'kappa'/'kappa_s' :          cm^-1
Gray number emissivity 'eta0'   :     cm^-3 s^-1
Gray energy emissivity 'eta1'   : MeV cm^-3 s^-1
Gray number opacity 'kappa0'    :          cm^-1
Gray energy opacity 'kappa1'    :          cm^-1
Gray scattering opacity 'scat1' :          cm^-1

=============================================================================
  PROFILING SUMMARY  (N = 100 iterations)
=============================================================================
Function                                     Total (ms)     Avg (us)  % total
-----------------------------------------------------------------------------
GaussLegendre                                     0.320        3.202    0.42%
NuEquilibriumParams                               0.165        1.648    0.22%
ComputeM1DensitiesEq                              0.186        1.863    0.25%
ComputeSpectralOpacities [eq]                     6.047       60.471    7.98%
    Pair+Brem integration                         2.620       26.199    3.46%
    NEPS integration                              1.931       19.313    2.55%
    Beta processes                                0.207        2.067    0.27%
    Iso scattering                                0.193        1.930    0.25%
ComputeM1Opacities [eq]                          31.443      314.428   41.51%
    Iso scattering (1D)                           2.146       21.459    2.83%
    Beta processes (1D)                           1.254       12.544    1.66%
    Pair+Brem (2D)                               13.843      138.432   18.28%
    NEPS / Inelastic scatt (2D)                  12.968      129.683   17.12%
    Assembly                                      0.180        1.801    0.24%
CalculateDistrParamsFromM1                        0.243        2.430    0.32%
ComputeSpectralOpacities [m1]                     6.107       61.066    8.06%
    Pair+Brem integration                         2.666       26.663    3.52%
    NEPS integration                              1.958       19.583    2.59%
    Beta processes                                0.204        2.039    0.27%
    Iso scattering                                0.194        1.939    0.26%
ComputeM1Opacities [m1]                          31.231      312.312   41.23%
    Iso scattering (1D)                           2.112       21.118    2.79%
    Beta processes (1D)                           1.268       12.676    1.67%
    Pair+Brem (2D)                               13.765      137.651   18.17%
    NEPS / Inelastic scatt (2D)                  12.885      128.852   17.01%
    Assembly                                      0.168        1.678    0.22%
-----------------------------------------------------------------------------
TOTAL                                            75.742      757.421  100.00%
=============================================================================
```