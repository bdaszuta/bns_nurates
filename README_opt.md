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
GaussLegendre                                     0.288        2.882    0.13%
NuEquilibriumParams                               0.155        1.548    0.07%
ComputeM1DensitiesEq                              0.172        1.717    0.08%
ComputeSpectralOpacities [eq]                    53.368      533.676   24.46%
    Pair+Brem integration                         0.000        0.000    0.00%
    NEPS integration                              0.000        0.000    0.00%
    Beta processes                                0.000        0.000    0.00%
    Iso scattering                                0.000        0.000    0.00%
ComputeM1Opacities [eq]                          55.125      551.248   25.26%
    Iso scattering (1D)                           0.000        0.000    0.00%
    Beta processes (1D)                           0.000        0.000    0.00%
    Pair+Brem (2D)                                0.000        0.000    0.00%
    NEPS / Inelastic scatt (2D)                   0.000        0.000    0.00%
    Assembly                                      0.000        0.000    0.00%
CalculateDistrParamsFromM1                        0.240        2.402    0.11%
ComputeSpectralOpacities [m1]                    53.876      538.759   24.69%
    Pair+Brem integration                         0.000        0.000    0.00%
    NEPS integration                              0.000        0.000    0.00%
    Beta processes                                0.000        0.000    0.00%
    Iso scattering                                0.000        0.000    0.00%
ComputeM1Opacities [m1]                          54.981      549.814   25.20%
    Iso scattering (1D)                           0.000        0.000    0.00%
    Beta processes (1D)                           0.000        0.000    0.00%
    Pair+Brem (2D)                                0.000        0.000    0.00%
    NEPS / Inelastic scatt (2D)                   0.000        0.000    0.00%
    Assembly                                      0.000        0.000    0.00%
-----------------------------------------------------------------------------
TOTAL                                           218.205     2182.045  100.00%
=============================================================================
```

## Pos:
```bash
# ...
```