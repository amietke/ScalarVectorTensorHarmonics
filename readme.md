## Library and projection tools for scalar, vector and tensor spherical harmonics

<div style="text-align: justify">

Library based on real-valued scalar spherical harmonics Y_lm without the Condon-Shortley phase. 
From these, vector spherical harmonics, and symmetric trace-free tensor 
spherical harmonics as defined in [1] are constructed and can all together be called from `SVTH.m`
on an arbitrary set of points (including the two coordinate poles) on the unit sphere. 

Run `PlotHarmonics.m` to plot all the corresponding fields for some mode (l,m) on a spherical surface.

Run `HarmonicProjectionExamples.m` to see examples of transformations between mode representations
and real space representations for each field. For given harmonic mode coefficients, real-space fields
are simply constructed as linear combinations of the fields coming from `SVTH.m`. Inverse transforms, from an
arbitrary real-space field on the sphere to its harmonic mode coefficients, are implemented as least-square regressions that
generalize the approach from [2] to vector and tensor spherical harmonics.


[1] V. D. Sandberg. Tensor spherical harmonics on S2 and S3 as eigenvalue problems. J. Math. Phys., 19(12): 2441–2446 (1978) [https://doi.org/10.1063/1.523649](https://doi.org/10.1063/1.523649) <br>
[2] N. Sneeuw. Global spherical harmonic analysis by least-squares and numerical quadrature methods in historical perspective. Geophys. J. Int., 118(3): 707–716 (1994) [https://doi.org/10.1111/j.1365-246X.1994.tb03995.x](https://doi.org/10.1111/j.1365-246X.1994.tb03995.x)

</div>

