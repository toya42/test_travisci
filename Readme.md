## Test_travisci [![Build Status](https://travis-ci.org/toya42/test_travisci.svg?branch=master)](https://travis-ci.org/toya42/test_travisci)

Governing equation 
   Two-dimensional vorticity equation
Boundary condition
   doubly-eriodic

Spacial discretization 
   Fourier-Galerkin
Time discretization

​    nonlinear term: Traditional 4 stage 4th order Runge-Kutta method

​    linear term: Integrating factor technique

## Dependeicy

   Fortran2008
   Intel MKL
   cmake(3.14.5 or later)

## License

This software is released under the MIT License, see LICENSE.

## Authors

Twitter:

@toya42_fortran

## References

C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Zand, *Specral Methods Fundamentals in Single Domain*, Springer-Verlag (2006).

C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Zand, *Specral Methods Evolution to Complex Geometries and Applications to Fluid Dynamics*, Springer-Verlag (2007).