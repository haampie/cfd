# Exercise 4

(a) When k goes to 0 the numerical solution does not seem to converge to the analytical solution. The problem is in fact ill-posed, since for k = 0, the ODE reduces to first-order and cannot meet both boundary conditions. However, if k goes to 0, one would expect the numerical solution to converge to solution of the first-order ODE, with a jump/discontinuity in 1 to meet the r.h.s. boundary condition.

This is not the case however, due to artificial diffusion.

(b) If upwind discretization is used on the operator Ly = y' - ky'', then it is equivalent to applying central discretization to the operator My = y' - (k + h / 2)y''. So if k is small while h is fixed, the numerical solution of the upwind scheme coincides with the analytical solution to Ly = S for k = h / 2 = 1 / (2 * N) = 1 / 40 = 0.025.

(c) Now k = 1 / (2 * N) = 1 / 80 = 1 / 80 = 0.0125.

