## Exercise 2

(a) We must have $P \equiv \frac{uh}{k} = \frac{u L}{k N} \le 2$
So $N \ge \frac{U L}{2k} = 250 000$
(b) Yes, since the off-diagonal elements are strictly < 0. If U >= 0, then the lower diagonal element is -P - 1 < 0, the upper diagonal element is -1 < 0, and the diagonal element is P + 2 > 0. If U < 0, then P < 0, so it is completely similar.
With regard to the boundary layer thickness, we see in `exc2c.m` that the boundary layer thickness grows linearly with `h`.

Via lemma 1.2.2 we see that for large (numerical) P we have phi_i = P^(i - N) if the boundary conditions are phi(0) = 0 and phi(1) = 1. So phi_(N-i) = 1 / P^i approximately.

(c) No difference in wiggles, since only the last coefficient in the diagonal has changed from 2 to 1 + P / 2 (still positive).