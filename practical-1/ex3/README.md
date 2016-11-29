# Exercise 3

(a) run `ex3a`
(b) Let $e(h) := x(h) - x$ be the discrete error, where $x(h)$ is the numerical solution and $x$ the true solution. Suppose $h_2 < h_1$, then

$$|x(h1) - x(h2)| = |e(h1) - e(h2)| = O(e(h1))$$

Assume $e(h) = O(h^p)$, then $$|x(h1) - x(h2)| / |x(h2) - x(h3)| ~= (h1 / h2)^p$$

So the rate of convergence $p$ is asymptotically

$$p = log(|x(h1) - x(h2)| / |x(h2) - x(h3)|) / log(h1 / h2)$$

 