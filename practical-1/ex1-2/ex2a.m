function ex2a()

  [x, y] = convection_diffusion_dirichlet(50, 10, 1e-3, 250000, 0, 20, 'central');
  plot(x, y);

end