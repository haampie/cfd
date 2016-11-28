function ex2b()

  [x, y] = convection_diffusion_dirichlet(50, 10, 1, 250, 0, 20, 'upwind');
  plot(x, y);

end