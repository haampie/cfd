function ex1a()

  % [x, y] = convection_diffusion_dirichlet(10, 1, 1, 1000, 0, 1, 'b3');

  my_grid = 1 - (2 .^ (5 : -0.1 : 0) - 1) / 31;

  [xx, yy] = convection_diffusion_dirichlet(100, 1, 1, 50, 0, 1, 'central', my_grid);
  plot(xx, yy, '*');
end