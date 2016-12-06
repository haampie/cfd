function ex1a()
  close all
  
  U = 10;
  L = 1;
  k = 1;
  Pe = U * L / k;

  [x, y] = convection_diffusion_dirichlet(U, L, k, 500, 0, 1, 'b3');
  my_grid = 1 - (2 .^ (5 : -0.1 : 0) - 1) / 31;
  
  [xx, yy] = convection_diffusion_dirichlet(U, L, k, 50, 0, 1, 'central', my_grid);
  y_exact = (1 - exp(Pe * xx' ./ L)) ./ (1 - exp(Pe));
  
  figure;
  plot(xx, y_exact, '-'); hold on 
  plot(x, y, '-+'); hold on;
  plot(xx, yy, '-*'); hold off;
  legend('exact', 'central', 'b3');
end