function ex2c()

  N = 2500;
  [x, y_upwind] = convection_diffusion_mixed(50, 1, 1, N, 0, 10, 'upwind');
  [~, y_central] = convection_diffusion_mixed(50, 1, 1, N, 0, 10, 'central');
  
  plot(x, y_upwind);
  hold on;
  plot(x, y_central);
  hold off;
end