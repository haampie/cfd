function ex3a()
  close all;

  U = 50;
  L = 1;
  k = 0.1;
  Pe = U * L / k;

  phi_0 = 0;
  phi_L = 20;

  N_values = floor(10 .^ (1 : 0.1 : 5));
  h_values = 1 ./ N_values;
  error_values = zeros(2, length(N_values));

  for idx = 1 : length(N_values)
    [x, y_central] = convection_diffusion_dirichlet(U, L, k, N_values(idx), phi_0, phi_L, 'central');
    [~, y_upwind] = convection_diffusion_dirichlet(U, L, k, N_values(idx), phi_0, phi_L, 'upwind');
    y_exact = phi_0 + (phi_L - phi_0) * (1 - exp(Pe * x' ./ L)) ./ (1 - exp(Pe));
    
    error_values(1, idx) = norm(y_central - y_exact) / sqrt(N_values(idx));
    error_values(2, idx) = norm(y_upwind - y_exact) / sqrt(N_values(idx));
  end

  figure;
  loglog(h_values, error_values(1, :), 'b-*');
  hold on;
  loglog(h_values, error_values(2, :), 'r-+');
end