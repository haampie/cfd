function ex2a()
  close all;

  % Fix U, L and k
  U = 50;
  L = 1;
  k = 0.1;
  Pe = U * L / k;

  % Dirichlet boundary conditions.
  phi_0 = 0;
  phi_L = 20;

  % But vary h
  N_values = floor(10 .^ (4 : -0.1 : 1));
  h_values = L ./ N_values;
  
  % Holds the error for different values of N / h
  error_values = zeros(1, length(N_values));

  for idx = 1 : length(N_values)
    [x, y_central] = convection_diffusion_dirichlet(U, L, k, ...
        N_values(idx), phi_0, phi_L, 'b3');
    y_exact = phi_0 + (phi_L - phi_0) * (1 - exp(Pe * x' ./ L)) ./ (1 - exp(Pe));
    error_values(idx) = norm(y_central - y_exact) * sqrt(h_values(idx));
  end

  figure;
  loglog(h_values, error_values(1, :), 'b-*');
  title(sprintf('U = %d, k = %2.2f, L = %d, phi_0 = %d, phi_L = %d', ...
      U, k, L, phi_0, phi_L))
  xlabel('h')
  ylabel('error')
  legend('b3')
  
  % Dump the slope (~= 2)
  diff_e = log(error_values(2 : end)) - log(error_values(1 : end - 1));
  diff_h = log(h_values(2 : end)) - log(h_values(1 : end - 1));
  diff_e ./ diff_h
end