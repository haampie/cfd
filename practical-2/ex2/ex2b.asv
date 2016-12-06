function ex2b()
  close all;

  % Fix U, L and k
  U = 50;
  L = 1;
  k = 1;
  Pe = U * L / k;

  % Dirichlet boundary conditions.
  phi_0 = 0;
  phi_L = 20;

  % But vary h
  N_values = floor(10 .^ (2 : -0.05 : 1));
  alpha = 0.7;
  
  % Holds the error for different values of N / h
  error_values = zeros(1, length(N_values));

  for idx = 1 : length(N_values)
    grid = (-alpha.^(2 : N_values(idx) + 1) + alpha) / (1 - alpha);
    grid = grid / grid(end);
    grid = [0 grid];
        
    [x, y_central] = convection_diffusion_dirichlet(U, L, k, N_values(idx), phi_0, phi_L, 'central', grid);
    y_exact = phi_0 + (phi_L - phi_0) * (1 - exp(Pe * x' ./ L)) ./ (1 - exp(Pe));
    
    plot(x, y_central); hold on
    
    h_values = grid(2 : end) - grid(1 : end - 1);
    error_values(idx) = sqrt(sum(h_values .* (y_central(1 : end - 1)' - y_exact(1 : end - 1)').^2));
  end
  
  figure;
  loglog(N_values, error_values(1, :), 'b-*');
  legend('central')
  
  diff_e = log(error_values(2 : end)) - log(error_values(1 : end - 1));
  diff_N = log(N_values(2 : end)) - log(N_values(1 : end - 1));
  diff_e ./ diff_N
end