function ex2b()
  close all;

  % Fix U, L and k
  U = 50;
  L = 1;
  k = 0.1;
  Pe = U * L / k;

  % Dirichlet boundary conditions.
  phi_0 = 0;
  phi_L = 1;

  % But vary h
  N_values = unique(floor(10 .^ (1 : 0.03 : 2)), 'sorted');
  alpha_values = [0.5 0.7 0.95];
  
  % Holds the error for different values of N / h
  error_values = zeros(length(alpha_values), length(N_values));

  % Calculate the error for multiple values of alpha
  for alpha_idx = 1 : length(alpha_values)
      
      % And multiple values of N
      for N_idx = 1 : length(N_values)
          
        % Construct an exponentially refined mesh
        the_grid = (-alpha_values(alpha_idx).^(2 : N_values(N_idx) + 1) + alpha_values(alpha_idx)) / (1 - alpha_values(alpha_idx));
        the_grid = the_grid / the_grid(end);
        the_grid = [0 the_grid];

        % Compute the numerical solution
        [x, y_central] = convection_diffusion_dirichlet(U, L, k, N_values(N_idx), phi_0, phi_L, 'central', the_grid);
        
        % Compute the exact solution
        y_exact = phi_0 + (phi_L - phi_0) * (1 - exp(Pe * x' ./ L)) ./ (1 - exp(Pe));

        % Compute the error (vectorized for speed)
        h_values = the_grid(2 : end) - the_grid(1 : end - 1);
        error_values(alpha_idx, N_idx) = sqrt(sum(h_values .* (y_central(1 : end - 1)' - y_exact(1 : end - 1)').^2));
      end
  end
  
  % Plot the errors
  for row_idx = 1 : length(alpha_values)
      figure;
      loglog(N_values, error_values(row_idx, :), '-*'); hold on;
      grid on
      title(sprintf('U = %d, k = %2.2f, L = %d, alpha = %2.2f', U, k, L, alpha_values(row_idx)));
      xlabel('N')
      ylabel('Error approximation')
      
      diff_e = log(error_values(row_idx, 2 : end)) - log(error_values(row_idx, 1 : end - 1));
      diff_N = log(N_values(2 : end)) - log(N_values(1 : end - 1));
      diff_e ./ diff_N
  end
end