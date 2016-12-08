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
  N = unique(floor(10 .^ (1 : 0.05 : 2)));
  alpha = [0.5 0.6 0.7 0.8 0.95 1];
  
  % Holds the error for different values of N / h
  errors = zeros(length(alpha), length(N));

  % Calculate the error for multiple values of alpha
  for alpha_idx = 1 : length(alpha)
      
      % And multiple values of N
      for N_idx = 1 : length(N)
          
        % Compute the numerical solution
        if alpha(alpha_idx) == 1
            % Equidistant grid.
            the_grid = linspace(0, 1, 1 + N(N_idx));
            
            % Compute numerical solution
            [x, y_central] = convection_diffusion_dirichlet(U, L, k, ...
                N(N_idx), phi_0, phi_L, 'central');
        else
            % Construct an exponentially refined mesh
            the_grid = (-alpha(alpha_idx).^(2 : N(N_idx) + 1) ...
                + alpha(alpha_idx)) / (1 - alpha(alpha_idx));
            the_grid = the_grid / the_grid(end);
            the_grid = [0 the_grid];
            
            % Compute numerical solution
            [x, y_central] = convection_diffusion_dirichlet(U, L, k, ...
                N(N_idx), phi_0, phi_L, 'central', the_grid);
        end
        
        % Compute the exact solution
        y_exact = phi_0 + (phi_L - phi_0) * (1 - exp(Pe * x' ./ L)) ./ (1 - exp(Pe));

        % Compute the error (vectorized)
        h_values = the_grid(2 : end) - the_grid(1 : end - 1);
        errors(alpha_idx, N_idx) = sqrt(sum(h_values .* ...
            (y_central(1 : end - 1)' - y_exact(1 : end - 1)').^2));
      end
  end
  
  % Plot the errors
  for idx = 1 : length(alpha)
      
      % Plot N vs error.
      figure;
      loglog(N, errors(idx, :), '-*');
      grid on
      title(sprintf('U = %d, k = %2.2f, L = %d, alpha = %2.2f', U, k, ...
          L, alpha(idx)));
      xlabel('N')
      ylabel('Error approximation')
      
      % Show the slope
      fprintf('Alpha = %f', alpha(idx))
      diff_e = log(errors(idx, 2 : end)) - log(errors(idx, 1 : end - 1));
      diff_N = log(N(2 : end)) - log(N(1 : end - 1));
      diff_e ./ diff_N
  end
end