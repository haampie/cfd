function [grid, solution] = parabolic_dirichlet(U, L, k, N, ...
    T, dt, initial_condition, phi_left, phi_right, method, varargin)

  h = L / N;
  P = U * h / k;
  k_t_over_x_x = k * dt / h^2;

  grid = linspace(0, L, N + 1);

  subsub_diagonal = zeros(N - 1, 1);
  sub_diagonal = zeros(N - 1, 1);
  diagonal = zeros(N - 1, 1);
  super_diagonal = zeros(N - 1, 1);
  
  switch method
    case 'central'
      if isempty(varargin)
        
        % Equidistant grid
        low = -P / 2 - 1;
        mid = 2;
        up = P / 2 - 1;

        sub_diagonal = low * ones(N - 1, 1);
        diagonal = mid * ones(N - 1, 1);
        super_diagonal = up * ones(N - 1, 1);
      else
        
        % Non-equidistant grid
        grid = varargin{1};
        N = length(grid) - 1;
        h = grid(2 : N + 1) - grid(1 : N);
        for i = 1 : N - 1
          sub_diagonal(i) = -h(i) * h(i + 1) * U - 2 * k * h(i + 1);
          diagonal(i) = 2 * k * (h(i) + h(i + 1));
          super_diagonal(i) = h(i) * h(i + 1) * U - 2 * k * h(i);
        end
      end

    case 'upwind'
      low = -P - 1;
      mid = P + 2;
      up = -1;

      sub_diagonal = low * ones(N - 1, 1);
      diagonal = mid * ones(N - 1, 1);
      super_diagonal = up * ones(N - 1, 1);

    case 'b3'
      lowlow = P / 2;
      low = -1 - 2 * P;
      mid = 3 * P / 2 + 2;
      up = -1;

      subsub_diagonal = lowlow * ones(N - 1, 1);
      sub_diagonal = low * ones(N - 1, 1);
      diagonal = mid * ones(N - 1, 1);
      super_diagonal = up * ones(N - 1, 1);

      % Central discretization at the left boundary 
      diagonal(1) = 2;
      super_diagonal(1) = -1 + P / 2;
    otherwise
      error('central/upwind/b3 are only supported');
  end
  
  solution = zeros(floor(T / dt), N + 1);
  solution(1, :) = initial_condition;
  
  for time_idx = 2 : size(solution, 1)
      solution(time_idx, 1) = phi_left;
      solution(time_idx, 2) = k_t_over_x_x * sub_diagonal(1) * solution(time_idx - 1, 1) + ...
          (1 + k_t_over_x_x * diagonal(2)) * solution(time_idx - 1, 2) + ...
          k_t_over_x_x * super_diagonal(3) * solution(time_idx - 1, 3);

      solution(time_idx, N + 1) = phi_right;
      
      % Compute the rest of the solution.
      for sol_idx = 3 : N
          solution(time_idx, sol_idx) = ...
              k_t_over_x_x * subsub_diagonal(sol_idx - 1) * solution(time_idx - 1, sol_idx - 2) + ...
              k_t_over_x_x * sub_diagonal(sol_idx - 1) * solution(time_idx - 1, sol_idx - 1) + ...
              (1 + k_t_over_x_x * diagonal(sol_idx - 1)) * solution(time_idx - 1, sol_idx) + ...
              k_t_over_x_x * super_diagonal(sol_idx - 1) * solution(time_idx - 1, sol_idx + 1);
      end
  end
end