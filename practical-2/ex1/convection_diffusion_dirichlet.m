function [grid, solution] = convection_diffision_dirichlet(U, L, k, N, phi_left, phi_right, method, varargin)

  h = L / N;
  P = U * h / k;
  % physical_peclet = U * L / k;

  rhs = zeros(N - 1, 1);
  grid = linspace(0, L, N + 1);

  subsub_diagonal = zeros(N - 1, 1);
  sub_diagonal = zeros(N - 1, 1);
  diagonal = zeros(N - 1, 1);
  super_diagonal = zeros(N - 1, 1);
  supersuper_diagonal = zeros(N - 1, 1);

  switch method
    case 'central'
      if isempty(varargin)
        
        % Equidistant grid
        low = -P / 2 - 1;
        mid = 2;
        up = P / 2 - 1;

        % Same discretization at the boundary
        rhs(1) = -low * phi_left;
        rhs(N - 1) = -up * phi_right;

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

        rhs(1) = -sub_diagonal(1) * phi_left;
        rhs(N - 1) = -super_diagonal(N - 1) * phi_right;
      end

    case 'upwind'
      low = -P - 1;
      mid = P + 2;
      up = -1;

      % Same discretization at the boundary
      rhs(1) = -low * phi_left;
      rhs(N - 1) = -up * phi_right;

      sub_diagonal = low * ones(N - 1, 1);
      diagonal = mid * ones(N - 1, 1);
      super_diagonal = up * ones(N - 1, 1);

    case 'b3'
      lowlow = P / 2;
      low = -1 - 2 * P;
      mid = 3 * P / 2 + 2;
      up = -1;

      % Central discretization at the left boundary
      rhs(1) = (1 + P / 2) * phi_left;
      rhs(2) = -lowlow * phi_left;
      rhs(N - 1) = -up * phi_right;

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

  solution = [
    phi_left; 
    pdma(subsub_diagonal, sub_diagonal, diagonal, super_diagonal, supersuper_diagonal, rhs, N - 1); 
    phi_right
  ];

  
end