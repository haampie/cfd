function [grid, solution] = convection_diffusion_mixed(U, L, k, N, phi_left, neumann_right, method)

  h = L / N;
  P = U * h / k;

  assert(U > 0)

  switch method
    case 'central'
      low = -P / 2 - 1;
      mid = 2;
      up = P / 2 - 1;
    case 'upwind'
      low = -P - 1;
      mid = P + 2;
      up = -1;
    otherwise
      error('central and upwind are only supported');
  end

  rhs = zeros(N, 1);
  rhs(1) = -low * phi_left;
  rhs(N) = h * neumann_right;

  sub_diagonal = low * ones(N, 1);
  diagonal = mid * ones(N, 1);
  super_diagonal = up * ones(N, 1);
  
  sub_diagonal(N) = -1;
  diagonal(N) = 1;

  solution = [
    phi_left; 
    tdma(sub_diagonal, diagonal, super_diagonal, rhs, N); 
  ];

  grid = linspace(0, L, N + 1);

end