function convection_diffusion_mixed(U, L, k, N, phi_left, neumann_right, method)

  h = L / N;
  P = U * h / k;

  switch method
    case 'central'
      low = -P / 2 - 1;
      mid = 2;
      up = P / 2 - 1;
    case 'upwind'
      if U >= 0
        low = -P - 1;
        mid = P + 2;
        up = -1;
      else
        low = -1;
        mid = 2 - P;
        up = P - 1;
      end
    otherwise
      error('central and upwind are only supported');
  end

  rhs = zeros(N - 1, 1);
  rhs(1) = -low * phi_left;
  rhs(N - 1) = -up * phi_right;

  sub_diagonal = low * ones(N - 1, 1);
  diagonal = mid * ones(N - 1, 1);
  super_diagonal = up * ones(N - 1, 1);

  sol = [
    phi_left; 
    tdma(sub_diagonal, diagonal, super_diagonal, rhs, N - 1); 
    phi_right
  ];

  plot(linspace(0, L, N + 1), sol)

end