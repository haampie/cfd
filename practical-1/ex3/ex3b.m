function ex3b()
  close all;

  U = 50;
  L = 1;
  k = 0.1;
  Pe = U * L / k;

  phi_0 = 0;
  phi_L = 20;

  N_values = 2 .^ (2 : 20);
  h_values = L ./ N_values;
  errors = zeros(2, length(N_values) - 3);

  for idx = 1 : length(N_values) - 3
    [x1, y1] = convection_diffusion_dirichlet(U, L, k, N_values(idx + 0), phi_0, phi_L, 'upwind');
    [x2, y2] = convection_diffusion_dirichlet(U, L, k, N_values(idx + 1), phi_0, phi_L, 'upwind');
    [x3, y3] = convection_diffusion_dirichlet(U, L, k, N_values(idx + 2), phi_0, phi_L, 'upwind');

    [~, y1_C] = convection_diffusion_dirichlet(U, L, k, N_values(idx + 0), phi_0, phi_L, 'central');
    [~, y2_C] = convection_diffusion_dirichlet(U, L, k, N_values(idx + 1), phi_0, phi_L, 'central');
    [~, y3_C] = convection_diffusion_dirichlet(U, L, k, N_values(idx + 2), phi_0, phi_L, 'central');

    h1 = h_values(idx);
    h2 = h_values(idx + 1);
    h3 = h_values(idx + 2);

    errors(1, idx) = log(l2_norm(y1 - y2(1 : 2 : end), h1) / l2_norm(y2 - y3(1 : 2 : end), h2)) / log(h1 / h2);
    errors(2, idx) = log(l2_norm(y1_C - y2_C(1 : 2 : end), h1) / l2_norm(y2_C - y3_C(1 : 2 : end), h2)) / log(h1 / h2);
  end

  plot(errors(1, :));
  hold on;
  plot(errors(2, :));
  legend('upwind', 'central');
  ylabel('rate of convergence')
  xlabel('steps')
end

function nx = l2_norm(x, h)
  nx = norm(x) * sqrt(h);
end