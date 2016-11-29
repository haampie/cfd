function ex2b()

  close all

  N_values = 10 : 5 : 1000;
  BL_values = zeros(size(N_values));
  BL_start = 0.1;

  figure;
  for idx = 1 : length(N_values)
    [x, y] = convection_diffusion_dirichlet(1000, 1, 1, N_values(idx), 0, 1, 'upwind');

    blnums = find(y > BL_start);
    blpos = blnums(1);

    x1 = x(blpos - 1);
    x2 = x(blpos);
    y1 = y(blpos - 1);
    y2 = y(blpos);

    slope = (y2 - y1) / (x2 - x1);
    bl_point = x1 + (BL_start - y1) / slope;
    BL_values(idx) = 1 - bl_point;

    plot(x, y); hold on;
  end


  figure;
  plot(1 ./ N_values, BL_values, '+');
  legend('experimental')
  xlabel('h-values')
  ylabel('BL-thickness')
end