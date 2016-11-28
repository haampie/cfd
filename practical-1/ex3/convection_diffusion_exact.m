function yexact = convection_diffusion_exact(x, k)

  pe = 1 / k;

  if pe > 100
     empe = 0;
  else
     empe = exp(-pe);
  end

  rl = 5;
  x1 = 0.4;

  a1 = -0.5 * rl;
  b1 = (0.4 * rl + 2 * k * a1);
  f1 = a1 * x1 * x1 + b1 * x1 - (2 * x1 * a1 + b1) / pe;
  g1 = (2 * a1 * x1 + b1) / pe;
  d2 = (1 - g1 * exp( - pe * x1) - f1) / (1 - empe);
  c2 = 1 - d2;
  c1 = c2 - f1;
  d1 = d2 * exp(pe * (x1 - 1)) - g1;


   if x < x1
      yexact = x * (a1 * x + b1) + c1 + d1 * exp(pe * (x - x1));
   else
      yexact = c2 + d2 * exp(pe * (x - 1));
   end
end