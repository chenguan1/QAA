function [u] = rrs2u(r)
  g0 = 0.089;
  g1 = 0.1245;
  u = (-g0 + sqrt(g0*g0 + 4*g0*g1*r) ) / (2 * g1);
end