function [r] = Rrs2rrs(R)
  r = R ./ (0.52+1.7.*R);
end