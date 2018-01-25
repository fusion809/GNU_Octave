function dy = simpen(y,t)
  dy    = zeros(2,1);
  g     = 9.8;
  l     = 1;
  dy(1) = y(2);
  dy(2) = -g/l*cos(y(1));