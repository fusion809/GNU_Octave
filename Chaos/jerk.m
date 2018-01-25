function dx = jerk(x,t)
  dx    = zeros(3,1);
  A     = 3/5;
  dx(1) = x(2);
  dx(2) = x(3);
  dx(3) = abs(x(1)) - 1 - A * x(3) - x(2);
