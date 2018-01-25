# Variables
N   = 10000;
x0  = 0;
x1  = 10;
Y0  = 0;
y0  = 1;
dy0 = 1;

function DY = f(y, x)
  
  C      = 1;
  lambda = 1;
  n      = 10;
  DY     = zeros (2,1);

  DY(1)  = y(2);
  DY(2)  = (C^2)/(y(1)^2)-lambda^2*(y(1))^(n-1);

endfunction

x = linspace(x0, x1, N+1)';
Y0 = [y0 dy0]';
Y = lsode("f", Y0, x);

figure(1)
plot(x,Y(:,1))

figure(2)
plot(x,Y(:,2))