g = 9.81;
r = 1;
k = 1:100;
N = 100000;
theta0 = 0;
thetadot0 = 0;
x0 = [theta0; thetadot0];
period = 4*sqrt(r/g)*(pi/2 + 2*pi*sum(exp(k*pi)./(1+exp(2*k*pi))));
t = linspace(0, period, N+1);
lsode_options(1)=1e-20;
y = lsode(@f, x0, t);
y = [y -g/r*cos(y(:,1))];
figure(1);
plot(t, y(:,1), 'r');
title("theta");
figure(2);
plot(t, y(:,2), 'g');
title("dtheta/dt");