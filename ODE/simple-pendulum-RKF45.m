clear all;
N         = 10000;
g         = 9.81;
l         = 1;
k         = 1:100;
period    = 4*sqrt(l/g)*(pi/2 + 2*pi*sum(exp(k*pi)./(1+exp(2*k*pi))));
t0        = 0;
tf        = period;
theta0    = 0;
dtheta0   = 0;
dt        = [(tf-t0)/N];
t         = [t0];
epsil     = 1e-16;
theta     = [theta0];
dtheta    = [dtheta0];
dtheta(1) = dtheta0;
i         = 1;
while t(i)<tf;
    K1 = dt(i)*simpen([theta(i); dtheta(i)], t(i));
    K2 = dt(i)*simpen([theta(i); dtheta(i)] + 1/4*K1, t(i)+1/4*dt(i));
    K3 = dt(i)*simpen([theta(i); dtheta(i)] + 3/32*K1+9/32*K2, t(i)+3/8*dt(i));
    K4 = dt(i)*simpen([theta(i); dtheta(i)] + 1932/2197*K1 - 7200/2197*K2 + 7296/2197*K3, t(i) + 12/13*dt(i));
    K5 = dt(i)*simpen([theta(i); dtheta(i)] + 439/216*K1 - 8*K2 + 3680/513*K3 - 845/4104*K4, t(i)+dt(i));
    K6 = dt(i)*simpen([theta(i); dtheta(i)] - 8/27*K1 + 2*K2 - 3544/2565*K3 + 1859/4104*K4 - 11/40*K5, t(i)+dt(i));

    theta1  = theta(i) + 25/216 * K1(1) + 1408/2565*K3(1) + 2197/4104*K4(1) - 1/5*K5(1);
    dtheta1 = dtheta(i) + 25/216 * K1(2) + 1408/2565*K3(2) + 2197/4104*K4(2) - 1/5*K5(2);
    theta2  = theta(i) + 16/135 * K1(1) + 6656/12825*K3(1) + 28561/56430*K4(1) - 9/50 * K5(1) + 2/55 * K6(1);
    dtheta2 = dtheta(i) + 16/135 * K1(2) + 6656/12825*K3(2) + 28561/56430*K4(2) - 9/50 * K5(2) + 2/55 * K6(2);
    theta   = [theta; theta1];
    dtheta  = [dtheta; dtheta1];
    TE = max(abs(-1/360 * K1 + 128/4275 * K3 + 2197/75240 * K4 - 1/50 * K5 - 2/55 * K6));
    s = 0.9*(epsil/TE)^(1/5);
    if (s*dt(i)+t(i) < tf)
        dt = [dt; s*dt(i)];
    else
        dt = [dt; tf - t(i)];
    end
    t = [t; t(i)+dt(i)];
    i = i + 1;
end
printf("Number of t values used in the analysis = %d\n", length(t))
printf("Minimum theta = %d\n", min(theta))
figure(1)
plot(t,theta)
title("Simple pendulum: \\theta vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("\\theta", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 0/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([0 period]);
ylim([-pi 0]);
print -dpng "Figure 1 Simple pendulum theta vs t.png"
figure(2)
plot(t, dtheta)
title("Simple pendulum: d\\theta/dt vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("d\\theta/dt", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1.5/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([0 period]);
ylim([-sqrt(2*g/l) sqrt(2*g/l)]);
print -dpng "Figure 2 Simple pendulum dtheta vs t.png"

% Calculate residuals and plot them; commented out 
% NN = length(t)-1;
% tch = (tf-t0)/2*(cos(pi*(0:length(t)-1)'/(length(t)-1))+1);
% dthetach = spline(t, dtheta, tch);
% d2thetach = 2/(tf-t0) * cheb(NN)*dthetach;
% d2theta = spline(tch, d2thetach, t);
% residual = abs(d2theta + g/l*cos(theta));
% figure(3)
% semilogy(t, residual)
% title("Simple pendulum: residual vs t", "fontsize", 25)
% xlabel("t", "fontsize", 20);
% h=ylabel("Residual", "fontsize", 20, "rotation", 0);
% pos = get(h, 'position');            % Get current [x, y, z] position
% new_pos = pos;                       
% new_pos(1) = new_pos(1) - 1.5/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
% set(h, 'position', new_pos);         % Apply new position
% set(gca, 'fontsize', 16)
% xlim([0 period]);
% ylim([min(residual) max(residual)]);
% print -dpng "Figure 3 Simple pendulum residual vs t.png"
% 1st derivative residual
residual = abs(dtheta - sign(dtheta).*sqrt(dtheta0^2 + 2*g/l*(sin(theta0)-sin(theta))));
figure(3)
semilogy(t, residual)
title("Simple pendulum: residual vs t", "fontsize", 25)
xlabel("t", "fontsize", 20);
h=ylabel("Residual", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1.5/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([0 period]);
ylim([min(residual) max(residual)]);
print -dpng "Figure 3 Simple pendulum residual vs t.png"

% Animation
x = l * cos(theta);
y = l * sin(theta);

% Set up figure
figure(4); clf;
axis equal off;
axis([-2.2 2.2 -2.2 2.2]);
hold on;

% Output GIF filename
gifname = "Figure 4 Simple pendulum.gif";

for i = 1:length(t)
  % Clear and plot
  cla;
  plot([0 x(i)], [0 y(i)], 'b-', 'LineWidth', 2);      % Rod 1
  plot(x(i), y(i), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
  title(sprintf("t = %.2f s", t(i)));

  drawnow;
  frame = getframe(gcf);
  img = frame2im(frame);
  [A, map] = rgb2ind(img);

  if i == 1
    imwrite(A, map, gifname, "gif", "LoopCount", Inf, "DelayTime", dt(i));
  else
    imwrite(A, map, gifname, "gif", "WriteMode", "append", "DelayTime", dt(i));
  end
end

disp("GIF saved as 'Figure 4 Simple pendulum.gif'");