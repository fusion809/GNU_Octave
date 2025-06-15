clear all;
N         = 10000;
g         = 9.81;
l         = 1;
params    = [g l];
k         = 1:100;
period    = 4*sqrt(l/g)*(pi/2 + 2*pi*sum(exp(k*pi)./(1+exp(2*k*pi))));
t0        = 0;
tf        = period;
theta0    = 0;
dtheta0   = 0;
dtInit    = (tf-t0)/N;
epsil     = 1e-16;
rhs       = @(y,t) simpen(params, y, t);
sol       = RKF45(rhs, [t0 tf], [theta0 dtheta0], dtInit, epsil);
t         = sol(:,1);
dt        = sol(:,2);
theta     = sol(:,3);
dtheta    = sol(:,4);
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
xlim([t0 tf]);
ylim([min(theta) max(theta)]);
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
xlim([t0 tf]);
ylim([min(dtheta) max(dtheta)]);
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