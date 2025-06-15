clear all;
N          = 10000;
g          = 9.81;
r1         = 1;
r2         = 1;
m1         = 1;
m2         = 1;
params     = [g; r1; r2; m1; m2];
t0         = 0;
tf         = 10;
theta10    = 0;
theta20    = 0;
dtheta10   = 0;
dtheta20   = 0;
dt         = [(tf-t0)/N];
t          = [t0];
epsil      = 1e-10;
i          = 1;
rhs        = @(y, t) DP(params, y, t);
sol        = RKF45(rhs, [t0 tf], [theta10 theta20 dtheta10 dtheta20], (tf-t0)/N, epsil);
t          = sol(:,1);
dt         = sol(:,2);
X          = sol(:,3:end);
theta1     = X(:,1);
theta2     = X(:,2);
dtheta1    = X(:,3);
dtheta2    = X(:,4);
printf("Number of t values used in the analysis = %d\n", length(t))
printf("Minimum theta1 = %d\n", min(theta1))
printf("Minimum theta2 = %d\n", min(theta2))
figure(1)
plot(t,theta1)
title("Double pendulum: \\theta_1 vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("\\theta_1", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(theta1) max(theta1)]);
print -dpng "Figure 1 Double pendulum theta1 vs t.png"
figure(2)
plot(t, dtheta1)
title("Double pendulum: d\\theta_1/dt vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("d\\theta_1/dt", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 2/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(dtheta1) max(dtheta1)]);
print -dpng "Figure 2 Double pendulum dtheta1 vs t.png"
figure(3)
plot(t,theta2)
title("Double pendulum: \\theta_2 vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("\\theta_2", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(theta2) max(theta2)]);
print -dpng "Figure 3 Double pendulum theta2 vs t.png"
figure(4)
plot(t, dtheta2)
title("Double pendulum: d\\theta_2/dt vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("d\\theta_2/dt", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 2/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(dtheta2) max(dtheta2)]);
print -dpng "Figure 4 Double pendulum dtheta2 vs t.png"

% Animation
x1 = r1 * cos(theta1);
y1 = r1 * sin(theta1);
x2 = x1 + r2 * cos(theta2);
y2 = y1 + r2 * sin(theta2);

% Set up figure
figure(5); clf;
axis equal off;
axis([-2.2 2.2 -2.2 2.2]);
hold on;

% Output GIF filename
gifname = "Figure 5 Double pendulum.gif";

for i = 1:length(t)
  % Clear and plot
  cla;
  plot([0 x1(i)], [0 y1(i)], 'b-', 'LineWidth', 2);      % Rod 1
  plot([x1(i) x2(i)], [y1(i) y2(i)], 'r-', 'LineWidth', 2);  % Rod 2
  plot(x1(i), y1(i), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
  plot(x2(i), y2(i), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
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

disp("GIF saved as 'Figure 5 Double pendulum.gif'");