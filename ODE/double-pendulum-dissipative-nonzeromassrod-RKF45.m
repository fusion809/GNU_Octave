clear all;
N          = 10000;
g          = 9.81;
r1         = 1;
r2         = 1;
m1r        = 1;
m2r        = 1;
m1b        = 1;
m2b        = 1;
b1r        = 0.10;
b1b        = 0.10;
b2r        = 0.10;
b2b        = 0.10;
c1r        = 0.04;
c2r        = 0.04;
c1b        = 0.04;
c2b        = 0.04;
delay      = 1;
params     = struct("g", g, "r1", r1, "r2", r2, "m1r", m1r, "m1b", m1b, ...
"m2r", m2r, "m2b", m2b, "b1r", b1r, "b1b", b1b, "b2r", b2r, "b2b", b2b, ...
"c1r", c1r, "c1b", c1b, "c2r", c2r, "c2b", c2b);
t0         = 0;
tf         = 10;
descriptor = ["masses = ", num2str(m1r), ", bcoefs = ", num2str(b1r), ", ccoefs = ", num2str(c1r), ", tf = ", num2str(tf)];
theta10    = 0;
theta20    = 0;
dtheta10   = 0;
dtheta20   = 0;
dt         = [(tf-t0)/N];
t          = [t0];
epsil      = 1e-7;
theta1     = [theta10];
theta2     = [theta20];
dtheta1    = [dtheta10];
dtheta2    = [dtheta20];
i          = 1;
while t(i)<tf;
    K1 = dt(i)*DP2(params, [theta1(i); dtheta1(i); theta2(i); dtheta2(i)], t(i));
    K2 = dt(i)*DP2(params, [theta1(i); dtheta1(i); theta2(i); dtheta2(i)] + 1/4*K1, t(i)+1/4*dt(i));
    K3 = dt(i)*DP2(params, [theta1(i); dtheta1(i); theta2(i); dtheta2(i)] + 3/32*K1+9/32*K2, t(i)+3/8*dt(i));
    K4 = dt(i)*DP2(params, [theta1(i); dtheta1(i); theta2(i); dtheta2(i)] + 1932/2197*K1 - 7200/2197*K2 + 7296/2197*K3, t(i) + 12/13*dt(i));
    K5 = dt(i)*DP2(params, [theta1(i); dtheta1(i); theta2(i); dtheta2(i)] + 439/216*K1 - 8*K2 + 3680/513*K3 - 845/4104*K4, t(i)+dt(i));
    K6 = dt(i)*DP2(params, [theta1(i); dtheta1(i); theta2(i); dtheta2(i)] - 8/27*K1 + 2*K2 - 3544/2565*K3 + 1859/4104*K4 - 11/40*K5, t(i)+dt(i));

    theta1_X1  = theta1(i) + 25/216 * K1(1) + 1408/2565*K3(1) + 2197/4104*K4(1) - 1/5*K5(1);
    dtheta1_X1 = dtheta1(i) + 25/216 * K1(2) + 1408/2565*K3(2) + 2197/4104*K4(2) - 1/5*K5(2);
    theta1_X2  = theta1(i) + 16/135 * K1(1) + 6656/12825*K3(1) + 28561/56430*K4(1) - 9/50 * K5(1) + 2/55 * K6(1);
    dtheta1_X2 = dtheta1(i) + 16/135 * K1(2) + 6656/12825*K3(2) + 28561/56430*K4(2) - 9/50 * K5(2) + 2/55 * K6(2);
    theta2_X1  = theta2(i) + 25/216 * K1(3) + 1408/2565*K3(3) + 2197/4104*K4(3) - 1/5*K5(3);
    dtheta2_X1 = dtheta2(i) + 25/216 * K1(4) + 1408/2565*K3(4) + 2197/4104*K4(4) - 1/5*K5(4);
    theta2_X2  = theta2(i) + 16/135 * K1(3) + 6656/12825*K3(3) + 28561/56430*K4(3) - 9/50 * K5(3) + 2/55 * K6(3);
    dtheta2_X2 = dtheta2(i) + 16/135 * K1(4) + 6656/12825*K3(4) + 28561/56430*K4(4) - 9/50 * K5(4) + 2/55 * K6(4);
    theta1   = [theta1; theta1_X1];
    theta2   = [theta2; theta2_X1];
    dtheta1  = [dtheta1; dtheta1_X1];
    dtheta2  = [dtheta2; dtheta2_X1];
    TE = max(abs(-1/360 * K1 + 128/4275 * K3 + 2197/75240 * K4 - 1/50 * K5 - 2/55 * K6));
    s = 0.9*(epsil/TE)^(1/5);
    if (s*dt(i)+t(i) < tf)
        dt = [dt; s*dt(i)];
    else
        dt = [dt; tf-t(i)];
    end
    t = [t; t(i)+dt(i)];
    i = i + 1;
end
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
gifname = ["Figure 5 Double pendulum ", num2str(delay), " delay ", descriptor, ".gif"];

frame(length(1:length(t))) = struct('cdata', [], 'colormap', []);
k = 1;
for i = 1:length(t)
    cla;
  plot([0 x1(i)], [0 y1(i)], 'b-', 'LineWidth', 2);      % Rod 1
  plot([x1(i) x2(i)], [y1(i) y2(i)], 'r-', 'LineWidth', 2);  % Rod 2
  plot(x1(i), y1(i), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
  plot(x2(i), y2(i), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
  title(sprintf("t = %.2f s", t(i)));

  drawnow;
  frame(k) = getframe(gcf);
  k = k+1;
end


for k=1:delay:length(frame)
  img = frame2im(frame(k));
  [A, map] = rgb2ind(img);

  if k == 1
    imwrite(A, map, gifname, "gif", "LoopCount", Inf, "DelayTime", dt(k));
  else
    imwrite(A, map, gifname, "gif", "WriteMode", "append", "DelayTime", dt(k));
  end
end

disp(["GIF saved as ", gifname]);