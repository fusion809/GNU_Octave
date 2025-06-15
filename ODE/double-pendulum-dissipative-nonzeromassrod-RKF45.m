clear all;
params     = struct();
params.g   = 9.81;
params.r1  = 2;
params.r2  = 2;
params.m1r = 9;
params.m2r = 9;
params.m1b = 1;
params.m2b = 1;
params.b1r = 0.10;
params.b1b = 0.10;
params.b2r = 0.10;
params.b2b = 0.10;
params.c1r = 0.04;
params.c2r = 0.04;
params.c1b = 0.04;
params.c2b = 0.04;
delay      = 1;
t0         = 0;
tf         = 20;
theta10    = 0;
theta20    = 0;
dtheta10   = 0;
dtheta20   = 0;
dt         = 1e-3;
t          = [t0];
epsil      = 1e-12;
X0         = [theta10 theta20 dtheta10 dtheta20];
i          = 1;
rhs        = @(y, t) DP2(params, y, t);
sol        = RKF45(rhs, [t0 tf], X0, dt, epsil);
t          = sol(:,1);
dt         = sol(:,2);
X          = sol(:,3:end);
theta1     = X(:,1);
theta2     = X(:,3);
dtheta1    = X(:,2);
dtheta2    = X(:,4);
% print values relating to the solution
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

% Descriptor
descriptor = ["rod masses = ", num2str(params.m1r), ", bob masses = "];
descriptor = [descriptor, num2str(params.m1b), ", bcoefs = "];
descriptor = [descriptor, num2str(params.b1r), ", ccoefs = "];
descriptor = [descriptor, num2str(params.c1r), ", tf = ", num2str(tf)];
%% Figure no is 5
animate(params, [params.r1 params.r2], [theta1 theta2], t, dt, 5, delay, descriptor);