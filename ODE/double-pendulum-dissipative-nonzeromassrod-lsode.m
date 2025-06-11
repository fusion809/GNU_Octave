% Angle measured relative to negative y-axis based on ChatGPT code. Not sure if the dissipative forces are counted right.
params = struct('m1', 1.0, 'm2', 1.0, ...
                'M1', 1.0, 'M2', 1.0, ...
                'L1', 1.0, 'L2', 1.0, ...
                'g', 9.81, ...
                'gamma1', -0.1, 'gamma2', -0.1, ...
                'c1', -0.1, 'c2', -0.1);

t0 = 0;
tf = 100;
theta0 = [pi/2; pi/2; 0; 0];  % initial angles and angular velocities
t = linspace(t0, tf, 10000);    % time vector
rhs = @(y, t) double_pendulum_rhs(y, t, params);
sol = lsode(rhs, theta0, t);
figure(1);
plot(t, sol(:,1));
title("Double complex pendulum: \\theta_1 vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("\\theta_1", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(sol(:,1)) max(sol(:,1))]);
print -dpng "Figure 1 Double complex pendulum theta1 vs t.png"
figure(2);
plot(t, sol(:,2));
title("Double complex pendulum: \\theta_2 vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("\\theta_2", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(sol(:,2)) max(sol(:,2))]);
print -dpng "Figure 1 Double complex pendulum theta2 vs t.png"
figure(3);
plot(t, sol(:,3));
title("Double complex pendulum: d\\theta_1/dt vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("d\\theta_1/dt", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(sol(:,3)) max(sol(:,3))]);
print -dpng "Figure 3 Double complex pendulum dtheta1 vs t.png"
figure(4);
plot(t, sol(:,4));
title("Double complex pendulum: d\\theta_2/dt vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("d\\theta_2/dt", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(sol(:,4)) max(sol(:,4))]);
print -dpng "Figure 4 Double complex pendulum dtheta2 vs t.png"