% Parameters: [m1, m2, g, k1, k2, l1, l2]
params = [1, 1, 9.81, 10, 10, 1, 1];

% Initial conditions
y0 = [1.2; 0; pi/2; 0; 1.0; 0; pi/2 + 0.1; 0];

% Time span
tspan = [0 20];

% Solve
[t, y] = ode45(@(t, y) double_elastic_pendulum(t, y, params), tspan, y0);

% Plot result
plot(t, y(:,1), 'r', t, y(:,5), 'b');
xlabel('Time (s)');
ylabel('r_1 and r_2');
legend('r_1', 'r_2');
title('Elastic Double Pendulum');