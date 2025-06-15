clear all;
N         = 10000;
g         = 9.81;
r         = 1;
mr        = 1;
mb        = 1;
gamma     = 0.1;
c         = 0.01;
params    = [g; r; mr; mb; gamma; c];
t0        = 0;
tf        = 10;
theta0    = 0;
dtheta0   = 0;
dtInit    = (tf-t0)/N;
epsil     = 1e-16;
rhs       = @(y,t) sinpen(params, y, t);
sol       = RKF45(rhs, [t0 tf], [theta0 dtheta0], dtInit, epsil);
t         = sol(:,1);
dt        = sol(:,2);
theta     = sol(:,3);
dtheta    = sol(:,4);
printf("Number of t values used in the analysis = %d\n", length(t))
printf("Minimum theta = %d\n", min(theta))
figure(1)
plot(t,theta)
title("Single pendulum: \\theta vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("\\theta", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 0/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(theta) max(theta)]);
print -dpng "Figure 1 Single pendulum theta vs t.png"
figure(2)
plot(t, dtheta)
title("Single pendulum: d\\theta/dt vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("d\\theta/dt", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1.5/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(dtheta) max(dtheta)]);
print -dpng "Figure 2 Single pendulum dtheta vs t.png"
% Animation
x = r * cos(theta);
y = r * sin(theta);

% Set up figure
% figure(3); clf;
% axis equal off;
% axis([-2.2 2.2 -2.2 2.2]);
% hold on;

% % Output GIF filename
% gifname = "Figure 3 Single pendulum.gif";

% for i = 1:length(t)
%   % Clear and plot
%   cla;
%   plot([0 x(i)], [0 y(i)], 'b-', 'LineWidth', 2);      % Rod 1
%   plot(x(i), y(i), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
%   title(sprintf("t = %.2f s", t(i)));

%   drawnow;
%   frame = getframe(gcf);
%   img = frame2im(frame);
%   [A, map] = rgb2ind(img);

%   if i == 1
%     imwrite(A, map, gifname, "gif", "LoopCount", Inf, "DelayTime", dt(i));
%   else
%     imwrite(A, map, gifname, "gif", "WriteMode", "append", "DelayTime", dt(i));
%   end
% end

% disp("GIF saved as 'Figure 3 Single pendulum.gif'");