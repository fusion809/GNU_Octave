clear all;
N          = 10000;
g          = 9.81;
l0         = 1;
k          = 1;
b          = 1;
c          = 1;
params     = struct('g', g, 'l0', l0, "k", k, "b", b, "c", c);
t0         = 0;
tf         = 10;
theta0     = 0;
z0         = 0;
dtheta0    = 0;
dz0        = 0;
dt         = [(tf-t0)/N];
t          = [t0];
epsil      = 1e-8;
theta      = [theta0];
z          = [z0];
dtheta     = [dtheta0];
dz         = [dz0];
i          = 1;
t          = linspace(t0, tf, N+1);
EPtht      = @(th, t) EP(params, th, t);
sol        = lsode(EPtht, [theta0 z0 dtheta0 dz0], t);
theta      = sol(:,1);
z          = sol(:,2);
dtheta     = sol(:,3);
dz         = sol(:,4);
printf("Number of t values used in the analysis = %d\n", length(t))
printf("Minimum theta = %d\n", min(theta))
printf("Minimum z = %d\n", min(z))
figure(1)
plot(t,theta)
title("Elastic pendulum: \\theta vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("\\theta", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(theta) max(theta)]);
print -dpng "Figure 1 Elastic pendulum theta vs t.png"
figure(2)
plot(t, dtheta)
title("Elastic pendulum: d\\theta/dt vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("d\\theta/dt", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 2/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(dtheta) max(dtheta)]);
print -dpng "Figure 2 Elastic pendulum dtheta vs t.png"
figure(3)
plot(t,z)
title("Elastic pendulum: z vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("z", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(z) max(z)]);
print -dpng "Figure 3 Elastic pendulum z vs t.png"
figure(4)
plot(t, dz)
title("Elastic pendulum: dz/dt vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
h=ylabel("dz/dt", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 2/10;      % Shift left by 1 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
xlim([t0 tf]);
ylim([min(dz) max(dz)]);
print -dpng "Figure 4 Elastic pendulum dz vs t.png"

% Animation
% x = (l0+z) .* cos(theta);
% y = (l0+z) .* sin(theta);

% % Set up figure
% figure(5); clf;
% axis equal off;
% axis([-2.2 2.2 -2.2 2.2]);
% hold on;

% % Output GIF filename
% gifname = "Figure 5 Elastic pendulum.gif";

% for i = 1:length(t)
%   % Clear and plot
%   cla;
%   plot([0 x(i)], [0 y(i)], 'b-', 'LineWidth', 2);      % Rod 1
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

% disp("GIF saved as 'Figure 5 Elastic pendulum.gif'");