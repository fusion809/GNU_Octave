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
while t(i)<tf;
    K1 = dt(i)*EP(params, [theta(i); z(i); dtheta(i); dz(i)], t(i));
    K2 = dt(i)*EP(params, [theta(i); z(i); dtheta(i); dz(i)] + 1/4*K1, t(i)+1/4*dt(i));
    K3 = dt(i)*EP(params, [theta(i); z(i); dtheta(i); dz(i)] + 3/32*K1+9/32*K2, t(i)+3/8*dt(i));
    K4 = dt(i)*EP(params, [theta(i); z(i); dtheta(i); dz(i)] + 1932/2197*K1 - 7200/2197*K2 + 7296/2197*K3, t(i) + 12/13*dt(i));
    K5 = dt(i)*EP(params, [theta(i); z(i); dtheta(i); dz(i)] + 439/216*K1 - 8*K2 + 3680/513*K3 - 845/4104*K4, t(i)+dt(i));
    K6 = dt(i)*EP(params, [theta(i); z(i); dtheta(i); dz(i)] - 8/27*K1 + 2*K2 - 3544/2565*K3 + 1859/4104*K4 - 11/40*K5, t(i)+dt(i));

    theta_X1  = theta(i) + 25/216 * K1(1) + 1408/2565*K3(1) + 2197/4104*K4(1) - 1/5*K5(1);
    dtheta_X1 = dtheta(i) + 25/216 * K1(3) + 1408/2565*K3(3) + 2197/4104*K4(3) - 1/5*K5(3);
    theta_X2  = theta(i) + 16/135 * K1(1) + 6656/12825*K3(1) + 28561/56430*K4(1) - 9/50 * K5(1) + 2/55 * K6(1);
    dtheta_X2 = dtheta(i) + 16/135 * K1(3) + 6656/12825*K3(3) + 28561/56430*K4(3) - 9/50 * K5(3) + 2/55 * K6(3);
    z_X1      = z(i) + 25/216 * K1(2) + 1408/2565*K3(2) + 2197/4104*K4(2) - 1/5*K5(2);
    dz_X1     = dz(i) + 25/216 * K1(4) + 1408/2565*K3(4) + 2197/4104*K4(4) - 1/5*K5(4);
    z_X2      = z(i) + 16/135 * K1(2) + 6656/12825*K3(2) + 28561/56430*K4(2) - 9/50 * K5(2) + 2/55 * K6(2);
    dz_X2     = dz(i) + 16/135 * K1(4) + 6656/12825*K3(4) + 28561/56430*K4(4) - 9/50 * K5(4) + 2/55 * K6(4);
    theta     = [theta; theta_X1];
    z         = [z; z_X1];
    dtheta    = [dtheta; dtheta_X1];
    dz        = [dz; dz_X1];
    TE        = max(abs(-1/360 * K1 + 128/4275 * K3 + 2197/75240 * K4 - 1/50 * K5 - 2/55 * K6));
    s         = 0.9*(epsil/TE)^(1/5);
    if (s*dt(i)+t(i) < tf)
        dt    = [dt; s*dt(i)];
    else
        dt    = [dt; tf-t(i)];
    end
    t = [t; t(i)+dt(i)];
    i = i + 1;
end
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