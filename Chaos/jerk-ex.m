clear all
N    = 1000000;
tmax = 1000;
t    = linspace(0, tmax, N+1)';
lsode_options("relative tolerance", 1e-15);
# With x0, dx0, d2x0 = 1 the solution is non-chaotic
# while with x0, dx0, d2x0 = 0 the solution is chaotic
x    = lsode("jerk",  [ 0.0; 0.0; 0.0]', t);

figure(1)

subplot(221)
plot(x(:,1),x(:,2),'-');
xlabel('x(t)','FontSize',16);
ylabel('dx(t)','FontSize',16);

subplot(222)
plot(x(:,1),x(:,3),'-');
xlabel('x(t)','FontSize',16);
ylabel('d2x(t)','FontSize',16);

subplot(223)
plot(x(:,2),x(:,3),'-');
xlabel('dx(t)','FontSize',16);
ylabel('d2x(t)','FontSize',16);

figure(2)
plot3(x(:,1),x(:,2),x(:,3))
xlabel('x(t)','FontSize',16);
ylabel('dx(t)','FontSize',16);
zlabel('d2x(t)','FontSize',16);
