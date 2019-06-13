clear all
N=2000000;
x0=1;
y0=1;
z0=1;
tend=2000;
t=linspace(0,tend,N+1)';
X=lsode("lorenz",[x0; y0; z0]',t);
figure(1)
plot3(X(:,1),X(:,2),X(:,3))
xlabel('x(t)','FontSize',18);
ylabel('y(t)','FontSize',18);
zlabel('z(t)','FontSize',18);

figure(2)
subplot(221)
plot(X(:,1),X(:,2),'-');
xlabel('x(t)','FontSize',16);
ylabel('y(t)','FontSize',16);

subplot(222)
plot(X(:,1),X(:,3),'-');
xlabel('x(t)','FontSize',16);
ylabel('z(t)','FontSize',16);

subplot(223)
plot(X(:,2),X(:,3),'-');
xlabel('y(t)','FontSize',16);
ylabel('z(t)','FontSize',16);