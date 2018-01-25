clear all
t=linspace(0,350,100001); 
lsode_options("relative tolerance", 1e-15);
y = lsode("lorenz",[ 1.0; 1.0; 1.0]',t);

figure(1)

subplot(221)
plot(y(:,1),y(:,2),'-');
xlabel('x(t)','FontSize',16);
ylabel('y(t)','FontSize',16);

subplot(222)
plot(y(:,1),y(:,3),'-');
xlabel('x(t)','FontSize',16);
ylabel('z(t)','FontSize',16);

subplot(223)
plot(y(:,2),y(:,3),'-');
xlabel('y(t)','FontSize',16);
ylabel('z(t)','FontSize',16);