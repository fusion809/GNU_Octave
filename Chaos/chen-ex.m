clear all
NN=1e6;
tf=1e2;
t=linspace(0,tf,NN+1); 
lsode_options("relative tolerance", 1e-15);
X = lsode("chen",[ -0.1; 0.5; -0.6]',t);
x = X(:,1); y=X(:,2); z=X(:,3);

figure(1)
plot3(x,y,z)
xlabel('x(t)','FontSize',18);
ylabel('y(t)','FontSize',18);
zlabel('z(t)','FontSize',18);

figure(2)
subplot(221)
plot(x,y,'-');
xlabel('x(t)','FontSize',16);
ylabel('y(t)','FontSize',16);

subplot(222)
plot(x,z,'-');
xlabel('x(t)','FontSize',16);
ylabel('z(t)','FontSize',16);

subplot(223)
plot(y,z,'-');
xlabel('y(t)','FontSize',16);
ylabel('z(t)','FontSize',16);
