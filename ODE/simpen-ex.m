clear all
t=linspace(0,5,1000000001); 
lsode_options("absolute tolerance", 1e-30);
y = lsode("simpen",[ 0; 0]',t);
theta = y(:,1);
dtheta = y(:,2);

figure(1)

subplot(221)
plot(t,theta,'-');
title('Angle from x axis vs. time');
xlabel('t','FontSize',16);
ylabel('\theta(t)','FontSize',16);

subplot(222)
plot(t,dtheta,'-');
title('Rate of angle change, from the x axis, vs. time');
xlabel('t','FontSize',16);
ylabel('d\theta(t)/dt','FontSize',16);

subplot(223)
plot(theta,dtheta,'-');
title('Phase space plot');
xlabel('\theta(t)','FontSize',16);
ylabel('d\theta(t)/dt','FontSize',16);