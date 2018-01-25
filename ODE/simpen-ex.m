clear all
t=linspace(0,10,300001); 
lsode_options("relative tolerance", 1e-16);
y = lsode("simpen",[ 0; 0]',t);

figure(1)

subplot(221)
plot(t,y(:,1),'-');
title('Angle from x axis vs. time');
xlabel('t','FontSize',16);
ylabel('\theta(t)','FontSize',16);

subplot(222)
plot(t,y(:,2),'-');
title('Rate of angle change, from the x axis, vs. time');
xlabel('t','FontSize',16);
ylabel('d\theta(t)/dt','FontSize',16);

subplot(223)
plot(y(:,1),y(:,2),'-');
title('Phase space plot');
xlabel('\theta(t)','FontSize',16);
ylabel('d\theta(t)/dt','FontSize',16);