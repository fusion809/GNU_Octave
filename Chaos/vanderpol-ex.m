clear all
NN=1e6;
tf=3e3;
t=linspace(0,tf,NN+1); 
lsode_options("relative tolerance", 1e-15);
x = lsode("vanderpol",[ 0; 1.0]',t);

figure(1)
plot(x(:,1),x(:,2),'-');
xlabel('x(t)','FontSize',16);
ylabel('dx(t)','FontSize',16);