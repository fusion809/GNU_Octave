clear all
N=100000;
y0=0;
dy0=1;
tend=1000;
x=linspace(0,tend,N+1)';
y=lsode("f",[y0; dy0]',x);
figure(1)
plot(x,y)
figure(2)
plot(y(:,1),y(:,2))