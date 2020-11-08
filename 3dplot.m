clear all
N=1000;
t=linspace(0,100,N+1)';
f=[1+3*t, -1+2*t, -4+5*t];
f2=[4+3*t, -3-2*t, 2-t];
plot3(f(:,1),f(:,2),f(:,3),"linewidth", 4)
hold on
plot3(f2(:,1),f2(:,2),f2(:,3),"linewidth", 4)