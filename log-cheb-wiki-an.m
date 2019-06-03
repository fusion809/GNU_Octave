clear all
N = 20000;
t = linspace(-1,1,N+1)';
n = 0:N;
nsub = 1:N; nsub = nsub';
T = cos(acos(t)*n);
a = [-log(2); -2*((-1).^nsub)./nsub];
y = T*a;
plot(t,y)