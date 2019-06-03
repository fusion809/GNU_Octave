clear all
N = 10;
n = 0:N;
x = -cos(pi*n/N)';
xsub = x(2:N);
T = cos(acos(x)*n);
Tsub = T(2:N,2:N);
y = log(1+x);
ysub = y(2:N);

a = Tsub\ysub;