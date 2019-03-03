N     = 10000000;
t0    = 0;
t1    = 100;
x0    = -0.1;
dx0   = 0;

function Der = f0(X,t)
	alpha = 0.02;
	f     = 10;
	w     = 1;
	Der(1) = X(2);
	Der(2) = - alpha * X(2) + 0.5*X(1)*(1-(X(1))^2) + f*cos(w*t);
endfunction


t  = linspace(t0, t1, N+1)';
X0 = [x0; dx0];
X  = lsode("f0", X0, t);

figure(1)
plot(X(:,1), X(:,2), 'k');
figure(2)
plot(t, X(:,1), 'r');
figure(3)
plot(t, X(:,2), 'g');
