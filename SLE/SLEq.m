% Solves the time-independent Schrodinger equation for
% Potential energy V(x) = x
clear all
format long g
N=2000;
a=0; b=80;
n=0:N;
x=cos(pi*n'/N); xsub=x(2:N); % Set up chebyshev grid
Q=(1-xsub).^2;
y=(b-a)/2*x+(a+b)/2; % Set up standard grid
ysub=y(2:N);
T=cos(acos(x)*n); % Set up Chebyshev functions
Usub=diag(1./sqrt(1-xsub.^2))*sin(acos(xsub)*n);
dT=[((-1).^n).*n.^2; Usub*diag(n); n.^2];
D1=dT/T; D2=D1^2;
E1=D1(2:N,2:N);
E2=D2(2:N,2:N);
d1=2/(b-a)*E1;
d2=4/(b-a)^2*E2;
H=-d2+diag(ysub);
[Psisub, LAM]=eig(H);
Lam=diag(LAM);
[lam IX]=sort(Lam, 'ascend');
Psisub=Psisub(:,IX);
Psi=[zeros(1,N-1); Psisub; zeros(1,N-1)];
Amp=Psi.^2;
figure(1); plot(y,Amp(:,1)); title('Amplitude first vector')
figure(2); plot(y,Amp(:,10))
figure(3); plot(y,Amp(:,20))
figure(4); plot(y,Amp(:,40))
