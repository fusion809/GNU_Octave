clear all
N         = 1e3;
NN        = 1e5;
Iter      = 100;
b         = 3;
g         = 9.8;
l         = 1;
theta0    = 0;
dtheta0   = 0;
n         = 0:N;
parat     = pi+n'*pi/N;
x         = cos(parat);
t         = b/2*(x+1);
T         = cos(parat*n);
# Now for arrays that do not include the endpoints
xsub      = x(2:N);
paratsub  = parat(2:N);
Tsub      = T(2:N,:);
Usub      = diag(1./sqrt(1-xsub.^2))*sin(paratsub*n);
dTsub     = Usub*diag(n);
dT        = [-(n.^2).*(-1).^(n); dTsub ; (n).^2];
D1        = 2/b*dT/T;
D2        = D1^2;
thetainit = -pi/2+pi/2*cos(pi/1.184*t);
theta     = zeros(N+1,Iter); 
theta(:,1)= thetainit;
a(:,1)    = T\thetainit;

for i=2:Iter
  dtheta=D1*theta(:,i-1);
  RHS=-D2*theta(:,i-1)-g/l*cos(theta(:,i-1));
  LHS=D2-g/l*diag(sin(theta(:,i-1)));
  LHS(1,1)=1;
  LHS(1,2:N+1)=0;
  LHS(N+1,:)=D1(1,:);
  RHS(1)=theta0-theta(1,i-1);
  RHS(N+1)=dtheta0-dtheta(1);
  delta(:,i-1)=LHS\RHS;
  theta(:,i)=theta(:,i-1)+delta(:,i-1);
  a(:,i)=T\theta(:,i);
  deltarms(i-1)=sqrt(delta(:,i-1)'*delta(:,i-1)/(N+1));
  clear RHS LHS dtheta;
endfor

aend=a(:,end);
xx=linspace(-1,1,NN+1)';
tt=b/2*(xx+1);
TT=cos(acos(xx)*n);
thetalin=TT*a;
figure(1)
plot(t,thetainit);
figure(2)
plot(tt,thetalin(:,end))