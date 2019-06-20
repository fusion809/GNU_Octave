clear all
N         = 2e3;
NN        = 3e5;
Iter      = 5;
b         = 10;
period    = 1.1845;
g         = 9.8;
l         = 1;
theta0    = 0;
dtheta0   = 0;
n         = 0:N;
parat     = pi+n'*pi/N;
Chx       = cos(parat);
t         = b/2*(Chx+1);
T         = cos(parat*n);
# Now for arrays that do not include the endpoints
Chxsub    = Chx(2:N);
paratsub  = parat(2:N);
Tsub      = T(2:N,:);
Usub      = diag(1./sqrt(1-Chxsub.^2))*sin(paratsub*n);
dTsub     = Usub*diag(n);
dT        = [-(n.^2).*(-1).^(n); dTsub ; (n).^2];
D1        = 2/b*dT/T;
D2        = D1^2;
thetainit = -pi/2+pi/2*cos(pi/period*t);
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
  rcon(i-1)=rcond(LHS);
  # Initial conditions
  # designed to reverse any creeking in init. conds.
  RHS(1)=theta0-theta(1,i-1);
  RHS(N+1)=dtheta0-dtheta(1);
  delta(:,i-1)=LHS\RHS;
  theta(:,i)=theta(:,i-1)+delta(:,i-1);
  a(:,i)=T\theta(:,i);
  deltarms(i-1,1)=sqrt(delta(:,i-1)'*delta(:,i-1)/(N+1));
  clear RHS LHS dtheta;
endfor

aend=a(:,end);
Chxx=linspace(-1,1,NN+1)';
tt=b/2*(ChxChx+1);
TT=cos(acos(ChxChx)*n);
thetalin=TT*a;
thetaodelin=lsode("simpen",[theta0 dtheta0],tt);
diff=abs(thetalin(:,end)-thetaodelin(:,1));
# the following equals 6.7666681e-4 for N=2e3; NN=3e5; b=3; Iter=100
diffrms=sqrt(diff'*diff/(NN+1));
figure(1)
plot(t,thetainit);
figure(2)
plot(tt,thetalin(:,end))
figure(3)
semilogy(tt,diff)
figure(4)
semilogy(deltarms)