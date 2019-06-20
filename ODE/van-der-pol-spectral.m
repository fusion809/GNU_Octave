clear all
N         = 1e3;
NN        = 1e4;
Iter      = 10;
b         = 10;
mu        = 4;
x0        = 0;
dx0       = 1;
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
Chxx      = linspace(-1,1,NN+1)';
tt        = b/2*(Chxx+1);
TT        = cos(acos(Chxx)*n);
xodelin   = lsode("vanderpol",[x0 dx0],t);
xinit     = xodelin(:,1);
x         = zeros(N+1,Iter); 
x(:,1)    = xinit;
a(:,1)    = T\xinit;

for i=2:Iter
  dx=D1*x(:,i-1);
  RHS=-D2*x(:,i-1)+mu*(1-x(:,i-1).^2).*dx-x(:,i-1);
  LHS=D2-2*mu*diag(x(:,i-1).*dx)-mu*diag(x(:,i-1).^2)*D1+diag(ones(N+1,1));
  LHS(1,1)=1;
  LHS(1,2:N+1)=0;
  LHS(N+1,:)=D1(1,:);
  rcon(i-1)=rcond(LHS);
  # Initial conditions
  # designed to reverse any creeking in init. conds.
  RHS(1)=x0-x(1,i-1);
  RHS(N+1)=dx0-dx(1);
  delta(:,i-1)=LHS\RHS;
  x(:,i)=x(:,i-1)+delta(:,i-1);
  a(:,i)=T\x(:,i);
  deltarms(i-1,1)=sqrt(delta(:,i-1)'*delta(:,i-1)/(N+1));
  clear RHS LHS dx;
endfor
xlin      = TT*a;
aend=a(:,end);
xlin=TT*a;
diff=abs(xlin(:,end)-xlin(:,1));
# the following equals 6.7666681e-4 for N=2e3; NN=3e5; b=3; Iter=100
diffrms=sqrt(diff'*diff/(NN+1));
figure(1)
plot(t,xinit);
figure(2)
plot(tt,xlin(:,end))
figure(3)
semilogy(tt,diff)
figure(4)
semilogy(deltarms)