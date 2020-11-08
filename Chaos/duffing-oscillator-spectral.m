clear all
N         = 3e3;
NN        = 1e5;
Iter      = 15;
b         = 5;
delta     = 0.02;
gamma     = 8;
omega     = 0.5;
alpha     = 1;
beta      = 5;
x0        = 0;
dx0       = 0;
n         = 0:N;
parat     = pi+n'*pi/N;
chx       = cos(parat);
t         = b/2*(chx+1);
T         = cos(parat*n);
# Now for arrays that do not include the endpoints
chxsub    = chx(2:N);
paratsub  = parat(2:N);
Tsub      = T(2:N,:);
Usub      = diag(1./sqrt(1-chxsub.^2))*sin(paratsub*n);
dTsub     = Usub*diag(n);
dT        = [-(n.^2).*(-1).^(n); dTsub ; (n).^2];
D1        = 2/b*dT/T;
D2        = D1^2;
xinit     = lsode("duffing",[x0 dx0],t);
xinit     = xinit(:,1);
x         = zeros(N+1,Iter); 
x(:,1)    = xinit;
a(:,1)    = T\xinit;

for i=2:Iter
  dx=D1*x(:,i-1);
  RHS=-D2*x(:,i-1)-delta*dx-beta*x(:,i-1).^3+gamma*cos(omega*t);
  LHS=D2+delta*D1+3*beta*diag(x(:,i-1).^2);
# The different recipe for the last iter. relates to the fact the least squares solve leaves the init conds.
# not strictly adhered to
# Regardless of whether this if conditional is used, or whether least squares is used the result is still
# wrong (as applying the equation to the result gives erroneous values at the endpoints, or the init conds are
# not fully met.
  if ( i==Iter )
    LHS(1,1)=1;
    LHS(1,2:N+1)=0;
    LHS(N+1,:)=D1(1,:);
    RHS(1)=x0-x(1,i-1);
    RHS(N+1)=dx0-dx(1);
  else
    LHS(N+2,1)=1;
    LHS(N+2,2:N+1)=0;
    LHS(N+3,:)=D1(1,:);
    #  rcon(i-1)=rcond(LHS);
    # Initial conditions
    # designed to reverse any creeking in init. conds.
    RHS(N+2)=x0-x(1,i-1);
    RHS(N+3)=dx0-dx(1);
  endif
  
  Delta(:,i-1)=LHS\RHS;
  x(:,i)=x(:,i-1)+Delta(:,i-1);
  a(:,i)=T\x(:,i);
  Deltarms(i-1,1)=sqrt(Delta(:,i-1)'*Delta(:,i-1)/(N+1));
  clear RHS LHS dx;
endfor

aend=a(:,end);
chxchx=linspace(-1,1,NN+1)';
tt=b/2*(chxchx+1);
TT=cos(acos(chxchx)*n);
xlin=TT*a;
diff=abs(xlin(:,end)-xlin(:,1));
# the following equals 6.7666681e-4 for N=2e3; NN=3e5; b=3; Iter=100
diffrms=sqrt(diff'*diff/(NN+1));
figure(1)
plot(tt,xlin(:,1));
title("xlin(:,1) vs tt")
figure(2)
plot(tt,xlin(:,end))
title("xlin(:,end) vs tt")
figure(3)
semilogy(tt,diff)
title("xlin(:,end)-xlin(:,1) vs tt")
figure(4)
semilogy(Deltarms)
title("Logarithmic graph of root mean square of Delta across iterations")
figure(5)
plot(tt,xlin(:,1),'r',tt,xlin(:,end),'g')
title("xlin(:,1) (red) and xlin(:,end) (green) vs tt")
