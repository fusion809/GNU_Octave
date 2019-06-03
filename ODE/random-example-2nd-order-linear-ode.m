# Equation to be solved:
# y'' + y' + x^2 y = 0
# y(0)  = 0
# y'(0) = 1

clear all
N = 2000;
b = 20;
#x = linspace(a,b,N+1)';
#y = lsode("randex",[0; 1]',x);

#plot(x,y(:,1),'r',x,y(:,2),'g')

n         = 0:N;
x         = -cos(pi*n'/N);
T         = cos(acos(x)*n);
# Now for arrays that do not include the endpoints
xsub      = x(2:N);
Tsub      = T(2:N,:);
Usub      = diag(1./sqrt(1-xsub.^2))*sin(acos(xsub)*n);
dTsub     = Usub*diag(n);
dT        = [-(n.^2).*(-1).^(n); dTsub ; (n).^2];
d2Tsub    = diag(1./(1-xsub.^2))*(diag(xsub)*Usub-Tsub*diag(n))*diag(n);

Hsub      = 4/(b^2)*d2Tsub+2/b*dTsub+(b^2)/4*diag((xsub+ones(N-1,1))).^2*Tsub;
H         = [Tsub(1,:); Hsub; 2/b*dTsub(1,:)];
RHS       = zeros(N+1,1);
#RHS(1)    = 0; # y(0)
RHS(N+1)  = 1; # y'(0)

xx        = b/2*(x+1); # Real x on [0, b]
a         = H\RHS;
y         = T*a;
dy        = 2/b*dT*a;

plot(xx,y,'r',xx,dy,'g')