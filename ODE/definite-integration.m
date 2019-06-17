clear all
format long g
N         = 1000;
xx0       = 0;
xxf       = 50;
n         = 0:N;
x         = -cos(pi*n'/N);
# Transformation to [xx0, xxf]
xx        = (xxf-xx0)/2*(x+1);
T         = cos(acos(x)*n);
# Now for arrays that do not include the endpoints
xsub      = x(2:N);
Tsub      = T(2:N,:);
Usub      = diag(1./sqrt(1-xsub.^2))*sin(acos(xsub)*n);
dTsub     = Usub*diag(n);
# Add the endpoints
dT        = [-(n.^2).*(-1).^(n); dTsub ; (n).^2];
# LHS without a (coefficient vector)
H         = 2/(xxf-xx0)*dT;
# RHS, essentially what we're trying to integrate
F         = cos(xx).*exp(-xx.^2);
# Initial condition; with definite integrals y(xx0) = 0
H(1,:)    = T(1,:);
F(1)      = 0;
# coefficient vector
a         = H\F;
# solution vector
y         = T*a;
y         = y-y(1);
yexact    = sqrt(pi)*(erf(xx-1/2*i)+erf(xx+1/2*i))/(4*exp(1/4))
yexact    = yexact-yexact(1);
err       = abs(yexact-y);
errrms    = sqrt(err'*err)
yquad     = quad("g", 0, inf);
errquad   = abs(yexact(end)-yquad)