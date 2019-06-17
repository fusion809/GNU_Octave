# Clear all variables
clear all
# Show numbers in a long, yet intuitive form
format long g
# Number of collocation points
N         = 1000;
NN        = 100000;
# Integration interval
xx0       = 0;
xxf       = 50;
# Vector
n         = 0:N;
# Extrema grid
x         = -cos(pi*n'/N);
# Transformation to [xx0, xxf]
xtrans    = (xxf-xx0)/2*(x+1);
# linearly spaced grid
xlin      = linspace(-1,1,NN+1)';
Tlin      = cos(acos(xlin)*n);
xtranslin = (xlin+1)*(xxf-xx0)/2;
# Chebyshev polys of the first kind
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
F         = (cos(xtrans)).^2.*exp(-xtrans.^2);
# Initial condition; with definite integrals y(xx0) = 0
H(1,:)    = T(1,:);
F(1)      = 0;
# coefficient vector
a         = H\F;
# solution to ODE vector
y         = T*a;
# definite integral
y         = y-y(1);
# y lin
ylin      = Tlin*a;
# exact indefinite solution per Wolfram Alpha
yexact    = sqrt(pi)*(2*e*erf(xtranslin)+i*(erfi(1-i*xtranslin)-erfi(1+i*xtranslin)))/(8*exp(1));
yexact    = real(yexact);
# exact definite solution
yexact    = yexact-yexact(1);
# error to Chebyshev approximation
err       = abs(yexact-ylin);
errrms    = sqrt(err'*err/(N+1))
# y approximated by quadratic integration function
yquad     = quad("g", 0, inf);
# error in approximation
errquad   = abs(yexact(end)-yquad)