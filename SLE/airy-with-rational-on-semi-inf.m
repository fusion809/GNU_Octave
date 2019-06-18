clear all
format long e;
N=5000;
perc=0.266
Nfrag=round(perc*N);
n=0:N;
# 75 at N=1k; 162 for N=3k; 225 for N=5k;
L=226;
t=pi+n'*pi/N;
x=cos(t);
y=L*(cot(t/2)).^2;
ysub=y(2:N);
ylintrans = L/2*(x+1);
# Chebyshev polys of the first kind
T         = cos(t*n);
# Now for arrays that do not include the endpoints
xsub      = x(2:N);
Tsub      = T(2:N,:);
Usub      = diag(1./sqrt(1-xsub.^2))*sin(acos(xsub)*n);
dTsub     = Usub*diag(n);
dT        = [-(n.^2).*(-1).^(n); dTsub ; (n).^2];
dTL       = diag(2*L./((y+L).^2))*dT;
dTlin     = 2/L*dT;
D1        = dTL/T;
D1lin     = dTlin/T;
D2        = D1^2;
D2lin     = D1lin^2;
H         = -D2+diag(y);
Hlin      = -D2lin+diag(ylintrans);
H         = H(2:N,2:N);
Hlin      = Hlin(2:N,2:N);
[Y, Lam]  = eig(H);
[Ylin, Lamlin] = eig(Hlin);
Lam       = diag(Lam);
Lamlin    = diag(Lamlin);
[Lam, IX] = sort(Lam, 'ascend');
[Lamlin, IXlin] = sort(Lamlin, 'ascend');
Y         = Y(:,IX);
Ylin      = Ylin(:,IXlin);
aizerochk = abs(airy(0,-Lam));
aizerochklin = abs(airy(0,-Lamlin));
L
errrms    = sqrt(aizerochk(1:Nfrag)'*aizerochk(1:Nfrag)/(Nfrag+1))
errrmslin    = sqrt(aizerochklin(1:Nfrag)'*aizerochklin(1:Nfrag)/(Nfrag+1))

#figure(1)
#plot(ysub(1:3*N/5), Y(1:3*N/5,1));
#figure(2)
#plot(ysub(1:3*N/5), Y(1:3*N/5,2));
#figure(3)
#plot(ysub(1:3*N/5), Y(1:3*N/5,Nfrag));