% This script gives an incorrect result for theta
N = 1000;
n = 0:N;
r = 1;
g = 9.81;
chebt = pi*n'/N;
x = cos(chebt);

D1base = cheb(N);
D2base = D1base^2;
k = 1:100;
period = 4*sqrt(r/g)*(pi/2 + 2*pi*sum(exp(k*pi)./(1+exp(2*k*pi))));
t = period/4 * (x+1);
D1 = 4/period * D1base;
E1 = D1(2:N, 2:N);
D2 = 16/period^2 * D2base;
E2 = D2(2:N, 2:N);
theta0 = 0;
thetadot0 = 0;
thetainit = -pi/2 * (1-cos(2*pi*t/period));
theta = thetainit;
% dthetainitan = -pi^2/period * sin(2*pi*t/period);
% d2thetainitan = -2*pi^3/period^2 * cos(2*pi*t/period);
% dthetainit = D1*thetainit;
% d2thetainit = D2*thetainit;
% dthetainiterr = dthetainit - dthetainitan;
% d2thetainiterr = d2thetainit - d2thetainitan;
i = 0;
delta = 5*ones(N-1,1);
thetamid = theta(2:N);
while (sum(abs(delta)) > 1e-10)
%while (i < 2)
    %RHS = sqrt(thetadot0^2 + 2*g/r*(sin(theta0)-sin(thetamid)));
    % = E1 + diag(g/r*cos(thetamid)./RHS);
    H = E2 - diag(g/r*sin(thetamid));
    %dtheta = E1*thetamid;
    %F = -dtheta + RHS;
    F = -E2*thetamid - g/r*cos(thetamid);
    %H(1,:) = [1 zeros(1,N)];
    %H(N+1,:) = D1(1,:);
    %H(N+1,:) = [zeros(1,N) 1];
    %F(1) = -theta(1);
    %F(N+1) = -dtheta(1);
    %F(N+1) = -theta(N+1);
    delta = H\F;
    thetamid = thetamid + delta;
    i = i + 1;
end
