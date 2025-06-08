clear all
N         = 2.5e3;
Iter      = 5;
g         = 9.81;
l         = 1;
k         = 1:100;
period    = 4*sqrt(l/g)*(pi/2 + 2*pi*sum(exp(k*pi)./(1+exp(2*k*pi))));
theta0    = 0;
dtheta0   = 0;
n         = 0:N;
parat     = pi+n'*pi/N;
x         = cos(parat);
t         = period/2*(x+1);
T         = cos(parat*n);
% Now for arrays that do not include the endpoints
xsub      = x(2:N);
paratsub  = parat(2:N);
Tsub      = T(2:N,:);
Usub      = diag(1./sqrt(1-xsub.^2))*sin(paratsub*n);
dTsub     = Usub*diag(n);
dT        = [-(n.^2).*(-1).^(n); dTsub ; (n).^2];
D1        = 2/period*dT/T;
D2        = D1^2;
thetainit = -pi/2*(1-cos(2*pi/period*t));
theta     = zeros(N+1,Iter); 
theta(:,1)= thetainit;

for i=2:Iter
  dtheta=D1*theta(:,i-1);
  RHS=-D2*theta(:,i-1)-g/l*cos(theta(:,i-1));
  LHS=D2-g/l*diag(sin(theta(:,i-1)));
  LHS(1,1)=1;
  LHS(1,2:N+1)=0;
  LHS(N+1,:)=D1(1,:);
  % Initial conditions
  % designed to reverse any creeking in init. conds.
  RHS(1)=theta0-theta(1,i-1);
  RHS(N+1)=dtheta0-dtheta(1);
  delta(:,i-1)=LHS\RHS;
  theta(:,i)=theta(:,i-1)+delta(:,i-1);
  deltarms(i-1,1)=sqrt(delta(:,i-1)'*delta(:,i-1)/(N+1));
  clear RHS LHS dtheta;
endfor

%aend=a(:,end);
figure(1)
plot(t,thetainit);
title("Simple pendulum: \\theta^{(0)} vs t", "fontsize", 25);
xlabel("t", "fontsize", 20);
ylabel("\\theta^{(0)}", "fontsize", 20, "rotation", 0);
set(gca, 'fontsize', 16)
xlim([0 period]);
ylim([-pi 0]);
figure(2)
plot(t,theta(:,end))
title("Simple pendulum: \\theta^{(f)} vs t", "fontsize", 25)
xlabel("t", "fontsize", 20);
ylabel("\\theta^{(f)}", "fontsize", 20, "rotation", 0);
xlim([0 period]);
ylim([-pi 0]);
set(gca, 'fontsize', 16)
figure(3)
plot(t, D1*theta(:,end))
xlim([0 period])
ylim([-sqrt(2*g/l) sqrt(2*g/l)])
title("Simple pendulum: d\\theta^{(f)}/dt vs t", "fontsize", 25)
xlabel("t", "fontsize", 20);
h=ylabel("d\\theta^{(f)}/dt", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 1.5/10;      % Shift left by 1.5 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
figure(4)
plot(t, [theta(:,1) theta(:,end)])
xlim([0 period])
ylim([-pi 0])
title("Simple pendulum: \\theta^{(0/f)} vs t", "fontsize", 25)
xlabel("t", "fontsize", 20);
ylabel("\\theta^{(0/f)}", "fontsize", 20, "rotation", 0);
set(gca, 'fontsize', 16)
legend('\theta^{(0)}', '\theta^{(f)}')
figure(5)
residual = D2*[theta(:,1) theta(:,end)]+ g/l*cos([theta(:,1) theta(:,end)]);
semilogy(t, abs(residual))
xlim([0 period])
ylim([min(min(abs(log10(abs(residual))))) max(max(log10(abs(residual))))])
title("Simple pendulum: residual vs t", "fontsize", 25)
xlabel("t", "fontsize", 20);
h=ylabel("Residual", "fontsize", 20, "rotation", 0);
pos = get(h, 'position');            % Get current [x, y, z] position
new_pos = pos;                       
new_pos(1) = new_pos(1) - 3/10;      % Shift left by 1.5 cm (0.01 m = 1/10 "normalized units")
set(h, 'position', new_pos);         % Apply new position
set(gca, 'fontsize', 16)
legend('\theta^{(0)}', '\theta^{(f)}')