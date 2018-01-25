# Drag pendulum calc
N   = 1000;
t0  = 0;
t1  = 1.0356;
x0  = 0;
dx0 = 10;
y0  = 0;
dy0 = 10;

function Der = f0(X, t)
  alpha  = 0.1;
 # alpha  = 0; # Zero drag
  g      = 9.8;
  Der    = zeros(4,1);
  
  Der(1) = X(2);
  Der(2) = - alpha * ( X(2)^2 + X(4)^2 );
  Der(3) = X(4);
  Der(4) = - alpha * ( X(2)^2 + X(4)^2 ) - g * X(3);
endfunction

function Der = f1(X, t)
  alpha  = 0.05;
 # alpha  = 0; # Zero drag
  g      = 9.8;
  Der    = zeros(4,1);
  
  Der(1) = X(2);
  Der(2) = - alpha * ( X(2)^2 + X(4)^2 );
  Der(3) = X(4);
  Der(4) = - alpha * ( X(2)^2 + X(4)^2 ) - g * X(3);
endfunction

function Der = f2(X, t)
  alpha  = 0.01;
 # alpha  = 0; # Zero drag
  g      = 9.8;
  Der    = zeros(4,1);
  
  Der(1) = X(2);
  Der(2) = - alpha * ( X(2)^2 + X(4)^2 );
  Der(3) = X(4);
  Der(4) = - alpha * ( X(2)^2 + X(4)^2 ) - g * X(3);
endfunction

function Der = f3(X, t)
  alpha  = 0.005;
 # alpha  = 0; # Zero drag
  g      = 9.8;
  Der    = zeros(4,1);
  
  Der(1) = X(2);
  Der(2) = - alpha * ( X(2)^2 + X(4)^2 );
  Der(3) = X(4);
  Der(4) = - alpha * ( X(2)^2 + X(4)^2 ) - g * X(3);
endfunction

function Der = f4(X, t)
  alpha  = 0.001;
 # alpha  = 0; # Zero drag
  g      = 9.8;
  Der    = zeros(4,1);
  
  Der(1) = X(2);
  Der(2) = - alpha * ( X(2)^2 + X(4)^2 );
  Der(3) = X(4);
  Der(4) = - alpha * ( X(2)^2 + X(4)^2 ) - g * X(3);
endfunction

function Der = f5(X, t)
  alpha  = 0;
 # alpha  = 0; # Zero drag
  g      = 9.8;
  Der    = zeros(4,1);
  
  Der(1) = X(2);
  Der(2) = - alpha * ( X(2)^2 + X(4)^2 );
  Der(3) = X(4);
  Der(4) = - alpha * ( X(2)^2 + X(4)^2 ) - g * X(3);
endfunction

t  = linspace(t0, t1, N+1)';
X0 = [x0; dx0; y0; dy0];
X  = lsode("f0", X0, t);
X1 = lsode("f1", X0, t);
X2 = lsode("f2", X0, t);
X3 = lsode("f3", X0, t);
X4 = lsode("f4", X0, t);
X5 = lsode("f5", X0, t);

# x y plot
figure(1)
plot(X5(:,1), X5(:,3), 'k', X4(:,1), X4(:,3), 'r', X3(:,1), X3(:,3), 'g', X2(:,1), X2(:,3), 'm', X1(:,1), X1(:,3), 'b', X(:,1), X(:,3), 'c')
legend('alpha=0', 'alpha=0.001', 'alpha=0.005', 'alpha=0.01', 'alpha=0.05', 'alpha=0.1')