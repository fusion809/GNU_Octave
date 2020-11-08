function dx=duffing(x,t)
  delta=0.02;
  beta=5;
  alpha=1;
  omega=0.5;
  gamma=8;
  dx=zeros(2,1);
  dx(1)=x(2);
  dx(2)=-delta*x(2)-beta*x(1)^3+gamma*cos(omega*t);
endfunction
