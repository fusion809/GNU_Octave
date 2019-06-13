function dX=lorenz(X,t)
  dX=zeros(3,1);
  beta=8/3;
  rho=28;
  sigma=10;
  dX(1)=sigma*(X(2)-X(1));
  dX(2)=X(1)*(rho-X(3))-X(2);
  dX(3)=X(1)*X(2)-beta*X(3);
endfunction
