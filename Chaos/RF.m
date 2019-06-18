function dX=RF(X,t)
  # Rabinovich–Fabrikant equations
  alpha=1.1; gamma=0.87;
  dX=zeros(3,1);
  dX(1)=X(2)*(X(3)-1+X(1)^2)+gamma*X(1);
  dX(2)=X(1)*(3*X(3)+1-X(1)^2)+gamma*X(2);
  dX(3)=-2*X(3)*(alpha+X(1)*X(2));
endfunction
