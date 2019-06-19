function dX=rossler(X,t)
  a=0.36; b=0.1; c=30;
  dX=zeros(3,1);
  dX(1)=-X(2)-X(3);
  dX(2)=X(1)+a*X(2);
  dX(3)=b+X(3)*(X(1)-c);
endfunction
