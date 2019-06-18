function dX=chen(X,t)
  a=40; b=3; c=28;
  dX=zeros(3,1);
  dX(1)=a*(X(2)-X(1));
  dX(2)=(c-a)*X(1)-X(1)*X(3)+c*X(2);
  dX(3)=X(1)*X(2)-b*X(3);
endfunction
