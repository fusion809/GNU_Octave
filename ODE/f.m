function dy=f(y,x)
  dy=zeros(2,1);
  delta=0.02; alpha=1; beta=5; gamma=8; omega=0.5;
  dy(1)=y(2);
  dy(2)=-delta*y(2)-alpha*y(1)-beta*y(1)^3-gamma*cos(omega*x)
endfunction
