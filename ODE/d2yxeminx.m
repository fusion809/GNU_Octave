function dy=d2yxeminx(y,x)
	dy=zeros(2,1);
	dy(1)=y(2);
	dy(2)=(x.^2)./sqrt(1+x.^2)*y(1);
endfunction
