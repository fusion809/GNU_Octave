function y=fn(n,x)
  if (n==0)
    y=(((f(x)).*cos(acos(x)*n))./sqrt(1-x.^2))*(1/pi);
  else
    y=(((f(x)).*cos(acos(x)*n))./sqrt(1-x.^2))*(2/pi);
  endif
endfunction
