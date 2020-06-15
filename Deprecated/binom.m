function y=binom(n,p,x)
  y=(p^x*(1-p)^(n-x))*nchoosek(n,x);
endfunction
