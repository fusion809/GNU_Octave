function y=f(x)
  b=10;
  y=exp(trans(b,x)).*log(trans(b,x)+1).*airy(0,trans(b,x));
endfunction
