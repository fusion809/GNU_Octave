function y=f(x)
  b=10;
  y=exp(trans(b,x)).*airy(0,trans(b,x));
endfunction
