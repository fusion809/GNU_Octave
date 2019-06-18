function y=gsimp(x)
  y=airy(0,x).*exp(-x).*log(x+1);
endfunction
