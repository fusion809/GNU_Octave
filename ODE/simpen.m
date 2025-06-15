function dy = simpen(params, y,t)
  dy    = zeros(1,2);
  g     = params(1);
  l     = params(2);
  dy(1) = y(2);
  dy(2) = -g/l*cos(y(1));
endfunction
