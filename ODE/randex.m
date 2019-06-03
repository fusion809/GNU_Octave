function dy = randex(y,x)
  dy = zeros(2,1);
  dy(1) = y(2);
  dy(2) = -y(2)-x^2*y(1);
