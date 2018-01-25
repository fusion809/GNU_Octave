clear all

x=-2.34;
while abs(airy(0,x))>1e-10
  x=x-airy(0,x)/airy(1,x);
endwhile

format long g;
airy(0,x)