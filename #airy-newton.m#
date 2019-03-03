clear all

x=-2.34;
while abs(airy(0,x))>1e-10
  x=x-airy(0,x)/airy(1,x);
endwhile

error=airy(0,x);
format long g;
printf("x is %d\n", x)
printf("Ai(x) is %d\n", error)
