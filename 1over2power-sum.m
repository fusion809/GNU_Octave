clear all
format long g;
N=53;
i=1:N;
t=(1/2).^(i);
tot=sum(t);
err=abs(1-tot);
N
log10(err)
