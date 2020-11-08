clear all
N=100;
n=0:N;
x=linspace(0,pi,N+1)';
y=cos(x);
dy=-sin(x);
for i=1:length(x)-4
  dyapp(i,1)=(y(i,1)-8*y(i+1,1)+8*y(i+3,1)-y(i+4,1))/(12*(x(2)-x(1)));
  diff(i,1)=dyapp(i,1)-dy(i+2,1);
endfor

for i=1:length(x)-6
  dyapp6(i,1)=(-y(i,1)+9*y(i+1,1)-45*y(i+2,1)+45*y(i+4,1)-9*y(i+5,1)+y(i+6,1))/(60*(x(2)-x(1)));
  diff6(i,1)=dyapp6(i,1)-dy(i+3,1);
endfor

rms=sqrt(diff'*diff/length(diff));
rms6=sqrt(diff6'*diff6/length(diff6));