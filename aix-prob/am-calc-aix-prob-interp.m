clear all
format long e;
N=40;
NN=1000;
n=0:N;
b=10;
t=linspace(-1,1,NN+1)';
x=trans(b,t);
T=cos(acos(t)*n);
tol=[1e-12; 1e-13];
a=[quad("f0",-1,1,tol);
quad("f1",-1,1,tol);
quad("f2",-1,1,tol);
quad("f3",-1,1,tol);
quad("f4",-1,1,tol);
quad("f5",-1,1,tol);
quad("f6",-1,1,tol);
quad("f7",-1,1,tol);
quad("f8",-1,1,tol);
quad("f9",-1,1,tol);
quad("f10",-1,1,tol);
quad("f11",-1,1,tol);
quad("f12",-1,1,tol);
quad("f13",-1,1,tol);
quad("f14",-1,1,tol);
quad("f15",-1,1,tol);
quad("f16",-1,1,tol);
quad("f17",-1,1,tol);
quad("f18",-1,1,tol);
quad("f19",-1,1,tol);
quad("f20",-1,1,tol);
quad("f21",-1,1,tol);
quad("f22",-1,1,tol);
quad("f23",-1,1,tol);
quad("f24",-1,1,tol);
quad("f25",-1,1,tol);
quad("f26",-1,1,tol);
quad("f27",-1,1,tol);
quad("f28",-1,1,tol);
quad("f29",-1,1,tol);
quad("f30",-1,1,tol);
quad("f31",-1,1,tol);
quad("f32",-1,1,tol);
quad("f33",-1,1,tol);
quad("f34",-1,1,tol);
quad("f35",-1,1,tol);
quad("f36",-1,1,tol);
quad("f37",-1,1,tol);
quad("f38",-1,1,tol);
quad("f39",-1,1,tol);
quad("f40",-1,1,tol);
quad("f41",-1,1,tol);
quad("f42",-1,1,tol);
quad("f43",-1,1,tol);
quad("f44",-1,1,tol);
quad("f45",-1,1,tol);
quad("f46",-1,1,tol);
quad("f47",-1,1,tol);
quad("f48",-1,1,tol);
quad("f49",-1,1,tol);
quad("f50",-1,1,tol);
quad("f51",-1,1,tol);
quad("f52",-1,1,tol);
quad("f53",-1,1,tol);
quad("f54",-1,1,tol);
quad("f55",-1,1,tol);
quad("f56",-1,1,tol);
quad("f57",-1,1,tol);
quad("f58",-1,1,tol);
quad("f59",-1,1,tol);
quad("f60",-1,1,tol)];
y=T(:,1:N)*a(1:N);

# Accuracy to 9 decimal places exists for:
# Ai(x)
# 1/sqrt(1+x^2)
# tanh(x)/sqrt(1+x^2)
figure(1)
plot(x,y)
figure(2)
plot(x,f(t))

rms=sqrt((y-f(t))'*(y-f(t)))