clear all
N=50;
NN=1000;
n=0:N;
t=linspace(-1,1,NN+1)';
T=cos(acos(t)*n);
a=[quad("f0",-1,1);
quad("f1",-1,1);
quad("f2",-1,1);
quad("f3",-1,1);
quad("f4",-1,1);
quad("f5",-1,1);
quad("f6",-1,1);
quad("f7",-1,1);
quad("f8",-1,1);
quad("f9",-1,1);
quad("f10",-1,1);
quad("f11",-1,1);
quad("f12",-1,1);
quad("f13",-1,1);
quad("f14",-1,1);
quad("f15",-1,1);
quad("f16",-1,1);
quad("f17",-1,1);
quad("f18",-1,1);
quad("f19",-1,1);
quad("f20",-1,1);
quad("f21",-1,1);
quad("f22",-1,1);
quad("f23",-1,1);
quad("f24",-1,1);
quad("f25",-1,1);
quad("f26",-1,1);
quad("f27",-1,1);
quad("f28",-1,1);
quad("f29",-1,1);
quad("f30",-1,1);
quad("f31",-1,1);
quad("f32",-1,1);
quad("f33",-1,1);
quad("f34",-1,1);
quad("f35",-1,1);
quad("f36",-1,1);
quad("f37",-1,1);
quad("f38",-1,1);
quad("f39",-1,1);
quad("f40",-1,1);
quad("f41",-1,1);
quad("f42",-1,1);
quad("f43",-1,1);
quad("f44",-1,1);
quad("f45",-1,1);
quad("f46",-1,1);
quad("f47",-1,1);
quad("f48",-1,1);
quad("f49",-1,1);
quad("f50",-1,1)];
y=T(:,1:20)*a(1:20);

# Accuracy to 9 decimal places exists for:
# Ai(x)
# 1/sqrt(1+x^2)
# tanh(x)/sqrt(1+x^2)
figure(1)
plot(t,y)
figure(2)
plot(t,f(t))

rms=sqrt((y-f(t))'*(y-f(t)));