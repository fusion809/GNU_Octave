function dth = sinpen(params, th, t)
    dth = zeros(2,1);
    g = params(1);
    r = params(2);
    mr = params(3);
    mb = params(4); 
    gamma = params(5);
    c = params(6);
    dth(1) = th(2);
    dth(2) = -((1/2*mr+mb)*g*cos(th(1))+gamma*th(2)+c*abs(th(2))*th(2))/((1/3*mr + mb)*r);
endfunction