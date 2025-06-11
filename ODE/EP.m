function dth = EP(params, th, t)
    g      = params.g;
    l0     = params.l0;
    % In my working out, these are b', c' and k', 
    % which are b/m, c/m and k/m, respectively. 
    b      = params.b;
    c      = params.c;
    k      = params.k;
    theta  = th(1);
    z      = th(2);
    dtheta = th(3);
    dz     = th(4);
    dth(1) = th(3);
    dth(2) = th(4);
    dth(3) = -1/(z+l0)*(2*dtheta*dz + g*cos(theta)) - dtheta*(b+c*sqrt(dz^2+dtheta^2*(z+l0)^2));
    dth(4) = dtheta^2*(z+l0) - k*z - g*sin(theta) - dz*(b+c*sqrt(dz^2+dtheta^2*(z+l0)^2));
endfunction