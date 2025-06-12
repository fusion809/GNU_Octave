function dth = DP(params, TH, t)
    g = params(1);
    r1 = params(2);
    r2 = params(3);
    m1 = params(4);
    m2 = params(5);
    M1 = params(6);
    M2 = params(7);
    gamma1 = params(8);
    gamma2 = params(9);
    c1 = params(10);
    c2 = params(11);
    m = m1 + m2;
    dth = zeros(4,1);
    theta1 = TH(1);
    theta2 = TH(2);
    dtheta1 = TH(3);
    dtheta2 = TH(4); 
    dth(1) = dtheta1;
    dth(2) = dtheta2;
    % Based on ChatGPT equations
    D = -m2*r1*r2*dtheta2^2*sin(theta1-theta2)-m*g*r1*cos(theta1);
    E = m2*r1*r2*dtheta1^2*sin(theta1-theta2) - m2*g*r2*cos(theta2);
    Delta = m*m2*r1^2*r2^2-(m2*r1*r2*cos(theta1-theta2))^2;
    dth(3) = (m2*r2^2*D-m2*r1*r2*cos(theta1-theta2)*E)/Delta;
    dth(4) = (-m2*r1*r2*cos(theta1-theta2)*D+m*r1^2*E)/Delta;
endfunction