% Based on ChatGPT code
function dth = DP2(params, TH, t)
    g = params(1);
    r1 = params(2);
    r2 = params(3);
    mr1 = params(4);
    mr2 = params(5);
    mb1 = params(6);
    mb2 = params(7);
    gamma1 = params(8);
    gamma2 = params(9);
    c1 = params(10);
    c2 = params(11);
    m = m1 + m2;
    M = m1 + m2;
    dth = zeros(4,1);
    theta1 = TH(1);
    theta2 = TH(2);
    dtheta1 = TH(3);
    dtheta2 = TH(4); 
    dth(1) = dtheta1;
    dth(2) = dtheta2; % Moments of inertia
  % Total masses at each joint
  m1 = mb1 + mr1/2;    % mass at first hinge (bob + 1/2 rod mass)
  m2 = mb2 + mr2/2;    % mass at second hinge (bob + 1/2 rod mass)

  % Moment of inertia for uniform rods: (1/3)*m*L^2 if pivoted at one end
  I1 = (1/3)*mr1*r1^2;
  I2 = (1/3)*mr2*r2^2;

  % Distance from pivot to center of mass
  l1 = (mb1*r1 + (mr1/2)*(r1/2)) / (mb1 + mr1);  % composite CoM of first segment
  l2 = (mb2*r2 + (mr2/2)*(r2/2)) / (mb2 + mr2);  % composite CoM of second segment

  delta = theta2 - theta1;
  sin_delta = sin(delta);
  cos_delta = cos(delta);

  % Effective masses for equations of motion
  M11 = m1*r1^2 + m2*(r1^2 + r2^2 + 2*r1*r2*cos_delta) + I1 + I2;
  M12 = m2*(r2^2 + r1*r2*cos_delta) + I2;
  M21 = M12;
  M22 = m2*r2^2 + I2;

  % Right-hand side terms (including gravity and damping)
  RHS1 = -m2*r1*r2*sin_delta*dtheta2^2 - ...
         (m1*l1 + m2*r1)*g*cos(theta1) - m2*r2*g*cos(theta2) - ...
         gamma1*dtheta1 - c1*abs(dtheta1)*dtheta1;

  RHS2 = m2*r1*r2*sin_delta*dtheta1^2 - ...
         m2*r2*g*cos(theta2) - ...
         gamma2*dtheta2 - c2*abs(dtheta2)*dtheta2;

  % Solve linear system M * ddtheta = RHS
  M = [M11, M12; M21, M22];
  RHS = [RHS1; RHS2];

  dth(3:4) = M \ RHS;
endfunction