function dydt = double_elastic_pendulum(t, y, p)
  % Parameters
  m1 = p(1); m2 = p(2); g = p(3);
  k1 = p(4); k2 = p(5); l1 = p(6); l2 = p(7);

  % Unpack state
  r1 = y(1); dr1 = y(2);
  th1 = y(3); dth1 = y(4);
  r2 = y(5); dr2 = y(6);
  th2 = y(7); dth2 = y(8);

  % Angle difference
  delta = th2 - th1;

  % Initialise matrix A and rhs b
  A = eye(4);
  b = zeros(4,1);

  % External dissipative/generalized forces (assume zero here)
  Qr1 = 0; Qr2 = 0; Qt1 = 0; Qt2 = 0;

  % Equation for ddot{r1}
  b(1) = r1*dth1^2 - g*sin(th1) + ...
         (m2/(m1 + m2)) * (cos(delta)*r2*dth2^2 + sin(delta)*2*dr2*dth2) - ...
         k1*(r1 - l1)/(m1 + m2) + Qr1/(m1 + m2);
  A(1,2) = - (m2/(m1 + m2)) * cos(delta);         % for ddot{r2}
  A(1,4) = - (m2/(m1 + m2)) * r2 * sin(delta);    % for ddot{theta2}

  % Equation for ddot{r2}
  b(2) = - cos(delta)*r1*dth1^2 - sin(delta)*dr1*dth1 + ...
         r2*dth2^2 - g*sin(th2) + (Qr2 - k2*(r2 - l2))/m2;
  A(2,1) = -cos(delta);         % for ddot{r1}
  A(2,3) = -r1*sin(delta);      % for ddot{theta1}

  % Equation for ddot{theta1}
  b(3) = -2*dr1*dth1/r1 - g*cos(th1)/r1 - ...
         (m2/((m1 + m2)*r1)) * ...
         (cos(delta)*2*dr2*dth2 - sin(delta)*r2*dth2^2) + ...
         Qt1 / ((m1 + m2)*r1^2);
  A(3,4) = - (m2/((m1 + m2)*r1)) * cos(delta);   % for ddot{theta2}
  A(3,2) = - (m2/((m1 + m2)*r1)) * sin(delta);   % for ddot{r2}

  % Equation for ddot{theta2}
  b(4) = -2*dr2*dth2/r2 - (cos(delta)/r2)*(2*dr1*dth1) - ...
         (sin(delta)/r2)*(r1*dth1^2) - g*cos(th2)/r2 + ...
         Qt2 / (m2 * r2^2);
  A(4,3) = - (cos(delta)/r2) * r1;    % for ddot{theta1}
  A(4,1) =   sin(delta)/r2;          % for ddot{r1}

  % Solve for [ddr1; ddtheta1; ddr2; ddtheta2]
  x = A \ b;

  % Build dydt
  dydt = zeros(8,1);
  dydt(1) = dr1;
  dydt(2) = x(1);  % ddot{r1}
  dydt(3) = dth1;
  dydt(4) = x(2);  % ddot{theta1}
  dydt(5) = dr2;
  dydt(6) = x(3);  % ddot{r2}
  dydt(7) = dth2;
  dydt(8) = x(4);  % ddot{theta2}
end