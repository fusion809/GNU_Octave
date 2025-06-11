% Based on ChatGPT code
function dtheta = double_pendulum_rhs(theta, t, params)
  % theta = [theta1; theta2; dtheta1; dtheta2]
  % params: struct with fields:
  %   m1, m2 = rod masses
  %   M1, M2 = bob masses
  %   L1, L2 = rod lengths
  %   g      = gravity
  %   gamma1, gamma2 = linear damping constants
  %   c1, c2         = quadratic damping constants

  % Unpack parameters
  m1 = params.m1;
  m2 = params.m2;
  M1 = params.M1;
  M2 = params.M2;
  L1 = params.L1;
  L2 = params.L2;
  g  = params.g;
  gamma1 = params.gamma1;
  gamma2 = params.gamma2;
  c1 = params.c1;
  c2 = params.c2;

  % Unpack state
  theta1 = theta(1);
  theta2 = theta(2);
  dtheta1 = theta(3);
  dtheta2 = theta(4);

  % Moments of inertia
  I1 = (1/3) * m1 * L1^2;
  I2 = (1/3) * m2 * L2^2;

  % Composite masses
  m1T = m1 + M1 + m2 + M2;
  m2T = m2 + M2;

  % Position-dependent terms
  delta = theta1 - theta2;
  sin_delta = sin(delta);
  cos_delta = cos(delta);

  % Mass terms
  a1 = I1 + M1*L1^2 + m2*L1^2 + M2*L1^2;
  a2 = I2 + M2*L2^2;
  b = m2 * L1 * L2/2 + M2 * L1 * L2;

  % RHS terms
  num1 = b * dtheta2^2 * sin_delta ...
         - (m2*L2/2 + M2*L2) * g * sin(theta2) * cos_delta ...
         - m1T * g * L1 * sin(theta1) ...
         - gamma1 * dtheta1 - c1 * abs(dtheta1) * dtheta1 ...
         + gamma2 * dtheta2 * cos_delta + c2 * abs(dtheta2) * dtheta2 * cos_delta;

  num2 = -b * dtheta1^2 * sin_delta ...
         + m1T * g * L1 * sin(theta1) * cos_delta ...
         - (m2*L2/2 + M2*L2) * g * sin(theta2) ...
         - gamma2 * dtheta2 - c2 * abs(dtheta2) * dtheta2 ...
         + gamma1 * dtheta1 * cos_delta + c1 * abs(dtheta1) * dtheta1 * cos_delta;

  den = a1 * a2 - b^2 * cos_delta^2;

  ddtheta1 = (a2 * num1 + b * cos_delta * num2) / den;
  ddtheta2 = (a1 * num2 + b * cos_delta * num1) / den;

  % Return state derivative
  dtheta = [dtheta1; dtheta2; ddtheta1; ddtheta2];
end