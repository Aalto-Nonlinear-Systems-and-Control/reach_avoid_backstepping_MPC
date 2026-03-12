% example for rocket landing with bounded control inputs
clc; clear; close all;

% define the system dynamics using symbolic variables
% 7D state: [x1; x2; x3; x4; x5; x6; x7] = [x; y; theta; vx; vy; omega; m]
%   x1 = horizontal position x
%   x2 = vertical position y
%   x3 = tilt angle theta
%   x4 = horizontal velocity vx
%   x5 = vertical velocity vy
%   x6 = angular velocity omega
%   x7 = mass m
% 2D control input: [u1; u2] = [main thrust T; lateral thrust T_lat]

% physical parameters
g_val = 9.8; % gravitational acceleration (m/s^2)
L_val = 50.0; % rocket length (m)
Isp = 300.0; % specific impulse (s)

syms x1 x2 x3 x4 x5 x6 x7 u1 u2;
x_vars_sym = [x1; x2; x3; x4; x5; x6; x7];
u_vars_sym = [u1; u2];

% drift term f(X) = [vx; vy; omega; 0; -g; 0; 0]
fx_sym = [x4; x5; x6; 0; -g_val; 0; 0];

% input matrix g(X):
%   row 4 (dvx/dt):   [ cos(theta)/m,  sin(theta)/m ]
%   row 5 (dvy/dt):   [-sin(theta)/m,  cos(theta)/m ]
%   row 6 (domega/dt):[ -6/(m*L),      0            ]
%   row 7 (dm/dt):    [ 0,             -1/(Isp*g)   ]
gx_sym = [0, 0;
          0, 0;
          0, 0;
          cos(x3) / x7, sin(x3) / x7;
          -sin(x3) / x7, cos(x3) / x7;
          -6 / (x7 * L_val), 0;
          0, -1 / (Isp * g_val)];

% ── Output map h(x) via LQR on the linearised 6-state sub-system ─────────────
% The mass x7 does not appear in the linear model, so we linearise the
% 6D sub-system X6 = [x; y; theta; vx; vy; omega] around hover equilibrium:
%   X_eq = 0,  U_eq = [0; m_f * g]
m_f = 85000.0; % linearisation mass (dry mass at end-of-burn, kg)

% Linearised A (6x6)
A_lin = zeros(6, 6);
A_lin(1, 4) = 1.0; % dx/dt    = vx
A_lin(2, 5) = 1.0; % dy/dt    = vy
A_lin(3, 6) = 1.0; % dtheta/dt = omega
A_lin(4, 3) = g_val; % dvx/dt: gravity-coupling through tilt angle theta

% Linearised B (6x2)
B_lin = zeros(6, 2);
B_lin(4, 1) = 1.0 / m_f;
B_lin(5, 2) = 1.0 / m_f;
B_lin(6, 1) = -6.0 / (m_f * L_val);

% LQR weights (matching example_02_lqr.ipynb)
Q_lqr = diag([100, 100, 500, 10, 10, 50]);
R_lqr = diag([1e-6, 1e-6]);

% Solve continuous-time LQR: K = R^{-1} B^T P   (2 x 6 matrix)
[K_lqr, ~, ~] = lqr(A_lin, B_lin, Q_lqr, R_lqr);

fprintf('LQR output matrix C (2x6):\n');
disp(K_lqr);

% Output map: y = C * [x1; x2; x3; x4; x5; x6]  (virtual output for backstepping)
hx_sym = sym(K_lqr) * [x1; x2; x3; x4; x5; x6];
syms y1 y2;
y_vars_sym = [y1; y2];

% define the safe set and target set in the output space

% =========================================================================
% Find the largest safe polynomial superlevel set inscribed in the
% physical state constraints using SOSTOOLS (S-procedure / SOS relaxation)
%
% Safe set candidate: h_safe(X) = rho - (Y1^4 + Y2^4) >= 0
% where Y1, Y2 are the virtual outputs from the LQR output map.
%
% Containment certificate (S-procedure):
%   For each physical constraint p_i(X) >= 0, find SOS multiplier sigma_i
%   such that  p_i(X) - sigma_i(X) * h_safe(X)  is SOS.
%   This proves: { X : h_safe(X) >= 0 } ⊆ { X : p_i(X) >= 0 }.
% =========================================================================

% 1. Declare polynomial indeterminates for SOSTOOLS
%    Working in the 6D sub-state X6 = [x; y; theta; vx; vy; omega]
%    (mass m is not part of the LQR output map)
pvar x y theta vx vy omega;
X_sos = [x; y; theta; vx; vy; omega];

% 2. Physical constraints on the 6D sub-state
theta_max = 30 * pi / 180; % 30-deg tilt limit
p_ground = y; % altitude > 0
p_angle = theta_max ^ 2 - theta ^ 2; % attitude within 30 deg of vertical
p_glide = y ^ 2 - x ^ 2; % 45-deg glide-slope cone  (y > |x|)
p_vel = 100 ^ 2 - vy ^ 2; % vertical speed <= 100 m/s

PhysConstraints = {p_ground, p_angle, p_glide, p_vel};
nc = length(PhysConstraints);

% 3. Virtual outputs using LQR gain K_lqr (2x6), already computed above
Y1 = K_lqr(1, :) * X_sos;
Y2 = K_lqr(2, :) * X_sos;

% 4. Monomial basis for degree-2 SOS multipliers (degrees 0, 1, 2)
monos = monomials(X_sos, 0:2);
% monom2 = [monomials(X_sos, 0); monomials(X_sos, 1); monomials(X_sos, 2)];

% 5. Binary search for the largest rho such that the containment holds
rho_low = 0; rho_high = 200; rho_best = 0;

for iter = 1:10
    rho_test = (rho_low + rho_high) / 2;

    % Candidate safe set superlevel set: h_safe(X) = rho - (Y1^4 + Y2^4)
    h_safe = rho_test - (Y1 ^ 3 + Y2 ^ 3); % using degree 3 safe set for better fit than degree 4
    % h_safe = rho_test - (Y1 ^ 4 + Y2 ^ 4); % using degree 4 safe set for better fit than degree 3

    % Initialise SOSTOOLS feasibility program
    prog = sosprogram(X_sos);

    % S-procedure: add one SOS multiplier per physical constraint
    for i = 1:nc
        % sosvar declares sigma_i as an SOS polynomial (guaranteed >= 0)
        [prog, sigma_i] = sospolyvar(prog, monos);
        % Require: p_i(X) - sigma_i(X) * h_safe(X) is SOS
        prog = sosineq(prog, PhysConstraints{i} - sigma_i * h_safe);
    end

    % Solve the SOS feasibility program
    solver_opt.solver = 'mosek';
    % turn off verbose output
    solver_opt.verbose = 0;
    prog = sossolve(prog, solver_opt);

    % feasratio > 0  → feasible (certificate found)
    % feasratio < 0  → infeasible at this rho
    if prog.solinfo.info.feasratio > 0
        rho_best = rho_test;
        rho_low = rho_test;
    else
        rho_high = rho_test;
    end

end

fprintf('Largest inscribed safe region parameter: rho = %f\n', rho_best);
fprintf('Safe set: { X6 : %f - (Y1^4 + Y2^4) >= 0 }\n', rho_best);
