% demon script to find the reach-avoid backstepping controller with bounded control inputs
% using scenario optimization programming (SOP) and sum-of-squares (SOS) programming

clc;
clear;
close all;

% define the system dynamics using symbolic variables
% Planar rocket:  state = [x; z; theta; xdot; zdot; thetadot], input = [F; M]
%
%   xdot      = [xdot; zdot; thetadot; 0; -g; 0]
%             + [0 0; 0 0; 0 0; sin(theta)/m 0; cos(theta)/m 0; 0 -1/I] * [F; M]
%
syms px pz th vx vz om y1 y2;
x_vars_sym = [px; pz; th; vx; vz; om];

% physical parameters (SI units)
g_grav = 9.81; % gravitational acceleration [m/s^2]
m_mass = 1.0; % vehicle mass [kg]
I_inertia = 0.1; % moment of inertia [kg*m^2]

% drift term f(x)
fx_sym = [vx; vz; om; 0; -g_grav; 0];

% input matrix g(x)  (6 x 2)
gx_sym = [ ...
              0, 0; ...
              0, 0; ...
              0, 0; ...
              sin(th) / m_mass, 0; ...
              cos(th) / m_mass, 0; ...
              0, -1 / I_inertia ...
          ];

% geometry parameter
L = 1.0; % pendulum/nozzle arm length [m]

% output: tip position  y = [x - L*theta; z - L]
hx_sym = [px - L * sin(th); pz - L * cos(th)];

y_vars_sym = [y1; y2];

% define the safe set: tilted squircle — same boundary shape/position, peak at z_peak
% psi = squircle(y) * (1 + beta*(y2-z_c))
%   boundary : zero-level set identical to plain squircle (shape/position unchanged)
%   peak     : analytically placed at (x_c, z_peak) by computing beta from z_peak
%   formula  : beta = 4*C*u^3 / (W_x^4 - 5*C*u^4), u = z_peak-z_c, C = (W_x/W_z)^4
%   validity : |z_peak - z_c| < W_z * (1/5)^0.25  ≈  0.669 * W_z
W_x = 2; % working-space half-width in y1 direction [m]
W_z = 2; % working-space half-width in y2 direction [m]
x_c = 0.0; % squircle centre in y1 [m]
z_c = 1; % squircle centre in y2 [m]
C_coeff = (W_x / W_z) ^ 4;
safe_set_sym = (W_x ^ 4 - (y1 - x_c) ^ 4 - C_coeff * (y2 - z_c) ^ 4);

% define the target set: flat ellipse — tight in y2, relaxed in y1
% target position
x_target = 0.0; % target y1 [m]  (tip x at landing: px ≈ 0)
z_target = 1.0; % target y2 [m]  (tip z at landing: pz = z_target+L = 0.0)
r_y1 = 0.3; % semi-axis in y1 [m]  (larger → more relaxed horizontally)
r_y2 = 0.05; % semi-axis in y2 [m]  (smaller → tighter vertically, flat shape)
target_set_sym = (y1 - x_target) ^ 2 / r_y1 ^ 2 + (y2 - z_target) ^ 2 / r_y2 ^ 2 - 1; % target_set <= 0 inside ellipse

% synthesize the reach-avoid backstepping controller with bounded control inputs (symbolic)
[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym);
% here u is a symbolic vector [u1; u2], each entry is a symbolic expression of state variables x1,x2,x3,x4, and unknown parameters
% ers (mu values, lambda, k1 controler polynomial)

% define constraint for the control input bounds
% Ax = [1 0; -1 0; 0 1; 0 -1]; % example constraint matrix for control input bounds, here we want to enforce |u1| <= 100 and |u2| <= 100
% set the random seed for reproducibility
rng(42);
lb = [-80; -50]; % lower bounds for control inputs [F_min; M_min]  (widened from observed u1_min=-115, u2_min=-20)
ub = [350; 50]; % upper bounds for control inputs [F_max; M_max]  (widened from observed u1_max=23, u2_max=15)
ds = 4; % degree of the auxiliary SOS polynomials for the single-integrator system
dv = 4; % degree of the k1 controller polynomial

mu_val = 10; % example value for mu just for testing, HAVE TO DISCUSS THIS IN THE PAPER

samples_num = 1000; % number of random samples to find the valid samples that satisfy the control input bounds for the pseudo ux

% state order: [px; pz; th; vx; vz; om]
bound_min = [-2.0; 0; -pi / 6; -0.1; -0.1; -0.1]; % state order: [px; pz; th; vx; vz; om]  (pz centered on landing approach)
bound_max = [2.0; 4; pi / 6; 0.1; 0.1; 0.1]; % (pz_max=4 → y2_max=1.0 covers near-landing region)

% solve the bounded control inputs using scenario optimization programming (SOP) with SOS constraints
[u_opt, certificate_opt, valid_count] = solvesop_bounded_control(u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, x_vars_sym, y_vars_sym, ...
    hx_sym, safe_set_sym, target_set_sym, mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max);

% % compute the bounds of the obtained controller over zero superlevel set of the certificate
% [num_1, den_1] = numden(u_opt(1));
% % disp the obtained controller polynomial numerator and denominator for debugging
% disp('Numerator of the obtained controller u1:');
% disp(num_1);
% disp('Denominator of the obtained controller u1:');
% disp(den_1);
% [lb_1, ub_1] = compute_poly_bounds_sos(num_1, den_1, certificate_opt, ds, 1e-4);
% estimate the bounds of the obtained controller over zero superlevel set of the certificate using sampling (for verification)
[estimated_lb_1, estimated_ub_1] = compute_poly_bounds_sampling(x_vars_sym, u_opt(1), certificate_opt, 50000, bound_min, bound_max);
[estimated_lb_2, estimated_ub_2] = compute_poly_bounds_sampling(x_vars_sym, u_opt(2), certificate_opt, 50000, bound_min, bound_max);

disp('------------------------------------------------------------------------------------');

% disp the obtained controller after solving with bounded control inputs
disp('Obtained controller after solving with bounded control inputs:');
disp(u_opt);

% disp the obtained certificate after solving with bounded control inputs
disp('Obtained certificate after solving with bounded control inputs:');
disp(certificate_opt);

% disp the number of valid samples that satisfy the control input bounds for the pseudo ux
disp(['Number of valid samples that satisfy the control input bounds for the pseudo ux: ', num2str(valid_count), ' out of ', num2str(samples_num)]);

% disp(['Bounds for u1 over the zero superlevel set of the certificate: [', ...
%           num2str(lb_1), ', ', num2str(ub_1), ']']);

disp(['Estimated bounds for u1 over the zero superlevel set of the certificate using sampling: [', ...
          num2str(estimated_lb_1), ', ', num2str(estimated_ub_1), ']']);
disp(['Estimated bounds for u2 over the zero superlevel set of the certificate using sampling: [', ...
          num2str(estimated_lb_2), ', ', num2str(estimated_ub_2), ']']);

% export the results, also all setting parameters, to a python file for verification and testing
% first, use a struct to collect all the relevant parameters for exporting
params_for_export = struct();
params_for_export.fx_sym = fx_sym;
params_for_export.gx_sym = gx_sym;
params_for_export.hx_sym = hx_sym;
params_for_export.x_vars_sym = x_vars_sym;
params_for_export.y_vars_sym = y_vars_sym;
params_for_export.safe_set_sym = safe_set_sym;
params_for_export.target_set_sym = target_set_sym;
params_for_export.lb = lb;
params_for_export.ub = ub;
params_for_export.ds = ds;
params_for_export.dv = dv;
params_for_export.mu_val = mu_val;
params_for_export.samples_num = samples_num;
params_for_export.valid_count = valid_count;
params_for_export.bound_min = bound_min;
params_for_export.bound_max = bound_max;

export_to_python(u_opt, certificate_opt, params_for_export, 'sop_bounded_control_ex1_debug.py');
% export the computed controller and certificate to a python file for verification and testing
