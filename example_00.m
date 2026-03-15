% demon script to find the reach-avoid backstepping controller with bounded control inputs
% using scenario optimization programming (SOP) and sum-of-squares (SOS) programming

clc;
clear;
close all;

% define the system dynamics using symbolic variables
% Double integrator:  dx1/dt = x2,  dx2/dt = u1
syms x1 x2 y1 y2;
x_vars_sym = [x1; x2];

% drift term f(x) — no drift for a pure double integrator
fx_sym = [x2; 0];

% input matrix g(x) — control enters only the acceleration channel
gx_sym = [0 1; 1 0];

hx_sym = [x1; x2]; % output of the system, also the state for the single-integrator system in the backstepping design

y_vars_sym = [y1; y2];

% define the safe set
safe_set_sym = 1 - y1 ^ 2 - y2 ^ 2; % safe_set: safe_set >=0 inside safe set
% define the target set
% target_set_sym = y1 ^ 2 + y2 ^ 2 - 0.01; % target_set: target_set <=0 inside target set
target_set_sym = (2 * ((y2 - 0.1) / 2) ^ 2 + 3 * ((y1 + 0.4) / 3) ^ 4 ...
    + (3 * y1 + 0.3) ^ 2 * (4 * y2 + 0.2) ^ 2 - 0.01);

% synthesize the reach-avoid backstepping controller with bounded control inputs (symbolic)
[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym);
% here u is a symbolic vector [u1; u2], each entry is a symbolic expression of state variables x1,x2,x3,x4, and unknown parameters
% ers (mu values, lambda, k1 controler polynomial)

% define constraint for the control input bounds
% Ax = [1 0; -1 0; 0 1; 0 -1]; % example constraint matrix for control input bounds, here we want to enforce |u1| <= 100 and |u2| <= 100
lb = [-1; -1]; % lower bounds for control inputs
ub = [1; 1]; % upper bounds for control inputs
ds = 8; % degree of the auxiliary SOS polynomials for the single-integrator system
dv = 8; % degree of the k1 controller polynomial

mu_val = 0.001; % example value for mu just for testing, HAVE TO DISCUSS THIS IN THE PAPER

samples_num = 100; % number of random samples to find the valid samples that satisfy the control input bounds for the pseudo ux

bound_min = [-1.1; -1.1]; % lower bounds for sampling the state space for finding valid samples that satisfy the set (safe, target, vanilla reach-avoid certificate) constraints
bound_max = [1.1; 1.1]; % upper bounds for sampling the state space for finding valid samples that satisfy the set (safe, target, vanilla reach-avoid certificate) constraints

% solve the bounded control inputs using scenario optimization programming (SOP) with SOS constraints
[u_opt, certificate_opt, valid_count, k1_opt] = solvesop_bounded_control(u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, x_vars_sym, y_vars_sym, ...
    hx_sym, safe_set_sym, target_set_sym, mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max);

% disp the obtained controller after solving with bounded control inputs
disp('Obtained controller after solving with bounded control inputs:');
disp(u_opt);

% disp the obtained certificate after solving with bounded control inputs
disp('Obtained certificate after solving with bounded control inputs:');
disp(certificate_opt);

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

export_to_python(u_opt, certificate_opt, k1_opt, params_for_export, 'sop_bounded_control_ex1_debug.py');
% export the computed controller and certificate to a python file for verification and testing
