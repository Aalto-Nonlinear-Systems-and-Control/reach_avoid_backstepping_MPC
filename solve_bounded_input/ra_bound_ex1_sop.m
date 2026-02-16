% demon script to find the reach-avoid backstepping controller with bounded control inputs
% using scenario optimization programming (SOP) and sum-of-squares (SOS) programming

clc;
clear;
close all;

% define the system dynamics using symbolic variables
syms x1 x2 x3 x4 u1 u2 y1 y2;
x_vars_sym = [x1; x2; x3; x4];

fx_sym = [x2;
          - (x2 ^ 2) + x3 ^ 3 + x4 - 5 * x1;
          x4;
          2 * x2 * x3 + 2 * x3 * x4 - 7 * x1];

gx_sym = [0 0;
          2 * x2 ^ 2 + 5 * x3 * x4 - 9 * x1, 3 * x3 ^ 3 + 2 * x3 * x4 - 10 * x1;
          0 0;
          2 * x2 * x3 ^ 3 + 2 * x3 * x4 + 7 * x2, 7 * x2 ^ 2 * x3 + 2 * x3 * x4];

hx_sym = [x1; x3]; % output of the system, also the state for the single-integrator system in the backstepping design

y_vars_sym = [y1; y2];

% define the safe set
safe_set_sym = 1 - y1 ^ 2 - y2 ^ 2; % safe_set: safe_set >=0 inside safe set
% define the target set
target_set_sym = y1 ^ 2 + y2 ^ 2 - 0.01; % target_set: target_set <=0 inside target set

% synthesize the reach-avoid backstepping controller with bounded control inputs (symbolic)
[u, k1, mu, lambda, certificate, A_matrix, b_vector, ks] = reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, safe_set_sym);
% here u is a symbolic vector [u1; u2], each entry is a symbolic expression of state variables x1,x2,x3,x4, and unknown parameters (mu values, lambda, k1 controler polynomial)

% define constraint for the control input bounds
% Ax = [1 0; -1 0; 0 1; 0 -1]; % example constraint matrix for control input bounds, here we want to enforce |u1| <= 100 and |u2| <= 100
lb = [-100; -100]; % lower bounds for control inputs
ub = [100; 100]; % upper bounds for control inputs
ds = 2; % degree of the auxiliary SOS polynomials for the single-integrator system
dv = 2; % degree of the k1 controller polynomial

mu_val = 1; % example value for mu just for testing, HAVE TO DISCUSS THIS IN THE PAPER

samples_num = 10; % number of random samples to find the valid samples that satisfy the control input bounds for the pseudo ux

% solve the bounded control inputs using scenario optimization programming (SOP) with SOS constraints
[u_opt, certificate_opt] = solvesop_bounded_control(u, k1, mu, lambda, certificate, x_vars_sym, y_vars_sym, hx_sym, safe_set_sym, target_set_sym, mu_val, lb, ub, ds, dv, samples_num);

% disp the obtained controller after solving with bounded control inputs
disp('Obtained controller after solving with bounded control inputs:');
disp(u_opt);

% disp the obtained certificate after solving with bounded control inputs
disp('Obtained certificate after solving with bounded control inputs:');
disp(certificate_opt);
