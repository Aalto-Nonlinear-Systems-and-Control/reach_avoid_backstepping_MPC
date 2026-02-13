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

% before solving the bounded control, we have to fix parameter values (mu and lambda) with predefined numerical values
mu_val = 1; % example value for mu just for testing, HAVE TO DISCUSS THIS IN THE PAPER
lambda_val = 1e-8; % example value for lambda as a small positive value, HAVE TO DISCUSS THIS IN THE PAPER

% substitute mu values and lambda with predefined numerical values in the controller u
u_subs = u;

for i = 1:numel(mu)

    if ~isempty(mu{i})
        u_subs = subs(u_subs, mu{i}, mu_val);
    end

end

u_subs = subs(u_subs, lambda, lambda_val);

% substitute mu values and lambda with predefined numerical values in the certificate
certificate_subs = certificate;

for i = 1:numel(mu)

    if ~isempty(mu{i})
        certificate_subs = subs(certificate_subs, mu{i}, mu_val);
    end

end

certificate_subs = subs(certificate_subs, lambda, lambda_val);

% substitute the output mapping

% check all variables in the substituted controller expression
vars_in_u_subs = symvar(u_subs);
disp('Variables in the substituted controller expression:');
disp(vars_in_u_subs);

% print the obtained certificate
disp('Obtained certificate:');
disp(certificate_subs);

% print the obtained controller
disp('Obtained controller after substituting mu and lambda with predefined values:');
disp(u_subs);

% solve the bounded control inputs using scenario optimization programming (SOP) with SOS constraints
solvesop_bounded_control(u_subs, k1, certificate_subs, x_vars_sym, y_vars_sym, hx_sym, safe_set_sym, target_set_sym, -100, 100);
