% demon script to find the reach-avoid backstepping controller with bounded control inputs

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

% check all variables in the substituted controller expression
vars_in_u_subs = symvar(u_subs);
disp('Variables in the substituted controller expression:');
disp(vars_in_u_subs);

% solve for the reach-avoid backstepping controller with bounded control inputs
% here we only consider A_matrix is a constant identity matrix for simplicity
A_matrix = eye(2); % 2x2 identity matrix
% for simplicity, we set the bounds as constants like below
u_lb = [-100; -100]; % lower bound for each control input
u_ub = [100; 100]; % upper bound for each control input
u_opt = solvesos_bounded_control(u_subs, x_vars_sym, y_vars_sym, k1, hx_sym, safe_set_sym, target_set_sym, A_matrix, u_lb, u_ub);

% exit the script here for debugging
return;

% disp the synthesized controller (symbolic)
disp('The synthesized reach-avoid backstepping controller (symbolic):');
disp(u);
% disp the shape of the controller
disp('The shape of the controller:');
disp(size(u));

% determine the nominator and denominator of each control input
[num_u1, den_u1] = numden(u(1));
[num_u2, den_u2] = numden(u(2));
disp('Nominator and Denominator of u1:');
disp('Nominator of u1:');
disp(num_u1);
disp('Denominator of u1:');
disp(den_u1);
disp('Nominator and Denominator of u2:');
disp('Nominator of u2:');
disp(num_u2);
disp('Denominator of u2:');
disp(den_u2);

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% key step for considering the control input bounds in actual controller u(x)
% define the variables for solving the feasible k1 controller using SOS programming
pvar x1_k1 x2_k1 x3_k1 x4_k1 y1_k1 y2_k1;
dpvar delta_k1 lambda_k1;
vars_k1 = [x1_k1; x2_k1; x3_k1; x4_k1]; % state variables for k1 controller synthesis

xi = 1e-8;

ds = 6; % degree of auxiliary polynomials
d_k1 = 3; % degree of control polynomials

monos_k1 = monomials(vars_k1, 0:d_k1);

prog_k1 = sosprogram(vars_k1, [delta_k1, lambda_k1]);

[prog_k1, k11] = sospolyvar(prog_k1, monos_k1);
[prog_k1, k12] = sospolyvar(prog_k1, monos_k1);

% transform symbolic expressions to pvar polynomials for SOS programming
safe_set_pvar = sym2pvar(safe_set_sym, vars_sym, vars_k1);
target_set_pvar = sym2pvar(target_set_sym, vars_sym, vars_k1);
hx_pvar = [sym2pvar(hx_sym(1), vars_sym, vars_k1); sym2pvar(hx_sym(2), vars_sym, vars_k1)];

% compute Lfh for the single-integrator system \dot{y} = k1
Lfh_k1 = diff(safe_set_pvar, hx_pvar(1)) * k11 + diff(safe_set_pvar, hx_pvar(2)) * k12;

% auxiliary polynomials
monos_s_k1 = monomials(vars_k1, 0:ds);
[prog_k1, s0_k1] = sospolyvar(prog_k1, monos_s_k1);
[prog_k1, s1_k1] = sospolyvar(prog_k1, monos_s_k1);

% constraints
prog_k1 = sosineq(prog_k1, Lfh_k1 - lambda_k1 * safe_set_pvar + delta_k1 - s0_k1 * safe_set_pvar - s1_k1 * target_set_pvar);
prog_k1 = sosineq(prog_k1, delta_k1);
prog_k1 = sosineq(prog_k1, lambda_k1 - xi);
prog_k1 = sosineq(prog_k1, s0_k1);
prog_k1 = sosineq(prog_k1, s1_k1);

% additional constraints to enforce control input bounds >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% convert the reach_avoid backstepping controller to pvar expressions
% Note: num_u1, den_u1, etc. contain state variables AND symbolic k1 controller variables
% We need to map both sets of variables

% Build the complete variable mapping
% State variables: x1, x2, x3, x4 -> x1_k1, x2_k1, x3_k1, x4_k1
% k1 controller: k1(1), k1(2) -> k11, k12 (decision variables)
% define k11_pvar and k12_pvar as pvar variables for substitution
dpvar k11_pvar k12_pvar;
all_sym_vars = [vars_sym; k1]; % [x1; x2; x3; x4; k1_1; k1_2]
all_pvar_vars = [vars_k1; k11_pvar; k12_pvar]; % [x1_k1; x2_k1; x3_k1; x4_k1; k11_pvar; k12_pvar]

% Predefined parameter values for mu and lambda
mu_val = 1; % mu parameter value
lambda_val = 1e-5; % lambda parameter value (use a small positive value)

% Substitute mu values and lambda with predefined numerical values
num_u1_subs = num_u1;
num_u2_subs = num_u2;
den_u1_subs = den_u1;
den_u2_subs = den_u2;

% Substitute mu values
for i = 1:length(mu)

    if ~isempty(mu{i})

        for j = 1:length(mu{i})
            num_u1_subs = subs(num_u1_subs, mu{i}(j), mu_val);
            num_u2_subs = subs(num_u2_subs, mu{i}(j), mu_val);
            % den_u1_subs = subs(den_u1_subs, mu{i}(j), mu_val); % no need to substitute mu in denominator
            % den_u2_subs = subs(den_u2_subs, mu{i}(j), mu_val); % no need to substitute mu in denominator
        end

    end

end

% Substitute lambda with predefined value
num_u1_subs = subs(num_u1_subs, lambda, lambda_val);
num_u2_subs = subs(num_u2_subs, lambda, lambda_val);
% den_u1_subs = subs(den_u1_subs, lambda, lambda_val); % no need to substitute lambda in denominator
% den_u2_subs = subs(den_u2_subs, lambda, lambda_val); % no need to substitute lambda in denominator

% Now convert to pvar (all parameters are now numerical)
den_u1_pvar = sym2pvar(den_u1_subs, all_sym_vars, all_pvar_vars);
num_u1_pvar = sym2pvar(num_u1_subs, all_sym_vars, all_pvar_vars);
den_u2_pvar = sym2pvar(den_u2_subs, all_sym_vars, all_pvar_vars);
num_u2_pvar = sym2pvar(num_u2_subs, all_sym_vars, all_pvar_vars);

disp('Converted controller numerators and denominators to pvar expressions.');

% Debug: check if k11_pvar and k12_pvar appear in the expressions
disp('Variables in num_u1_pvar:');
disp(num_u1_pvar.varname);
disp('Variables in den_u1_pvar:');
disp(den_u1_pvar.varname);

% Check if the denominators actually contain k11_pvar or k12_pvar
% If not, they don't need substitution
has_k11_in_num1 = ismember('k11_pvar', num_u1_pvar.varname);
has_k12_in_num1 = ismember('k12_pvar', num_u1_pvar.varname);
has_k11_in_den1 = ismember('k11_pvar', den_u1_pvar.varname);
has_k12_in_den1 = ismember('k12_pvar', den_u1_pvar.varname);

disp(['num_u1 contains k11_pvar: ', num2str(has_k11_in_num1)]);
disp(['num_u1 contains k12_pvar: ', num2str(has_k12_in_num1)]);
disp(['den_u1 contains k11_pvar: ', num2str(has_k11_in_den1)]);
disp(['den_u1 contains k12_pvar: ', num2str(has_k12_in_den1)]);

% For now, skip the substitution step and just use the pvar expressions directly
% The SOS constraints will be formulated differently
disp('Skipping polynomial substitution - will formulate SOS constraints directly.');

% TODO: substitute k11 and k12 into the controller u expressions ???????

% % Substitute the k11_pvar and k12_pvar in num_u1_pvar and num_u2_pvar with constructed k11 and k12
% % Use polynomial composition since k11 and k12 are polynomial expressions with dpvar coefficients
% % For polynomial objects, we need to substitute one variable at a time
% num_u1_final = subs(num_u1_pvar, k11_pvar, k11);
% num_u1_final = subs(num_u1_final, k12_pvar, k12);
%
% num_u2_final = subs(num_u2_pvar, k11_pvar, k11);
% num_u2_final = subs(num_u2_final, k12_pvar, k12);
%
% den_u1_final = subs(den_u1_pvar, k11_pvar, k11);
% den_u1_final = subs(den_u1_final, k12_pvar, k12);
%
% den_u2_final = subs(den_u2_pvar, k11_pvar, k11);
% den_u2_final = subs(den_u2_final, k12_pvar, k12);
%
% disp('Substituted k1 controller polynomials into controller expressions.');

% set objective to minimize delta_k1
prog_k1 = sossetobj(prog_k1, delta_k1);
% solve the SOS program
solver_opt_k1.solver = 'mosek';
prog_k1 = sossolve(prog_k1, solver_opt_k1);
k1_1_value = sosgetsol(prog_k1, k11);
k1_2_value = sosgetsol(prog_k1, k12);

% get the solved lambda_k1 and delta_k1 (convert to numerical values)
lambda_k1_value = double(sosgetsol(prog_k1, lambda_k1));
delta_k1_value = double(sosgetsol(prog_k1, delta_k1));

% disp the solved lambda_k1 and delta_k1
disp(['The solved lambda_k1: ', num2str(lambda_k1_value)]);
disp(['The solved delta_k1: ', num2str(delta_k1_value)]);
disp('The synthesized k1 controller polynomials:');
disp('k1_1(x1,x2,x3,x4) = ');
disp(k1_1_value);
disp('k1_2(x1,x2,x3,x4) = ');
disp(k1_2_value);

% =====================================================
% Convert solved k1 controller to symbolic for simulation
% =====================================================
% convert pvar polynomials to symbolic expressions
k1_1_sym = pvar2sym(k1_1_value, vars_k1, vars_sym);
k1_2_sym = pvar2sym(k1_2_value, vars_k1, vars_sym);

disp('k1_1 as symbolic expression:');
disp(k1_1_sym);
disp('k1_2 as symbolic expression:');
disp(k1_2_sym);

% Substitute k1 controller into the reach-avoid backstepping controller u
u1_subs = subs(u(1), k1, [k1_1_sym; k1_2_sym]);
u2_subs = subs(u(2), k1, [k1_1_sym; k1_2_sym]);

disp('The reach-avoid backstepping controller after substituting k1 controller:');
disp('u1(x) = ');
disp(u1_subs);
disp('u2(x) = ');
disp(u2_subs);

% Construct the certificate with k1 substituted
certificate_subs = subs(certificate, k1, [k1_1_sym; k1_2_sym]);

% substitute mu values as mu_val and lambda as lambda_k1_value
certificate_final = certificate_subs;

for i = 1:length(mu)

    if ~isempty(mu{i})

        for j = 1:length(mu{i})
            certificate_final = subs(certificate_final, mu{i}(j), mu_val);
        end

    end

end

certificate_final = subs(certificate_final, lambda, lambda_k1_value);

disp('Final Certificate after substituting all parameters:');
disp(certificate_final);

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% =====================================================
% % solve an k1 controller, and substitute it into the controller u for debugging purpose
% pvar y1 y2;
% dpvar delta_k1 lambda_k1;
% vars_k1 = [y1; y2]; % state variables for k1 controller synthesis

% xi = 1e-8;

% ds = 8; % degree of auxiliary polynomials
% d_k1 = 3; % degree of control polynomials

% monos_k1 = monomials(vars_k1, 0:d_k1);

% prog_k1 = sosprogram(vars_k1, [lambda_k1, delta_k1]);

% [prog_k1, k11] = sospolyvar(prog_k1, monos_k1);
% [prog_k1, k12] = sospolyvar(prog_k1, monos_k1);

% % safe set and target set for k1 controller synthesis
% h_k1 = 1 - y1 ^ 2 - y2 ^ 2; % h(x)>0
% g_k1 = y1 ^ 2 + y2 ^ 2 -0.01; % g(x)<0

% lhs_k1 = diff(h_k1, vars_k1(1)) * k11 + diff(h_k1, vars_k1(2)) * k12;

% % auxiliary polynomials
% monos_s_k1 = monomials(vars_k1, 0:ds);
% [prog_k1, s0_k1] = sospolyvar(prog_k1, monos_s_k1);
% [prog_k1, s1_k1] = sospolyvar(prog_k1, monos_s_k1);

% % constraints
% prog_k1 = sosineq(prog_k1, lhs_k1 - lambda_k1 * h_k1 + delta_k1 - s0_k1 * h_k1 - s1_k1 * g_k1);
% prog_k1 = sosineq(prog_k1, delta_k1);
% prog_k1 = sosineq(prog_k1, lambda_k1 - xi);
% prog_k1 = sosineq(prog_k1, s0_k1);
% prog_k1 = sosineq(prog_k1, s1_k1);

% % set objective to minimize delta_k1
% prog_k1 = sossetobj(prog_k1, delta_k1);
% % solve the SOS program
% solver_opt_k1.solver = 'mosek';
% prog_k1 = sossolve(prog_k1, solver_opt_k1);

% k1_1_value = sosgetsol(prog_k1, k11);
% k1_2_value = sosgetsol(prog_k1, k12);
% % get the solved lambda_k1 and delta_k1 (convert to numerical values)
% lambda_k1_value = double(sosgetsol(prog_k1, lambda_k1));
% delta_k1_value = double(sosgetsol(prog_k1, delta_k1));

% % disp the solved lambda_k1 and delta_k1
% disp(['The solved lambda_k1: ', num2str(lambda_k1_value)]);
% disp(['The solved delta_k1: ', num2str(delta_k1_value)]);

% disp('The synthesized k1 controller polynomials:');
% disp('k1_1(y1,y2) = ');
% disp(k1_1_value);
% disp('k1_2(y1,y2) = ');
% disp(k1_2_value);

% % convert pvar polynomials to symbolic expressions
% syms y1_sym y2_sym real;
% k1_1_sym = pvar2sym(k1_1_value, [y1; y2], [y1_sym; y2_sym]);
% k1_2_sym = pvar2sym(k1_2_value, [y1; y2], [y1_sym; y2_sym]);

% disp('k1_1 as symbolic expression:');
% disp(k1_1_sym);
% disp('k1_2 as symbolic expression:');
% disp(k1_2_sym);

% % substitute the y=h(x) into k1 polynomials
% k1_1_subs = subs(k1_1_sym, [y1_sym, y2_sym], [hx_sym(1), hx_sym(2)]);
% k1_2_subs = subs(k1_2_sym, [y1_sym, y2_sym], [hx_sym(1), hx_sym(2)]);

% disp('k1_1 substituted (as function of x):');
% disp(k1_1_subs);
% disp('k1_2 substituted (as function of x):');
% disp(k1_2_subs);

% % DEBUG: substitute the k1 controller without considering control input bounds into the reach-avoid backstepping controller u
% % this is just for verify the reach-avoid backstepping controller synthesis works correctly
% u1_subs = subs(u(1), k1, [k1_1_subs; k1_2_subs]);
% u2_subs = subs(u(2), k1, [k1_1_subs; k1_2_subs]);

% disp('The reach-avoid backstepping controller after substituting k1 controller:');
% disp('u1(x) = ');
% disp(u1_subs);
% disp('u2(x) = ');
% disp(u2_subs);

% % substitute the reach-avoid backstepping controller to the system dynamics
% fx_closed_loop = fx_sym + gx_sym * [u1_subs; u2_subs];
% disp('The closed-loop system dynamics after substituting the reach-avoid backstepping controller:');
% disp('f_cl(x) = ');
% disp(fx_closed_loop);

% % Display and construct the reach-avoid backstepping certificate
% disp('=== Reach-Avoid Certificate (symbolic) ===');
% disp('Certificate (with symbolic parameters):');
% disp(certificate);

% % Substitute k1 values into the certificate
% certificate_subs = subs(certificate, k1, [k1_1_subs; k1_2_subs]);
% disp('Certificate after substituting k1 controller:');
% disp(certificate_subs);

% % substitute mu values as 1 for all mu parameters
% % mu is a cell array, so we need to iterate over each output's mu values
% certificate_final = certificate_subs;

% for i = 1:length(mu)

%     if ~isempty(mu{i})

%         for j = 1:length(mu{i})
%             certificate_final = subs(certificate_final, mu{i}(j), 1);
%         end

%     end

% end

% disp('Final Certificate after substituting mu values as 1:');
% disp(certificate_final);

% =====================================================
% SIMULATION: verify the reach-avoid backstepping controller via simulation

% now random sample points and keep only those be verified by the certificate
% Points must satisfy: certificate >= 0, safe_set >= 0, AND outside target set
num_samples = 50;
safe_points = [];

% Create lambdified functions for faster evaluation
safe_set_func = matlabFunction(safe_set_sym, 'Vars', [x1, x2, x3, x4]);
target_set_func = matlabFunction(target_set_sym, 'Vars', [x1, x2, x3, x4]);
cert_func = matlabFunction(certificate_final, 'Vars', [x1, x2, x3, x4]);

while size(safe_points, 1) < num_samples
    % sample a random point in the state space (same range as Python: [-1, 1])
    x_sample = -1 + 2 * rand(1, 4); % sample in [-1,1] for each state variable

    % Check all three conditions:
    % 1. Certificate >= 0
    cert_val = cert_func(x_sample(1), x_sample(2), x_sample(3), x_sample(4));
    % 2. Safe set >= 0 (inside safe set)
    safe_val = safe_set_func(x_sample(1), x_sample(2), x_sample(3), x_sample(4));
    % 3. Target set > 0 (outside target set, since target_set <= 0 is inside)
    target_val = target_set_func(x_sample(1), x_sample(2), x_sample(3), x_sample(4));

    if cert_val >= 0 && safe_val >= 0 && target_val > 0
        safe_points = [safe_points; x_sample]; % add to safe points
    end

end

disp(['Number of safe points verified by the certificate: ', num2str(size(safe_points, 1))]);

% Prepare control inputs with all parameters substituted
u1_final = u1_subs;
u2_final = u2_subs;
% substitute mu values as mu_val and lambda value
for i = 1:length(mu)

    if ~isempty(mu{i})

        for j = 1:length(mu{i})
            u1_final = subs(u1_final, mu{i}(j), mu_val);
            u2_final = subs(u2_final, mu{i}(j), mu_val);
        end

    end

end

u1_final = subs(u1_final, lambda, lambda_k1_value);
u2_final = subs(u2_final, lambda, lambda_k1_value);

% Convert symbolic expressions to MATLAB functions for faster evaluation
u1_func = matlabFunction(u1_final, 'Vars', [x1, x2, x3, x4]);
u2_func = matlabFunction(u2_final, 'Vars', [x1, x2, x3, x4]);
fx_func = matlabFunction(fx_sym, 'Vars', [x1, x2, x3, x4]);
gx_func = matlabFunction(gx_sym, 'Vars', [x1, x2, x3, x4]);

% Define closed-loop dynamics as a function for ode45
closed_loop_dyn = @(t, x) closed_loop_dynamics(x, u1_func, u2_func, fx_func, gx_func);

% simulate the closed-loop system from the safe points using ode45 (Runge-Kutta 4/5)
sim_time = 8; % simulation time
tspan = [0, sim_time];
trajectories = cell(size(safe_points, 1), 1);
time_vectors = cell(size(safe_points, 1), 1);

% Set ODE solver options for accuracy
ode_opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'MaxStep', 0.001);

for i = 1:size(safe_points, 1)
    x0 = safe_points(i, :)';
    [t_sol, x_sol] = ode45(closed_loop_dyn, tspan, x0, ode_opts);
    trajectories{i} = x_sol;
    time_vectors{i} = t_sol;
end

% plot the trajectories in the output space (h(x) space)
figure;
hold on;
theta = linspace(0, 2 * pi, 100);
% plot safe set boundary
plot(cos(theta), sin(theta), 'g--', 'LineWidth', 1.5);
% plot target set boundary
plot(0.1 * cos(theta), 0.1 * sin(theta), 'r--', 'LineWidth', 1.5);
% plot trajectories
for i = 1:length(trajectories)
    traj = trajectories{i};
    h1_traj = traj(:, 1);
    h2_traj = traj(:, 3);
    plot(h1_traj, h2_traj, 'b-', 'LineWidth', 0.8);
end

% plot start points
for i = 1:length(trajectories)
    traj = trajectories{i};
    scatter(traj(1, 1), traj(1, 3), 50, 'o', 'filled', 'MarkerFaceColor', [1 0.5 0], 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
end

% plot end points
for i = 1:length(trajectories)
    traj = trajectories{i};
    scatter(traj(end, 1), traj(end, 3), 50, 'x', 'MarkerEdgeColor', [0 0.6 0], 'LineWidth', 1.5);
end

xlabel('h_1 (x_1)');
ylabel('h_2 (x_3)');
title('Trajectories in Output Space under Reach-Avoid Backstepping Controller');
legend('Safe Set Boundary', 'Target Set Boundary', 'Trajectories', 'Start Points', 'End Points');
axis equal;
hold off;

% =====================================================

% =====================================================
% Helper function: convert pvar polynomial to symbolic
% =====================================================
function sym_expr = pvar2sym(pvar_poly, pvar_vars, sym_vars)
    % Convert a pvar polynomial to a symbolic expression
    %
    % INPUTS:
    % pvar_poly: polynomial object from SOSTOOLS
    % pvar_vars: pvar variables used in the polynomial (column vector)
    % sym_vars: symbolic variables to use in the result (column vector)
    %
    % OUTPUTS:
    % sym_expr: symbolic expression

    % Get coefficients and degree matrix
    coeffs = pvar_poly.coefficient;
    degmat = pvar_poly.degmat;

    % Initialize symbolic expression
    sym_expr = sym(0);

    % Build the symbolic expression term by term
    for i = 1:length(coeffs)
        % Start with the coefficient
        term = coeffs(i);

        % Multiply by each variable raised to its power
        for j = 1:length(sym_vars)

            if degmat(i, j) > 0
                term = term * sym_vars(j) ^ degmat(i, j);
            end

        end

        sym_expr = sym_expr + term;
    end

end

% =====================================================
% Helper function: convert symbolic to pvar polynomial
% =====================================================
function pvar_poly = sym2pvar(sym_expr, sym_vars, pvar_vars)
    % Convert a symbolic expression to a pvar polynomial
    %
    % INPUTS:
    % sym_expr: symbolic expression (must be a polynomial)
    % sym_vars: symbolic variables used in the expression (column vector)
    % pvar_vars: pvar variables to use in the result (column vector)
    %
    % OUTPUTS:
    % pvar_poly: pvar polynomial object

    % Expand the expression to ensure it's in polynomial form
    sym_expr = expand(sym_expr);

    % Initialize the pvar polynomial
    pvar_poly = polynomial(0);

    % Get all terms by converting to char and parsing, or use children
    % Get the terms of the polynomial (each additive term)
    terms = children(sym_expr);

    % If children returns empty or the expression itself, it's a single term
    if isempty(terms) || (length(terms) == 1 && isequal(terms{1}, sym_expr))
        terms = {sym_expr};
    end

    n_vars = length(sym_vars);

    for t = 1:length(terms)
        term_sym = terms{t};

        % Extract coefficient and powers for this term
        % Get the coefficient (the numeric part)
        [coeff_sym, mono_sym] = coeffs(term_sym, sym_vars);

        if isempty(coeff_sym)
            % It's a constant
            pvar_poly = pvar_poly + double(term_sym);
        else

            for k = 1:length(coeff_sym)
                coeff_val = double(coeff_sym(k));
                pvar_mono = coeff_val;

                % Get degree of each variable in this monomial
                for j = 1:n_vars
                    deg = polynomialDegree(mono_sym(k), sym_vars(j));

                    if deg > 0
                        pvar_mono = pvar_mono * pvar_vars(j) ^ deg;
                    end

                end

                pvar_poly = pvar_poly + pvar_mono;
            end

        end

    end

end

% =====================================================
% Helper function: closed-loop dynamics for ode45
% =====================================================
function x_dot = closed_loop_dynamics(x, u1_func, u2_func, fx_func, gx_func)
    % Compute closed-loop dynamics for ODE solver
    %
    % INPUTS:
    % x: state vector (4 x 1)
    % u1_func, u2_func: control function handles
    % fx_func, gx_func: dynamics function handles
    %
    % OUTPUTS:
    % x_dot: state derivative (4 x 1)

    % Compute control inputs
    u1 = u1_func(x(1), x(2), x(3), x(4));
    u2 = u2_func(x(1), x(2), x(3), x(4));

    % Compute dynamics
    fx = fx_func(x(1), x(2), x(3), x(4));
    gx = gx_func(x(1), x(2), x(3), x(4));

    % Closed-loop dynamics
    x_dot = fx + gx * [u1; u2];
end
