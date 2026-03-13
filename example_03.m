clear; clc; close all;
pvar y1 y2;

% define the obstacle set 0: a box to reprensent the unsafe region
obs_lower = [1.9; 0];
obs_upper = [2.1; 2];
% use a matrix to represent the box constraints: A * Y <= b
obs_0 = [y1 - obs_lower(1); obs_upper(1) - y1; y2 - obs_lower(2); obs_upper(2) - y2];

% define the second obstacle set 1: a box to reprensent the unsafe region
obs_lower_1 = [3.9; 1];
obs_upper_1 = [4.1; 3];
obs_1 = [y1 - obs_lower_1(1); obs_upper_1(1) - y1; y2 - obs_lower_1(2); obs_upper_1(2) - y2];

% working space bounds
ws_lower = [0; 0];
ws_upper = [6; 3];
% define the working space as box
ws_poly = [y1 - ws_lower(1); ws_upper(1) - y1; y2 - ws_lower(2); ws_upper(2) - y2];

% next define the a set of safe region to ensure the solved polynomial is positive in this region (circle)
safe_region_0 = (y1 - 0) ^ 2 + (y2 -1) ^ 2 - 0.01 ^ 2;
% another safe region (circle)
safe_region_1 = (y1 - 6) ^ 2 + (y2 - 1) ^ 2 - 0.01 ^ 2;

% define the monomials for the SOS program
ds = 7; % degree of the SOS polynomials
monomials_list = monomials([y1; y2], 0:ds);
prog = sosprogram([y1; y2]);

% define auxiliary SOS polynomials for the S-procedure
[prog, s_obs_00] = sospolyvar(prog, monomials_list);
[prog, s_obs_01] = sospolyvar(prog, monomials_list);
[prog, s_obs_02] = sospolyvar(prog, monomials_list);
[prog, s_obs_03] = sospolyvar(prog, monomials_list);
[prog, s_obs_10] = sospolyvar(prog, monomials_list);
[prog, s_obs_11] = sospolyvar(prog, monomials_list);
[prog, s_obs_12] = sospolyvar(prog, monomials_list);
[prog, s_obs_13] = sospolyvar(prog, monomials_list);
[prog, s_safe_0] = sospolyvar(prog, monomials_list);
[prog, s_safe_1] = sospolyvar(prog, monomials_list);
[prog, s_ws_0] = sospolyvar(prog, monomials_list);
[prog, s_ws_1] = sospolyvar(prog, monomials_list);
[prog, s_ws_2] = sospolyvar(prog, monomials_list);
[prog, s_ws_3] = sospolyvar(prog, monomials_list);

% define the decision variable for the polynomial we want to solve for
[prog, p] = sospolyvar(prog, monomials_list);

% construct the SOS constraints using the S-procedure
% p <= 0 inside obstacle 0: -p - s_obs_0' * obs_0 >= 0 (SOS everywhere)
prog = sosineq(prog, -p - [s_obs_00; s_obs_01; s_obs_02; s_obs_03]' * obs_0);
% p <= 0 inside obstacle 1: -p - s_obs_1' * obs_1 >= 0 (SOS everywhere)
prog = sosineq(prog, -p - [s_obs_10; s_obs_11; s_obs_12; s_obs_13]' * obs_1);
% p >= 1 inside safe circle 0: p + s_safe_0 * safe_region_0 - 1 >= 0 (SOS everywhere)
% (safe_region_0 <= 0 inside circle, so -s_safe_0*safe_region_0 >= 0, giving p >= 1 inside)
prog = sosineq(prog, p + s_safe_0 * safe_region_0);
% p >= 1 inside safe circle 1
prog = sosineq(prog, p + s_safe_1 * safe_region_1);

% negative outside the working space: -p - s_ws' * ws_poly >= 0 (SOS everywhere)
prog = sosineq(prog, -p - [s_ws_0; s_ws_1; s_ws_2; s_ws_3]' * ws_poly);

prog = sosineq(prog, s_obs_00);
prog = sosineq(prog, s_obs_01);
prog = sosineq(prog, s_obs_02);
prog = sosineq(prog, s_obs_03);
prog = sosineq(prog, s_obs_10);
prog = sosineq(prog, s_obs_11);
prog = sosineq(prog, s_obs_12);
prog = sosineq(prog, s_obs_13);
prog = sosineq(prog, s_safe_0);
prog = sosineq(prog, s_safe_1);
prog = sosineq(prog, s_ws_0);
prog = sosineq(prog, s_ws_1);
prog = sosineq(prog, s_ws_2);
prog = sosineq(prog, s_ws_3);

waypoints = [2.0, 2.5;
             2.5, 2;
             3.0, 1.5;
             3.5, 1;
             4.0, 0.5];

for i = 1:size(waypoints, 1)
    prog = sosineq(prog, subs(p, [y1; y2], [waypoints(i, 1); waypoints(i, 2)]) - 0.5);
end

% solve the SOS program
solver_opt.solver = 'mosek';

solver_opt.verbose = 0;

prog = sossolve(prog, solver_opt);
% extract the solution
p_sol = sosgetsol(prog, p);

% evaluate p_sol on a grid using subs
y1_range = linspace(-1, 7, 200);
y2_range = linspace(-1, 4, 200);
[Y1g, Y2g] = meshgrid(y1_range, y2_range);
Y1_flat = Y1g(:)';
Y2_flat = Y2g(:)';
p_flat = double(subs(p_sol, [y1; y2], [Y1_flat; Y2_flat]));
p_grid = reshape(p_flat, size(Y1g));
fprintf('p_sol range in domain: [%.4f, %.4f]\n', min(p_flat), max(p_flat));

% visualize the results
figure;
% two-tone fill: gray for p < 0, blue for p >= 0
contourf(Y1g, Y2g, p_grid, [0, 1e6], 'LineColor', 'none');
colormap([[0.85 0.85 0.85]; [0.5 0.7 1.0]]);
hold on;
contour(Y1g, Y2g, p_grid, [0 0], 'b-', 'LineWidth', 2); % zero level set boundary
% plot the obstacle sets
rectangle('Position', [obs_lower(1), obs_lower(2), obs_upper(1) - obs_lower(1), obs_upper(2) - obs_lower(2)], 'EdgeColor', 'r', 'LineWidth', 2);
rectangle('Position', [obs_lower_1(1), obs_lower_1(2), obs_upper_1(1) - obs_lower_1(1), obs_upper_1(2) - obs_lower_1(2)], 'EdgeColor', 'r', 'LineWidth', 2);
% plot the safe regions (zero level set of safe_region polynomials)
sr0_grid = reshape(double(subs(safe_region_0, [y1; y2], [Y1_flat; Y2_flat])), size(Y1g));
sr1_grid = reshape(double(subs(safe_region_1, [y1; y2], [Y1_flat; Y2_flat])), size(Y1g));
contour(Y1g, Y2g, sr0_grid, [0 0], 'g-', 'LineWidth', 2);
contour(Y1g, Y2g, sr1_grid, [0 0], 'g-', 'LineWidth', 2);
% plot the contour of the solved polynomial
xlabel('y1');
ylabel('y2');
axis equal;
title('SOS Polynomial with S-Procedure Constraints');
