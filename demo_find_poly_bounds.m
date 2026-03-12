% demo script to find the lower bound and upper bound of a polynomial using SOS programming
clc;
clear;
close all;

pvar x1 x2 x3 x4;
vars = [x1, x2, x3, x4];
ds = 6; % degree of auxiliary polynomials

poly = 63 * x1 * x2 ^ 2 * x3 ...
    - 20 * x1 * x2 * x3 ^ 3 ...
    - 70 * x1 * x2 ...
    - 2 * x1 * x3 * x4 ...
    + 6 * x2 ^ 4 * x3 ^ 3 ...
    - 14 * x2 ^ 4 * x3 ...
    + 21 * x2 ^ 4 ...
    + 6 * x2 ^ 3 * x3 * x4 ...
    - 35 * x2 ^ 2 * x3 ^ 2 * x4 ...
    - 4 * x2 ^ 2 * x3 * x4 ...
    + 4 * x2 * x3 ^ 4 * x4 ...
    + 14 * x2 * x3 * x4 ...
    - 6 * x3 ^ 2 * x4 ^ 2;

% try simpler poly
% poly = x1 ^ 2 + x2 ^ 2 + x3 ^ 2 + x4 ^ 2 -4;
% poly = x1 ^ 3 + x2 ^ 2 + x3 ^ 3 + x4 ^ 2 -4;

% define the safe set (here is unbounded in 4 dimensional space)
safe_set = 1 - x1 ^ 2 - x3 ^ 2; % safe_set: safe_set >=0 inside safe set

% define the target set (here is also unbounded in 4 dimensional space)
% target_set = 2 * ((x3 - 0.1) / 2) ^ 2 ...
%     + 3 * ((x1 + 0.4) / 3) ^ 4 ...
%     + (3 * x1 + 0.3) ^ 2 * (4 * x3 + 0.2) ^ 2 ... .
%     - 0.01;

% try another simpler target set
target_set = x1 ^ 2 + x3 ^ 2 - 0.01; % target_set: target_set <=0 inside target set

% define the lower bound and upper bound variables
dpvar lb ub;

% ------------------ Find the lower bound ------------------
disp('Finding the lower bound...');

% create the SOS program
prog = sosprogram(vars, lb);

% define the auxiliary polynomials
[prog, s1] = sospolyvar(prog, monomials(vars, 0:ds));
[prog, s2] = sospolyvar(prog, monomials(vars, 0:ds));
[prog, s3] = sospolyvar(prog, monomials(vars, 0:ds));
[prog, s4] = sospolyvar(prog, monomials(vars, 0:ds));

% additional constraints to avoid the unboundedness issue
prog = sosineq(prog, s3);
prog = sosineq(prog, s4);

% add the SOS constraints for lower bound
limit_val = 20; % CONSERVATIVE for finding the correct lower bound (10,15,20 works, larger limits output UNKNOWN)
prog = sosineq(prog, poly - lb ...
    - s1 * safe_set ...
    - s2 * target_set ...
    - s3 * (limit_val - x2 ^ 2 - x4 ^ 2) ... % additional constraints to avoid unboundedness
);
% - s3 * (limit_val - x2 ^ 2) ... % additional constraints to avoid unboundedness
% - s4 * (limit_val - x4 ^ 2) ... % additional constraints to avoid unboundedness
% );
% add constrains for auxiliary polynomials to be SOS
prog = sosineq(prog, s1);
prog = sosineq(prog, s2);
prog = sosineq(prog, s3);
prog = sosineq(prog, s4);

% set objective to maximize the lower bound
prog = sossetobj(prog, -lb);
% solve the SOS program
solver_opt.solver = 'mosek';
prog = sossolve(prog, solver_opt);
lb_value = sosgetsol(prog, lb);
disp("lower bound: " + poly2str(lb_value));

% ------------------ Find the upper bound ------------------
disp('Finding the upper bound...');
% create the SOS program
prog = sosprogram(vars, ub);

% define the auxiliary polynomials
[prog, s1] = sospolyvar(prog, monomials(vars, 0:ds));
[prog, s2] = sospolyvar(prog, monomials(vars, 0:ds));
[prog, s3] = sospolyvar(prog, monomials(vars, 0:ds));
[prog, s4] = sospolyvar(prog, monomials(vars, 0:ds));
% additional constraints to avoid the unboundedness issue

prog = sosineq(prog, s3);
prog = sosineq(prog, s4);

% add the SOS constraints for upper bound
prog = sosineq(prog, ub - poly ...
    - s1 * safe_set ...
    - s2 * (-target_set) ...
    - s3 * (limit_val - x2 ^ 2 - x4 ^ 2) ... % additional constraints to avoid unboundedness
);
% - s3 * (limit_val - x2 ^ 2) ... % additional constraints to avoid unboundedness
% - s4 * (limit_val - x4 ^ 2) ... % additional constraints to avoid unboundedness
% );
% add constrains for auxiliary polynomials to be SOS
prog = sosineq(prog, s1);
prog = sosineq(prog, s2);
prog = sosineq(prog, s3);
prog = sosineq(prog, s4);
% set objective to minimize the upper bound
prog = sossetobj(prog, ub);
% solve the SOS program
solver_opt.solver = 'mosek';
% turn off verbose output
solver_opt.verbose = 0;
prog = sossolve(prog, solver_opt);
ub_value = sosgetsol(prog, ub);
disp("upper bound: " + poly2str(ub_value));

% print all results together
disp("Final Results:");
disp("Lower Bound: " + poly2str(lb_value));
disp("Upper Bound: " + poly2str(ub_value));

function s = poly2str(p)
    c = char(p);
    s = c{1};
end
