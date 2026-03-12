% Test script for computing polynomial bounds
clc;
clear;
close all;

syms x1 x2;
num_poly = 0.01 - x1 ^ 2 - x2 ^ 2;
den_poly = x1 ^ 2 + x2 ^ 2 - 0.04;
superlevel_set_poly = 1 - x1 ^ 2 - x2 ^ 2;
ds_min = 2;

% Safety margin to exclude the singular ring where the denominator is zero
epsilon = 1e-9;

[lower_bound, upper_bound] = compute_poly_bounds_safe(num_poly, den_poly, superlevel_set_poly, ds_min, epsilon);

disp('=== SOS Bounds ===');
disp(['Computed lower bound: ', num2str(lower_bound)]);
disp(['Computed upper bound: ', num2str(upper_bound)]);

% --- Sampling validation ---
n_samples = 50000;
bounds = [-1, 1; -1, 1];
x1_samples = bounds(1, 1) + (bounds(1, 2) - bounds(1, 1)) * rand(n_samples, 1);
x2_samples = bounds(2, 1) + (bounds(2, 2) - bounds(2, 1)) * rand(n_samples, 1);

superlevel_set_func = matlabFunction(superlevel_set_poly, 'Vars', [x1, x2]);
num_poly_func = matlabFunction(num_poly, 'Vars', [x1, x2]);
den_poly_func = matlabFunction(den_poly, 'Vars', [x1, x2]);

% Important: the epsilon region must also be excluded during sampling for a fair comparison.
feasible_indices = find(superlevel_set_func(x1_samples, x2_samples) >= 0 ...
    & abs(den_poly_func(x1_samples, x2_samples)) >= epsilon);

rational_poly_values = num_poly_func(x1_samples(feasible_indices), x2_samples(feasible_indices)) ...
    ./ den_poly_func(x1_samples(feasible_indices), x2_samples(feasible_indices));

disp('=== Sampling Bounds (Excluding epsilon region) ===');
disp(['Min value at feasible samples: ', num2str(min(rational_poly_values))]);
disp(['Max value at feasible samples: ', num2str(max(rational_poly_values))]);

% =========================================================================
% Main Wrapper Function: compute lower and upper bounds of num_poly/den_poly
% over {superlevel_set_poly >= 0}, excluding the epsilon-singular region.
% =========================================================================
%
function [lower_bound, upper_bound] = compute_poly_bounds_safe(num_poly, den_poly, superlevel_set_poly, ds_min, epsilon)
    % Convert sym inputs to pvar
    if isa(num_poly, 'sym')
        sym_vars = symvar(num_poly + den_poly + superlevel_set_poly);
        var_names = arrayfun(@char, sym_vars, 'UniformOutput', false);
        eval(['pvar ' strjoin(var_names, ' ')]);
        pvar_vars = eval(['[' strjoin(var_names, ',') ']']).';
        num_poly = sym2pvar(num_poly, sym_vars, pvar_vars);
        den_poly = sym2pvar(den_poly, sym_vars, pvar_vars);
        superlevel_set_poly = sym2pvar(superlevel_set_poly, sym_vars, pvar_vars);
    end

    var_names = unique([num_poly.varname; den_poly.varname; superlevel_set_poly.varname]);
    eval(['pvar ' strjoin(var_names, ' ')]);
    vars = eval(['[' strjoin(var_names, ', ') ']']);

    % Four explicit calls; algebraic sign logic is encapsulated in the helper.
    lb_pos = solve_fraction_bound(vars, num_poly, den_poly, superlevel_set_poly, ds_min, epsilon, 1, 'lower');
    ub_pos = solve_fraction_bound(vars, num_poly, den_poly, superlevel_set_poly, ds_min, epsilon, 1, 'upper');

    lb_neg = solve_fraction_bound(vars, num_poly, den_poly, superlevel_set_poly, ds_min, epsilon, -1, 'lower');
    ub_neg = solve_fraction_bound(vars, num_poly, den_poly, superlevel_set_poly, ds_min, epsilon, -1, 'upper');

    lower_bound = min(lb_pos, lb_neg);
    upper_bound = max(ub_pos, ub_neg);
end

% =========================================================================
% Helper: solve one rational bound (lower or upper) for a fixed sign of Q
% sign_Q =  1 : denominator region  {den_poly >= epsilon}
% sign_Q = -1 : denominator region  {den_poly <= -epsilon}
% =========================================================================
function bound_val = solve_fraction_bound(poly_vars, num_poly, den_poly, S_poly, ds_min, epsilon, sign_Q, bound_type)
    ds = ds_min + mod(ds_min, 2);
    prog = sosprogram(poly_vars);

    dpvar b;
    prog = sosdecvar(prog, b);

    [prog, s_S] = sospolyvar(prog, monomials(poly_vars, 0:ds));
    [prog, s_Q] = sospolyvar(prog, monomials(poly_vars, 0:ds));

    % Archimedean ball: required for compactness so the SOS problem is well-posed
    % (radius 10 is sufficient to cover the unit disk)
    R2 = 10;
    ball_poly = R2 - sum(poly_vars .^ 2);
    [prog, s_ball] = sospolyvar(prog, monomials(poly_vars, 0:ds));

    % Core algebraic case analysis based on the sign of the denominator:
    if sign_Q == 1 % Region: den_poly >= epsilon  (denominator positive)
        Q_expr = den_poly - epsilon;

        if strcmp(bound_type, 'lower')
            % P/Q >= b  =>  P - b*Q >= 0  (Q positive, inequality preserved)
            main_eq = num_poly - b * den_poly;
            prog = sossetobj(prog, -b); % maximize lower bound
        else
            % P/Q <= b  =>  b*Q - P >= 0
            main_eq = b * den_poly - num_poly;
            prog = sossetobj(prog, b); % minimize upper bound
        end

    else % Region: den_poly <= -epsilon  (i.e. -den_poly >= epsilon, denominator negative)
        Q_expr = -den_poly - epsilon;

        if strcmp(bound_type, 'lower')
            % P/Q >= b  =>  b*Q - P >= 0  (Q negative, inequality flips!)
            main_eq = b * den_poly - num_poly;
            prog = sossetobj(prog, -b); % maximize lower bound
        else
            % P/Q <= b  =>  P - b*Q >= 0
            main_eq = num_poly - b * den_poly;
            prog = sossetobj(prog, b); % minimize upper bound
        end

    end

    % Assemble the S-procedure constraint
    prog = sosineq(prog, main_eq - s_S * S_poly - s_Q * Q_expr - s_ball * ball_poly);

    prog = sosineq(prog, s_S);
    prog = sosineq(prog, s_Q);
    prog = sosineq(prog, s_ball);

    solver_opt.solver = 'mosek';
    solver_opt.verbose = 0;
    prog = sossolve(prog, solver_opt);

    b_sol = sosgetsol(prog, b);

    if isempty(b_sol)
        bound_val = NaN;
    else
        bound_val = double(b_sol);
    end

end
