% function to solve bounded control for given reach-avoid backstepping controller
function [ux_opt, certificate_opt] = solvesop_bounded_control(ux, k1_sym, mu, lambda, certificate, x_vars, y_vars, hx, safe_set, target_set, mu_val, lb, ub, ds, dv, samples_num)

    % INPUTS:
    % ux: original controller expression obtained from backstepping design, symbolic vector of size (m x 1), m is the control input dimension
    % k1_sym: symbolic vector for the k1 controller in the original controller expression, symbolic vector of size (p x 1), p is the output dimension
    % certificate: symbolic expression for the certificate obtained from backstepping design, should be non-negative for the valid state space
    % x_vars: state variables of the original system, symbolic vector of size (n x 1)
    % y_vars: output variables of the original system, symbolic vector of size (p x 1)
    % hx: output mapping from state variables to output variables, symbolic vector of size (p x 1)
    % safe_set: symbolic expression for the safe set constraint, should be non-negative inside the safe set
    % target_set: symbolic expression for the target set constraint, should be non-positive out of the target set
    % Ax: matrix for linear constraints on control inputs, matrix of size (num_constraints x m)
    % lb: lower bounds for control inputs, vector of size (m x 1)
    % ub: upper bounds for control inputs, vector of size (m x 1)
    % ds: degree of the auxiliary SOS polynomials for the single-integrator system
    % dv: degree of the k1 controller polynomial

    [k1_y, k1_lambda, k1_delta] = solve_vanilla_k1_controller(y_vars, safe_set, target_set, dv, ds);

    % substitute the obtained k1 controller and output mapping into the certificate expression
    certificate_subs = subs(certificate, k1_sym, k1_y);
    certificate_subs = subs(certificate_subs, y_vars, hx);

    certificate_subs = substitute_mu_lambda(certificate_subs, mu, lambda, mu_val, k1_lambda); % use the computed value from the vanilla k1 controller design as the lambda value for solving the bounded control inputs

    % sample some random samples from the state space that satisfy the safe set, target set constraint and certificate constraint
    x_samples = get_random_samples(samples_num, x_vars, y_vars, hx, safe_set, target_set, certificate_subs);

    % substitute the obtained k1 controller and output mapping into the original controller expression to get pseudo ux for scenario optimization programming
    ux_pseudo = subs(ux, k1_sym, k1_y);
    ux_pseudo = subs(ux_pseudo, y_vars, hx);

    % also substitute the obtained mu_val and lambda_val into the original controller expression

    ux_pseudo = substitute_mu_lambda(ux_pseudo, mu, lambda, mu_val, k1_lambda); % use the computed value from the vanilla k1 controller design as the lambda value for solving the bounded control inputs

    % select the samples that satisfy the control input bounds for the pseudo ux
    % Preallocate assuming all samples might be valid (worst case)
    x_samples_valid = zeros(size(x_samples, 1), size(x_samples, 2));
    valid_count = 0;

    for i = 1:size(x_samples, 2)
        x_sample = x_samples(:, i);
        ux_pseudo_sample = double(subs(ux_pseudo, x_vars, x_sample));

        if all(ux_pseudo_sample >= lb) && all(ux_pseudo_sample <= ub)
            valid_count = valid_count + 1;
            x_samples_valid(:, valid_count) = x_sample;
        end

    end

    % Trim to actual number of valid samples
    x_samples_valid = x_samples_valid(:, 1:valid_count);

    % then we can use the valid samples to solve the scenario optimization programming (SOP) with SOS constraints to find a feasible k1 controller
    [k1_opt, k1_lambda_opt, k1_delta_opt] = solve_k1_controller_sop(ux_pseudo, k1_sym, certificate, x_vars, y_vars, hx, safe_set, target_set, x_samples_valid, lb, ub, ds, dv);

    % substitute the obtainted k1 controller into the original controller expression ux to get the final controller with bounded control inputs
    ux_opt = subs(ux, k1_sym, k1_opt);
    ux_opt = subs(ux_opt, y_vars, hx);
    % substitute the obtained k1_lambda_opt and k1_delta_opt into the final controller expression ux_opt
    ux_opt = substitute_mu_lambda(ux_opt, mu, lambda, mu_val, k1_lambda_opt); % use the computed value from the SOP design as the lambda value for the final controller expression

    % substitute the obtainted k1_opt controller into the original certificate expression to get the final certificate with bounded control inputs
    certificate_opt = subs(certificate, k1_sym, k1_opt);
    certificate_opt = subs(certificate_opt, y_vars, hx);

    certificate_opt = substitute_mu_lambda(certificate_opt, mu, lambda, mu_val, k1_lambda_opt); % use the computed value from the SOP design as the lambda value for the final certificate expression

end

% helper functikon to substitute mu values and lambda with predefined numerical values in the controller expression
function u_subs = substitute_mu_lambda(expr, mu, lambda, mu_val, lambda_val)
    expr_subs = expr;

    % substitute mu values
    for i = 1:numel(mu)

        if ~isempty(mu{i})
            expr_subs = subs(expr_subs, mu{i}, mu_val);
        end

    end

    % substitute lambda value
    u_subs = subs(expr_subs, lambda, lambda_val);

end

% helper function to get random samples from the state space that satisfy the safe set, target set constraint and certificate constraint
function x_samples = get_random_samples(num_samples, x_vars, y_vars, hx, safe_set, target_set, certificate)
    % num_samples: number of random samples to generate
    % x_vars: state variables, symbolic vector of size (n x 1)
    % y_vars: output variables, symbolic vector of size (p x 1)
    % safe_set: symbolic expression for the safe set constraint, should be non-negative inside the safe set
    % target_set: symbolic expression for the target set constraint, should be non-positive inside the target set
    % certificate: symbolic expression for the certificate, should be non-negative for the valid state space
    % x_samples: generated random samples for state variables, matrix of size (n x num_samples)
    % y_samples: generated random samples for output variables, matrix of size (p x num_samples)

    n = length(x_vars);
    p = length(y_vars);
    x_samples = [];
    y_samples = [];

    while size(x_samples, 2) < num_samples
        % generate random sample for state variables
        x_sample = randn(n, 1); % you can adjust the distribution and range of random samples as needed
        % evaluate the output variables using the output mapping hx (assuming hx is a function of x_vars)
        y_sample = double(subs(hx, x_vars, x_sample));

        % check if the sample satisfies the safe set constraint, target set constraint and certificate constraint
        if double(subs(safe_set, [x_vars; y_vars], [x_sample; y_sample])) >= 0 && ... % inside the safe set
                double(subs(target_set, [x_vars; y_vars], [x_sample; y_sample])) >= 0 && ... % outside the target set
                double(subs(certificate, [x_vars; y_vars], [x_sample; y_sample])) >= 0 % satisfies the certificate constraint
            x_samples = [x_samples, x_sample];
            y_samples = [y_samples, y_sample];
        end

    end

end
