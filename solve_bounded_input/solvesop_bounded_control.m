% function to solve bounded control for given reach-avoid backstepping controller
function [ux_opt] = solvesop_bounded_control(ux, k1_sym, certificate, x_vars, y_vars, hx, safe_set, target_set, lb, ub)
    [k1_y, k1_lambda, k1_delta] = solve_vanilla_k1_controller(y_vars, safe_set, target_set, 2, 2);

    % substitue the obtained k1 controller into the original controller expression
    % ux_subs = subs(ux, k1_sym, k1_y);

    % substitue the symbolic y variables with  the output mapping hx
    % ux_subs = subs(ux, y_vars, hx);

    % substitute the obtained k1 controller and output mapping into the certificate expression
    certificate_subs = subs(certificate, k1_sym, k1_y);
    certificate_subs = subs(certificate_subs, y_vars, hx);

    % sample 10 random samples from the state space that satisfy the safe set, target set constraint and certificate constraint
    x_samples = get_random_samples(10, x_vars, y_vars, hx, safe_set, target_set, certificate_subs);

    % substitute the obtained k1 controller and output mapping into the original controller expression to get pseudo ux for scenario optimization programming
    ux_pseudo = subs(ux, k1_sym, k1_y);
    ux_pseudo = subs(ux_pseudo, y_vars, hx);

    % select the samples that satisfy the control input bounds for the pseudo ux
    x_samples_valid = [];

    for i = 1:size(x_samples, 2)
        x_sample = x_samples(:, i);
        ux_pseudo_sample = double(subs(ux_pseudo, x_vars, x_sample));

        if all(ux_pseudo_sample >= lb) && all(ux_pseudo_sample <= ub)
            x_samples_valid = [x_samples_valid, x_sample];
        end

    end

    % then we can use the valid samples to solve the scenario optimization programming (SOP) with SOS constraints to find a feasible k1 controller

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
