% helper function to solve the k1 controller for the single-integrator system considering control input bounds
function [k1_opt, k1_lambda, k1_delta] = solve_k1_controller_sop(ux, k1_sym, certificate, x_vars, y_vars, hx, safe_set, target_set, x_samples_valid, lb, ub, ds, dv)
    % solve the k1 controller for the single-integrator system considering control input bounds using scenario optimization programming (SOP) with SOS constraints
    % the single-integrator system is defined as y_dot = k1, where y is the output of the original system
    % here we aims to find a k1 controller such that the single-integrator system can satisfy the reach-avoid specification defined by the safe set and target set while ensuring the control input bounds are satisfied

    % INPUTS:
    % ux: original controller expression obtained from backstepping design, symbolic vector of size (m x 1), m is the control input dimension
    % k1_sym: symbolic vector for the k1 controller in the original controller expression, symbolic vector of size (p x 1), p is the output dimension
    % certificate: symbolic expression for the certificate obtained from backstepping design, should be non-negative for the valid state space
    % x_vars: state variables of the original system, symbolic vector of size (n x 1)
    % y_vars: output variables of the original system, symbolic vector of size (p x 1)
    % hx: output mapping from state variables to output variables, symbolic vector of size (p x 1)
    % safe_set: symbolic expression for the safe set constraint, should be non-negative inside the safe set
    % target_set: symbolic expression for the target set constraint, should be non-positive out of the target set
    % x_samples_valid: random samples for state variables that satisfy the safe set, target set constraint and certificate constraint, matrix of size (n x num_valid_samples)
    % lb: lower bounds for control inputs, vector of size (m x 1)
    % ub: upper bounds for control inputs, vector of size (m x 1)
    % ds: degree of the auxiliary SOS polynomials for the single-integrator system
    % dv: degree of the k1 controller polynomial

    % create pvar variables  for the x_vars
    x_vars_pvar = [];

    for i = 1:length(x_vars)
        x_vars_pvar = [x_vars_pvar; pvar(strcat(char(x_vars(i)), '_pvar'))];
    end

    % get the number of outputs
    p = length(y_vars);
    % create pvar variables for the single-intgegrator system
    y_vars_pvar = [];

    for i = 1:p
        y_vars_pvar = [y_vars_pvar; pvar(strcat(char(y_vars(i)), '_pvar'))];
    end

    % define decision variable for solving k1 controller
    dpvar delta lambda;

    xi0 = 1e-8; % small positive value to ensure the SOS constraints are strictly feasible

    % define monomials for k1 controller and auxiliary SOS polynomials
    monom_k1 = monomials(y_vars_pvar, 0:dv);
    monom_sos = monomials(y_vars_pvar, 0:ds);

    % define the SOS program for solving k1 controller
    prog = sosprogram(y_vars_pvar, [delta; lambda]);

    % create polynomials for every single element in k1 controller
    k1_poly = [];

    for i = 1:p
        [prog, k1_i] = sospolyvar(prog, monom_k1);
        k1_poly = [k1_poly; k1_i];
    end

    % create pvar variables for the k1 controller
    k1_pvar = [];

    for i = 1:p
        k1_pvar = [k1_pvar; pvar(strcat(char(k1_sym(i)), '_pvar'))];
    end

    % create the pvar polynoimials for the system output mapping hx(x)
    hx_pvar = polynomial(zeros(length(hx), 1));

    for i = 1:length(hx)
        hx_pvar(i) = sym2pvar(hx(i), x_vars, x_vars_pvar);
    end

    % create the pvar version of the safe set and target set by substituting y_vars with y_vars_pvar
    safe_set_pvar = sym2pvar(safe_set, y_vars, y_vars_pvar);
    target_set_pvar = sym2pvar(target_set, y_vars, y_vars_pvar);

    % define auxiliary SOS polynomials
    [prog, s1] = sospolyvar(prog, monom_sos);
    [prog, s2] = sospolyvar(prog, monom_sos);

    % before adding constraints, compute the Lie derivative of the safe_set constraint along the single-integrator system: L_f(safe_set) = ∂safe_set/∂y * k1

    % Compute gradient of safe_set as a pvar polynomial
    grad_safe_set = [];

    for i = 1:p
        grad_safe_set = [grad_safe_set; diff(safe_set_pvar, y_vars_pvar(i))];
    end

    % Compute Lie derivative as dot product: ∇S · k1
    % This keeps it as a proper dpvar polynomial
    Lie_safe_set = grad_safe_set(1) * k1_poly(1);

    for i = 2:p
        Lie_safe_set = Lie_safe_set + grad_safe_set(i) * k1_poly(i);
    end

    % disp the Lie derivative to check the correctness of the computation
    % disp('Lie derivative of the safe set constraint along the single-integrator system:');
    % disp(Lie_safe_set);

    % Construct the SOS constraint expression
    sos_expr = Lie_safe_set - lambda * safe_set_pvar + delta - s1 * safe_set_pvar - s2 * target_set_pvar;

    % disp('The expression for the SOS constraint for k1 controller:');
    % disp(sos_expr);

    % add the SOS constraints for k1 controller
    prog = sosineq(prog, sos_expr);
    prog = sosineq(prog, delta);
    prog = sosineq(prog, lambda - xi0); % ensure lambda is positive
    prog = sosineq(prog, s1);
    prog = sosineq(prog, s2);

    % add additional constraints to ensure the control input bounds are satisfied for the sampled states

    for i = 1:length(ux)
        % get the numerator and denominator of the symbolic expression ux(i)
        [num, den] = numden(ux(i));

        % for all samples
        for j = 1:size(x_samples_valid, 2)
            x_sample = x_samples_valid(:, j);
            den_val = double(subs(den, x_vars, x_sample)); % evaluate the denominator for the sampled state

            if abs(den_val) < 1e-6 % if the denominator is close to zero, we skip this sample to avoid numerical issues
                continue;
            end

            % substitute the k1_sym with k1_pvar to make it a proper dpvar polynomial that can be evaluated for the sampled states and used in the SOS constraints
            num_pvar = sym2pvar(num, [x_vars; y_vars; k1_sym], [x_vars_pvar; y_vars_pvar; k1_pvar]);
            num_pvar = subs(num_pvar, k1_pvar, dpvar2poly(k1_poly));
            % substitute the y_vars_pvar variables with the output mapping hx
            num_pvar = subs(num_pvar, y_vars_pvar, hx_pvar);
            % now the expression only has x_vars_pvar variables
            num_pvar = subs(num_pvar, x_vars_pvar, x_sample); % evaluate the numerator for the sampled state

            % then add the constraint for the control input bounds: lb(i) <= num_pvar/den_val <= ub(i)
            prog = sosineq(prog, num_pvar / den_val - lb(i)); % ensure num_pvar/den_val >= lb(i) -> num_pvar/den_val - lb(i) >= 0
            prog = sosineq(prog, ub(i) - num_pvar / den_val); % ensure num_pvar/den_val <= ub(i) -> ub(i) - num_pvar/den_val >= 0

        end

    end

    disp('Added SOS constraints for k1 controller with bounded control inputs for the valid samples.');

    % set objective to minimize delta
    prog = sossetobj(prog, delta);

    % set solver options and solve the SOS program
    solver_opt.solver = 'mosek';
    prog = sossolve(prog, solver_opt);

    % extract the obtained k1 controller

    k1_opt = [];

    for i = 1:p
        k1_opt_i = sosgetsol(prog, k1_poly(i));
        k1_opt = [k1_opt; poly2sym(k1_opt_i, y_vars_pvar, y_vars)]; % vector of symbolic polynomials for k1 controller
    end

    % extract the obtained lambda and delta
    k1_lambda = double(sosgetsol(prog, lambda)); % extract lambda as a double value
    k1_delta = double(sosgetsol(prog, delta)); % extract delta as a double value
end
