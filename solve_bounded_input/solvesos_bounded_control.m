% function solve bounded control with given reach-avoid backstepping controller
function [ux_opt] = solvesos_bounded_control(ux, x_vars, y_vars, k1, hx, safe_set, target_set, bound_Ax, bound_lbx, bound_ubx)

    % Description:
    % This function solves for the optimal control input ux_opt that satisfies the bounded control constraints
    % using the provided reach-avoid backstepping controller ux_reach_avoid.
    % here support the bound defined as
    %  bound_lbx <= A(x) * u(x) <= bound_ubx
    % where A(x) is a symbolic polynomial matrix, and bound_lbx, bound_ubx are symbolic polynomial vectors.
    % u(x) = [u_1(x), ... , u_m(x)] here is rational polynomial vector defined as ux_i = num_ux_i / den_ux_i for each element.

    % INPUTS:
    % ux: rational polynomial vector representing the control input, symbolic polynomial vector of size (m x 1)
    % x_vars: symbolic variables representing the state variables, symbolic vector of size (n x 1)
    % hx: symbolic polynomial vector representing the system outputs, symbolic polynomial vector of size (p x 1)
    % bound_Ax: symbolic polynomial matrix A(x) defining the control bounds of size (p x m), also has no decision variables
    % bound_lbx: symbolic polynomial vector representing the lower bounds of size (p x 1), also has no decision variables
    % bound_ubx: symbolic polynomial vector representing the upper bounds of size (p x 1), also has no decision variables
    % safe_set: symbolic polynomial representing the safe set constraints
    % target_set: symbolic polynomial representing the target set constraints

    % OUTPUTS:
    % ux_opt: optimal control input satisfying the given bounded control constraints, also symbolic rational polynomial vector of size (m x 1)

    % extract the numerator and denominator of the rational polynomial control input for each element
    m = length(ux);
    % initialize numerator and denominator cell arrays
    num_ux = cell(m, 1);
    den_ux = cell(m, 1);

    for i = 1:m
        [num_ux{i}, den_ux{i}] = numden(ux(i));
    end

    % >>>>>>>>>>>>>>>>>>>>>>>>>> SOS program <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % define pvar variables for the SOS program based on the number of state variables

    % create pvar variables for the state variables
    n = length(x_vars);
    x_vars_pvar = [];

    for i = 1:n
        % create pvar variable by adding suffix '_pvar' to the symbolic variable name
        pvar_var_name = strcat(char(x_vars(i)), '_pvar');
        x_vars_pvar = [x_vars_pvar; pvar(pvar_var_name)];
    end

    % create pvar variables for the output variables
    p = length(y_vars);
    y_vars_pvar = [];

    for i = 1:p
        % create pvar variable by adding suffix '_pvar' to the symbolic variable name
        pvar_var_name = strcat(char(y_vars(i)), '_pvar');
        y_vars_pvar = [y_vars_pvar; pvar(pvar_var_name)];
    end

    % create pvar variables for k1 controller design
    k1_pvar = [];

    for i = 1:length(k1)
        pvar_var_name = strcat(char(k1(i)), '_pvar');
        k1_pvar = [k1_pvar; pvar(pvar_var_name)];
    end

    % declare dpvar for lambda and delta for SOS optimization
    dpvar lambda_dpvar delta_dpvar;

    % settings for the SOS program
    xi = 1e-8; % small positive constant for numerical stability
    deg_sos = 4; % degree for the SOS multipliers
    deg_k1 = 4; % degree for the k1 polynomial

    % create the SOS program for x variables
    prog_x = sosprogram(x_vars_pvar, [lambda_dpvar; delta_dpvar]);
    % create the SOS program for y variables
    prog_y = sosprogram(y_vars_pvar);

    % get the pvar polynomial for the system output mapping hx(x)
    hx_pvar = polynomial(zeros(length(hx), 1));

    for i = 1:length(hx)
        hx_pvar(i) = sym2pvar(hx(i), x_vars, x_vars_pvar);
    end

    % transform the safe set and target set from symbolic polynomial to pvar polynomial for the SOS program
    safe_set_y_pvar = sym2pvar(safe_set, y_vars, y_vars_pvar);
    target_set_y_pvar = sym2pvar(target_set, y_vars, y_vars_pvar);
    % substitute the output mapping hx(x) into the safe set and target set to transform them from y variables to x variables for the prog_x SOS program
    safe_set_x_pvar = subs(safe_set_y_pvar, y_vars_pvar, hx_pvar);
    target_set_x_pvar = subs(target_set_y_pvar, y_vars_pvar, hx_pvar);

    % transform the numeriator and denominator of the control input from symbolic polynomial to pvar polynomial for the SOS program
    num_ux_k1_pvar = cell(m, 1);
    den_ux_x_pvar = cell(m, 1);

    for i = 1:m
        % as for the polynomials here, k1 is still unknown, so we have to also substitute k1 with k1_pvar in the numerator and denominator of the control input
        num_ux_k1_pvar{i} = sym2pvar(num_ux{i}, [x_vars; k1], [x_vars_pvar; k1_pvar]);
        den_ux_x_pvar{i} = sym2pvar(den_ux{i}, [x_vars; k1], [x_vars_pvar; k1_pvar]);
    end

    % define monomials for the SOS multipliers for the prog_x SOS program
    monos_sos_x = monomials(x_vars_pvar, 0:deg_sos);

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>> SOS program <<<<<<<<<<<<<<<<<<<<<<<<<<<
    % add constraints to the SOS program as we did for the single-integrator system, but now we need to transform it to constraints of x

    % ++++++++++++++ Lie derivative constraint +++++++++++++++
    % compute the partial derivative of the safe set (function of y) with respect to system output y, which is the state for the single-integrator system in the backstepping design
    Dy_safe_set = jacobian(safe_set, y_vars);
    % then transform this partial derivative (vector of polynomials) from symbolic polynomial to pvar polynomial
    for i = 1:length(Dy_safe_set)
        Dy_safe_set_pvar(i) = sym2pvar(Dy_safe_set(i), y_vars, y_vars_pvar);
    end

    % define the k1(y) polynomial as a function of y_pvar variables for the prog_y SOS program
    monos_k1_y = monomials(y_vars_pvar, 0:deg_k1);

    % define all k1_i(y) polynomials for each control input as a function of y_pvar variables for the prog_y SOS program
    k1_y = [];

    for i = 1:p
        [prog_y, k1_i_y] = sospolyvar(prog_y, monos_k1_y);
        k1_y = [k1_y; k1_i_y];
    end

    % define k1_y_poly for the subsequent substitution
    k1_y_poly = [];

    for i = 1:length(k1_y)
        k1_y_poly = [k1_y_poly; dpvar2poly(k1_y(i))];
    end

    % get the lie derivative of the safe set along the single-integrator system dynamics defined by \dot{y} = k1(y)
    Ly_safe_set = Dy_safe_set_pvar * k1_y; % this is the Lie derivative of the safe set along the single-integrator system dynamics, which can be used for the prog_y SOS program

    % then substitute the Lie derivative with the output mapping hx(x) to transform the Lie derivative constraint from y variables to x variables for the prog_x SOS program
    Lx_safe_set = subs(Ly_safe_set, y_vars_pvar, hx_pvar); % now Lxy_safe_set is the Lie derivative of the safe set along the single-integrator system dynamics after substitution, which can be used for the prog_x SOS program

    % add the Lie derivative constraint to the prog_x SOS program
    % defien the SOS multiplier for the Lie derivative constraint
    [prog_x, sos_1] = sospolyvar(prog_x, monos_sos_x);
    [prog_x, sos_2] = sospolyvar(prog_x, monos_sos_x);

    % define additional variable bound for numerical feasibility of the SOS program
    feasible_space = 1e3;

    for i = 1:length(x_vars_pvar)
        feasible_space = feasible_space - x_vars_pvar(i) ^ 2;
    end

    [prog_x, sos_3] = sospolyvar(prog_x, monos_sos_x); % additional sos multiplier for bound all variables in the constraint for better numerical stability of the SOS solver

    % add the constraint: Lx_safe_set -lambda_dpvar*safe_set_x_pvar + delta_dpvar -sos_1*safe_set_x_pvar - sos_2*target_set_x_pvar is SOS
    prog_x = sosineq(prog_x, Lx_safe_set - lambda_dpvar * safe_set_x_pvar + delta_dpvar - sos_1 * safe_set_x_pvar - sos_2 * target_set_x_pvar - sos_3 * feasible_space);

    % ++++++++++++++ auxiliary constraints +++++++++++++++
    % add auxiliary constraints to ensure the non-negativity of lambda and delta
    prog_x = sosineq(prog_x, lambda_dpvar - xi); % lambda >= xi > 0
    prog_x = sosineq(prog_x, delta_dpvar); % delta >= 0
    prog_x = sosineq(prog_x, sos_1);
    prog_x = sosineq(prog_x, sos_2);

    % ++++++++++++++ Bounded control constraints +++++++++++++++
    % for each bound constraint defined in bound_Ax * ux <= bound_ubx and bound_Ax * ux >= bound_lbx

    % first we need to transform the obtained nominiators of the control input from y variables to x variables for the prog_x SOS program by substituting the output mapping hx(x) into the nominator of the control input
    num_ux_x_pvar = cell(m, 1);

    for i = 1:m
        % first substitute the k1_pvar with k1_y polynomial, then substitute the y_par with hx_pvar to make it works for the prog_x SOS program
        num_ux_y_pvar = subs(num_ux_k1_pvar{i}, k1_pvar, k1_y_poly);
        num_ux_x_pvar{i} = subs(num_ux_y_pvar, y_vars_pvar, hx_pvar);
        num_ux_x_pvar{i} = poly2dpvar(num_ux_x_pvar{i}); % convert the numerator of the control input to dpvar for the SOS program
    end

    % extract nonzero elements in bound_Ax symbolic polynomial matrix, and convert them to pvar polynomial matrix for the SOS program
    [row_idx, col_idx, val] = find(bound_Ax);

    % group the nonzero elements for each row of bound_Ax together for the subsequent SOS constraint construction
    bound_Ax_row = cell(size(bound_Ax, 1), 1);

    for i = 1:length(val)
        r = row_idx(i);
        c = col_idx(i);
        v = val(i);
        % if v is a scalr, then no need to convert it to pvar polynomial, just store the value; if v is a symbolic polynomial, then we need to convert it to pvar polynomial for the SOS program
        if isnumeric(v)
            bound_Ax_row{r} = [bound_Ax_row{r}; struct('col', c, 'posval', v, "negval", -1 * v)]; % store the nonzero element in the corresponding row group
        else
            % else, convert v to pvar polynomial
            v = sym2pvar(v, x_vars, x_vars_pvar);
            bound_Ax_row{r} = [bound_Ax_row{r}; struct('col', c, 'posval', v, "negval", -1 * v)]; % store the nonzero element in the corresponding row group

        end

    end

    % then we are ready to add the bounded control constraints to the prog_x SOS program
    % for each row of the bound_Ax matrix
    for i = 1:size(bound_Ax, 1)
        % for the lower bound constraint: bound_Ax(i, :) * ux >= bound_lbx(i)
        prog_x = add_bounded_control_constraint(prog_x, num_ux_x_pvar, den_ux_x_pvar, bound_Ax_row{i}, "posval", bound_lbx(i), monos_sos_x, safe_set_x_pvar, target_set_x_pvar, feasible_space);
        % for the upper bound constraint: bound_Ax(i, :) * ux <= bound_ubx(i)
        prog_x = add_bounded_control_constraint(prog_x, num_ux_x_pvar, den_ux_x_pvar, bound_Ax_row{i}, "negval", -1 * bound_ubx(i), monos_sos_x, safe_set_x_pvar, target_set_x_pvar, feasible_space);
    end

    % addtional constraints for the SOS program to ensure that all variables are bounded within a reasonable range for better numerical stability of the SOS solver, which can be tuned based on the specific problem settings
    % M = 1e6; % a large constant to define the reasonable range for the variables, which can be tuned based on the specific problem settings

    % for i = 1:length(x_vars_pvar)
    %     prog_x = sosineq(prog_x, M - x_vars_pvar(i)); % x_vars_pvar(i) <= M
    %     prog_x = sosineq(prog_x, M + x_vars_pvar(i)); % -x_vars_pvar(i) <= M
    % end

    % now the prog_x SOS program is ready to solve
    % set the options for the SOS solver
    sosoptions.solver = 'mosek';
    % set the objective for the SOS program, here we only care about to minimize the delta for the feasibility of the SOS program
    prog_x = sossetobj(prog_x, delta_dpvar);

    % solve the SOS program
    prog_x = sossolve(prog_x, sosoptions);

    % confirm the feasibility of the SOS program from the solver, like PRIMAL_AND_DUAL_FEASIBLE
    sol_info = prog_x.solinfo.info;

    if ~strcmp(sol_info, 'PRIMAL_AND_DUAL_FEASIBLE')
        error('The SOS program for bounded control is not feasible!');
    end

    % if the SOS program is feasible, then we can extract the optimal k1 polynomial and the optimal control input from the SOS program

    % get the optimal numerator polynomials from the SOS program
    num_ux_opt = cell(m, 1);

    for i = 1:m
        num_ux_opt{i} = sosgetsol(prog_x, num_ux_x_pvar{i});
    end

    % return the optimal control input as rational polynomial vector
    ux_opt = sym(zeros(m, 1));

    for i = 1:m
        ux_opt(i) = num_ux_opt{i} / den_ux_x_pvar{i};
    end

end

% helper function to handle each row of bounded control constraints and add them to the SOS program
function prog = add_bounded_control_constraint(prog, num_ux, den_ux, bound_Ax_row, bound_type, bound_val, monos_sos, safe_set, target_set, feasible_space)
    % add the bounded control constraints defined as
    % bound_Ax_row * ux >= bound_lbx_i, where ux = num_ux / den_ux over the domain (safe_set /target_set,means safe set minus target set)
    % this is equivalent to satisfying the following constraints simultaneously:

    sos_sum = 0; % sos multipliers for mataining the lower bound constraint

    % for each element in the row of bound_Ax
    for i = 1:length(bound_Ax_row)
        c = bound_Ax_row(i).col; % column index of the nonzero element in bound_Ax
        v = bound_Ax_row(i).(bound_type); % value of the nonzero element in bound_Ax, which can be a scalar or a polynomial

        % sos mulitplier for handling the lower bound constraint later
        [prog, s_i] = sospolyvar(prog, monos_sos);

        % add the rational polynomial constraint for this element to the SOS program
        prog = add_rational_poly_constraint(prog, v * num_ux{c}, den_ux{c}, s_i, monos_sos, safe_set, target_set, feasible_space);

        % store the sos multiplier for subsequent constraint construction
        sos_sum = sos_sum + s_i;

    end

    % add the constraint: sum(s_i) >= bound_val for all x in the domain (safe_set / target_set)
    % auxiliary SOS multiplier for this constraint
    [prog, sos_sum_s1] = sospolyvar(prog, monos_sos);
    [prog, sos_sum_s2] = sospolyvar(prog, monos_sos);
    prog = sosineq(prog, sos_sum - bound_val - sos_sum_s1 * safe_set - sos_sum_s2 * target_set);
    prog = sosineq(prog, sos_sum_s1);
    prog = sosineq(prog, sos_sum_s2);

end

% helper function to add rational polynomial constraint to the SOS program
function prog = add_rational_poly_constraint(prog, num_poly, den_poly, bound_poly, sos_monos, safe_set, target_set, feasible_space)
    % add the rational polynomial constraint: num_poly / den_poly >= bound_poly over the domain (safe_set / target_set)
    % CAUTION, based on the proof in the paper, only num_poly has the decision variables, while den_poly is some known polynomial without decision variables
    % this is equivalent to satisfying two separated constraints simultaneously:
    % 1) num_poly - den_poly * bound_poly >= 0 for all x in the domain that den_poly > 0
    % 2) num_poly - bound_poly * den_poly >= 0 for all x in the domain that den_poly < 0

    % sos multiplier for the first constraint
    [prog, sos_pos_1] = sospolyvar(prog, sos_monos);
    [prog, sos_pos_2] = sospolyvar(prog, sos_monos);
    [prog, sos_pos_3] = sospolyvar(prog, sos_monos);
    [prog, sos_pos_4] = sospolyvar(prog, sos_monos);
    [prog, sos_neg_1] = sospolyvar(prog, sos_monos);
    [prog, sos_neg_2] = sospolyvar(prog, sos_monos);
    [prog, sos_neg_3] = sospolyvar(prog, sos_monos);
    [prog, sos_neg_4] = sospolyvar(prog, sos_monos);

    % add the first constraint to the SOS program: num_poly - den_poly * bound_poly >= 0 for all x in the domain that den_poly > 0
    prog = sosineq(prog, num_poly - den_poly * bound_poly - sos_pos_1 * den_poly - sos_pos_2 * safe_set - sos_pos_3 * target_set - sos_pos_4 * feasible_space);
    prog = sosineq(prog, sos_pos_1);
    prog = sosineq(prog, sos_pos_2);
    prog = sosineq(prog, sos_pos_3);
    prog = sosineq(prog, sos_pos_4);

    % add the second constraint to the SOS program: num_poly - bound_poly * den_poly >= 0 for all x in the domain that den_poly < 0
    prog = sosineq(prog, num_poly - bound_poly * den_poly - sos_neg_1 * (-den_poly) - sos_neg_2 * safe_set - sos_neg_3 * target_set - sos_neg_4 * feasible_space);
    prog = sosineq(prog, sos_neg_1);
    prog = sosineq(prog, sos_neg_2);
    prog = sosineq(prog, sos_neg_3);
    prog = sosineq(prog, sos_neg_4);

end

% helper function to convert symbolic polynomial to pvar polynomial without dpvar
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
