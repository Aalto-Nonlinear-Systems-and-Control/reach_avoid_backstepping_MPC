global LMIN
S = .; R = '/Users/jianqiang/Downloads/reach_avoid_backstepping_MPC/matlab';
cd(S); addpath(R); addpath(S);   % S first on path: shadows the two SOS helpers
lmins = [1e-8, 1e-4, 1e-2, 1e-1, 1, 5];
syms x1 x2 th v y1 y2;
x_vars_sym = [x1; x2; th; v]; fx_sym = [v*cos(th); v*sin(th); 0; 0];
gx_sym = [0 0; 0 0; 1 0; 0 1]; hx_sym = [x1; x2]; y_vars_sym = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7)/0.5)^2 - 0.4;
alpha = 1e-3 * (-target_set_sym + 300); safe_set_sym = alpha * h_raw;
[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = ...
    reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym);
lb = [-5; -5]; ub = [5; 5]; ds = 4; dv = 4; mu_val = 0.1; samples_num = 1000;
bound_min = [-2; -2; 2*pi/3; 0.1]; bound_max = [2; 2; 4*pi/3; 1.0];
for i = 1:numel(lmins)
    LMIN = lmins(i); rng(42);
    fprintf('__CASE__,%g\n', LMIN);
    try
        t0 = tic;
        [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, dv, ds);
        c = double(coeffs(expand(sum(abs(k1y)))));  % rough size of k1
        fprintf('__VANRES__,%g,%g,%g,%g,%g\n', lam, del, del/lam, max(abs(double(coeffs(expand(k1y(1)))))), toc(t0));
        t0 = tic;
        uncon = sprintf('tmp_lmin_%d_uncon.py', i);
        [u_opt, cert_opt, valid_count, k1_opt] = solvesop_bounded_control(u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, ...
            x_vars_sym, y_vars_sym, hx_sym, safe_set_sym, target_set_sym, mu_val, lb, ub, ds, dv, samples_num, bound_min, bound_max, uncon);
        fprintf('__SOPRES__,%d,%g\n', valid_count, toc(t0));
        params = struct(); params.lb = lb; params.ub = ub; params.mu_val = mu_val; params.lambda = lam;
        export_to_python(u_opt, cert_opt, k1_opt, params, sprintf('tmp_lmin_%d_result.py', i));
    catch err
        fprintf('__ERR__,%s\n', strrep(err.message, newline, ' '));
    end
end
