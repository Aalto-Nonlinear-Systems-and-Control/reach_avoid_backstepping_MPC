global LMIN KMAX SCALE
cd(.); addpath('/Users/jianqiang/Downloads/reach_avoid_backstepping_MPC/matlab'); addpath(.);
syms x1 x2 th v y1 y2;
x_vars_sym = [x1; x2; th; v]; fx_sym = [v*cos(th); v*sin(th); 0; 0];
gx_sym = [0 0; 0 0; 1 0; 0 1]; hx_sym = [x1; x2]; y_vars_sym = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7)/0.5)^2 - 0.4;
safe_set_sym = 1e-3 * (-target_set_sym + 300) * h_raw;
[u, k1, J_k1, mu, lambda, certificate, cert_term_dict, A_matrix, b_vector, ks, p, r_deg] = ...
    reach_avoid_controller(fx_sym, gx_sym, hx_sym, x_vars_sym, y_vars_sym, safe_set_sym);
lb = [-5; -5]; ub = [5; 5]; samples_num = 3000;
bound_min = [-2; -2; 2*pi/3; 0.1]; bound_max = [2; 2; 4*pi/3; 1.0];
SCALE = 2; KMAX = 1; LMIN = 0.02; d = 8;
mus = [0.1, 1];
for i = 1:numel(mus)
    mu_val = mus(i); rng(42);
    fprintf('__CASE__,mu=%g\n', mu_val);
    try
        t0 = tic;
        [u_opt, cert_opt, valid_count, k1_opt] = solvesop_bounded_control(u, k1, J_k1, mu, lambda, certificate, cert_term_dict, p, r_deg, ...
            x_vars_sym, y_vars_sym, hx_sym, safe_set_sym, target_set_sym, mu_val, lb, ub, d, d, samples_num, bound_min, bound_max, ...
            sprintf('tmp_full_%d_uncon.py', i));
        fprintf('__SOPRES__,%d,%g\n', valid_count, toc(t0));
        params = struct(); params.lb = lb; params.ub = ub; params.mu_val = mu_val; params.lambda = LMIN; params.degree = d;
        export_to_python(u_opt, cert_opt, k1_opt, params, sprintf('tmp_full_%d_result.py', i));
    catch err
        fprintf('__ERR__,%s\n', strrep(err.message, newline, ' '));
    end
end
