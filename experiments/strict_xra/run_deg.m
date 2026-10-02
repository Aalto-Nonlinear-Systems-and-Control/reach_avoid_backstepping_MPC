global LMIN
cd(.); addpath('/Users/jianqiang/Downloads/reach_avoid_backstepping_MPC/matlab'); addpath(.);
syms y1 y2; y_vars_sym = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7)/0.5)^2 - 0.4;
safe_set_sym = 1e-3 * (-target_set_sym + 300) * h_raw;
for LMIN = [1e-8, 1]
  for d = [4 6 8]
    try
      t0 = tic; [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, d, d);
      fprintf('__DEG__,%g,%d,%g,%g,%g,%g\n', LMIN, d, lam, del, del/lam, toc(t0));
    catch err
      fprintf('__ERR__,%g,%d,%s\n', LMIN, d, strrep(err.message, newline, ' '));
    end
  end
end
