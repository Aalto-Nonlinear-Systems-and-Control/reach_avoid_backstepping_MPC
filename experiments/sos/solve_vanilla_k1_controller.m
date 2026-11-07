% [test] wrapper: solve the vanilla k1 SOS in scaled outputs ys = y / SCALE, map back
function [k1_opt, k1_lambda, k1_delta] = solve_vanilla_k1_controller(y_vars, safe_set, target_set, dv, ds)
    global SCALE KMAX
    sc = SCALE; if isempty(sc) || sc <= 0, sc = 1; end
    K0 = KMAX;
    if ~isempty(KMAX) && KMAX > 0, KMAX = KMAX / sc; end
    ss = expand(subs(safe_set, y_vars, sc * y_vars));
    ts = expand(subs(target_set, y_vars, sc * y_vars));
    [k1s, k1_lambda, k1_delta] = solve_vanilla_k1_core(y_vars, ss, ts, dv, ds);
    KMAX = K0;
    k1_opt = expand(sc * subs(k1s, y_vars, y_vars / sc));
    fprintf('__WRAP__,scale=%g,lambda=%.10g,delta=%.10g\n', sc, k1_lambda, k1_delta);
end
