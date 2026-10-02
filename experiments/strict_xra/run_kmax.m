global LMIN KMAX
here = fileparts(mfilename('fullpath')); cd(here); addpath(fullfile(here, '..', '..', 'matlab')); addpath(here);
syms y1 y2; y_vars_sym = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7)/0.5)^2 - 0.4;
safe_set_sym = 1e-3 * (-target_set_sym + 300) * h_raw;
KMAX = 1; d = 6;
for LMIN = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
    try
        t0 = tic; [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, d, d);
        fprintf('__KMAX__,%g,%g,%.10g,%.10g,%.10g,%g\n', KMAX, LMIN, lam, del, del/lam, toc(t0));
        fid = fopen(sprintf('out/kmax_l%g_k1.txt', LMIN), 'w');
        fprintf(fid, '%s\n%s\n%.17g\n%.17g\n', char(k1y(1)), char(k1y(2)), lam, del); fclose(fid);
    catch err
        fprintf('__ERR__,%g,%s\n', LMIN, strrep(err.message, newline, ' '));
    end
end
