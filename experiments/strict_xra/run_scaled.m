global LMIN KMAX
here = fileparts(mfilename('fullpath')); cd(here); addpath(fullfile(here, '..', '..', 'matlab')); addpath(here);
syms y1 y2; y_vars_sym = [y1; y2]; sc = 2;   % y = sc * ys, solve in ys (named y1,y2 here)
Y1 = sc*y1; Y2 = sc*y2;
h_raw = -(Y1^4 + Y2^4 - 16) * (Y1^4 + Y2^4 - 4);
target_set_sym = expand((Y2 - 0)^2 + ((Y1 + 1.7)/0.5)^2 - 0.4);
safe_set_sym = expand(1e-3 * (-target_set_sym + 300) * h_raw);
KMAX = 1/sc;
for d = [6 8]
for LMIN = [0.01, 0.02, 0.05, 0.1]
    try
        t0 = tic; [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, d, d);
        fprintf('__SCALED__,%d,%g,%.10g,%.10g,%.10g,%g\n', d, LMIN, lam, del, del/lam, toc(t0));
        fid = fopen(sprintf('out/scaled_d%d_l%g_k1.txt', d, LMIN), 'w');
        fprintf(fid, '%s\n%s\n%.17g\n%.17g\n', char(k1y(1)), char(k1y(2)), lam, del); fclose(fid);
    catch err
        fprintf('__ERR__,%g,%s\n', LMIN, strrep(err.message, newline, ' '));
    end
end
end
