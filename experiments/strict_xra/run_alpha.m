% [test] steeper weight alpha = (2.5 - y1)^2, psi normalised by 1/600, scaled rigorous k1
global LMIN KMAX SCALE
here = fileparts(mfilename('fullpath')); cd(here); addpath(fullfile(here, '..', '..', 'matlab')); addpath(here);
syms y1 y2; y_vars_sym = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7)/0.5)^2 - 0.4;
safe_set_sym = (5/2 - y1)^2 * h_raw / 600;
SCALE = 2; KMAX = 1;
for d = [8]
  for LMIN = [0.0005, 0.001, 0.002, 0.005, 0.01, 0.02]
    try
        t0 = tic; [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, d, d);
        fprintf('__ALPHA__,%d,%g,%.10g,%.10g,%.10g,%g\n', d, LMIN, lam, del, del/lam, toc(t0));
        fid = fopen(sprintf('out/alpha_d%d_l%g_k1.txt', d, LMIN), 'w');
        fprintf(fid, '%s\n%s\n%.17g\n%.17g\n', char(k1y(1)), char(k1y(2)), lam, del); fclose(fid);
    catch err
        fprintf('__ERR__,%d,%g,%s\n', d, LMIN, strrep(err.message, newline, ' '));
    end
  end
end
