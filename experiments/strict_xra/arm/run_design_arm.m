% k1 for the two-link arm: psi = ((1 + 0.6 (y2 + 2.5)) b - 0.02 |y|^4) / max, b the band of the paper.
% Solved in scaled outputs y / SCALE with lambda fixed (LMIN) and |k1_i| <= KMAX (output units), degree 8.
global LMIN KMAX SCALE
here = fileparts(mfilename('fullpath')); cd(here);
addpath(fullfile(here, '..', '..', '..', 'matlab')); addpath(fullfile(here, '..')); addpath(fullfile(here, '..', 'design'));
syms y1 y2; y_vars_sym = [y1; y2];
b = -((4 * (y1 - 2) - 2 * y2^3)^2) + (4/5) * y2^3 + 10;
psi_raw = (1 + (3/5) * (y2 + 5/2)) * b - (1/50) * (y1^2 + y2^2)^2;
target_set_sym = (y1 - 11/2)^2 / (36/25) + (y2 - 9/5)^2 / (4/25) - 2;
[G1, G2] = meshgrid(linspace(-6, 9, 1501), linspace(-4, 4, 801));
pf = matlabFunction(psi_raw, 'Vars', [y1 y2]); Pg = pf(G1, G2); pmax = max(Pg(:));
safe_set_sym = psi_raw / pmax;
lams = [0.02, 0.05, 0.1, 0.2];
if exist('DESIGN_LAMS', 'var'), lams = DESIGN_LAMS; end
SCALE = 4; KMAX = 1; d = 8; if exist("DESIGN_K", "var"), KMAX = DESIGN_K; end
for LMIN = lams
    try
        t0 = tic; [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, d, d);
        fprintf('__DESIGN__,ARM,%g,%.10g,%.10g,%.10g,%g\n', LMIN, lam, del, del/lam, toc(t0));
        fid = fopen(sprintf("k1/ARM_K%g_l%g.txt", KMAX, LMIN), 'w');
        fprintf(fid, '%s\n%s\n%.17g\n%.17g\n%s\n%.17g\n', char(k1y(1)), char(k1y(2)), lam, del, char(psi_raw), pmax); fclose(fid);
    catch err
        fprintf('__ERR__,ARM,%g,%s\n', LMIN, strrep(err.message, newline, ' '));
    end
end
