% k1 for one kitting task of the top-view arm: psi from tasks/<TASK>.json (normalised to max 1), target disc at the slot.
% Solved in scaled outputs y / SCALE with lambda fixed (LMIN) and |k1_i| <= KMAX on {psi >= 0}; degree DEG.
global LMIN KMAX SCALE
here = fileparts(mfilename('fullpath')); cd(here);
addpath(fullfile(here, '..', '..', '..', 'matlab')); addpath(fullfile(here, '..'));
if ~exist('TASK', 'var'), TASK = 'red_fwd'; end
if ~exist('DESIGN_LAMS', 'var'), DESIGN_LAMS = [0.02, 0.05, 0.1]; end
if ~exist('DEG', 'var'), DEG = 8; end
spec = jsondecode(fileread(fullfile('tasks', [TASK '.json'])));
syms y1 y2; y_vars_sym = [y1; y2];
psi_raw = str2sym(strrep(spec.psi, '**', '^')); pmax = spec.pmax;
safe_set_sym = psi_raw / pmax;
target_set_sym = (y1 - spec.target(1))^2 + (y2 - spec.target(2))^2 - spec.RT^2;
SCALE = 4; KMAX = 1; if exist('DESIGN_K', 'var'), KMAX = DESIGN_K; end
for LMIN = DESIGN_LAMS
    try
        t0 = tic; [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, DEG, DEG);
        fprintf('__DESIGN__,%s,%g,%g,%.10g,%.10g,%.10g,%g\n', TASK, DEG, LMIN, lam, del, del/lam, toc(t0));
        fid = fopen(sprintf('k1/%s_d%d_l%g.txt', TASK, DEG, LMIN), 'w');
        fprintf(fid, '%s\n%s\n%.17g\n%.17g\n%s\n%.17g\n', char(k1y(1)), char(k1y(2)), lam, del, char(psi_raw), pmax); fclose(fid);
    catch err
        fprintf('__ERR__,%s,%g,%g,%s\n', TASK, DEG, LMIN, strrep(err.message, newline, ' '));
    end
end
