% k1 of Example 2: psi from task.json (normalised to max 1), target disc at the slot. Solved in scaled outputs y / SCALE with
% lambda fixed (LMIN), |k1_i| <= KMAX on {psi >= 0} and degree DEG. Writes k1.txt. Writes the time of this script to data/sos_time.txt for Table I.
t_total = tic;
global LMIN KMAX SCALE
here = fileparts(mfilename('fullpath')); cd(here);
addpath(fullfile(here, '..', 'sos'));
spec = jsondecode(fileread('task.json'));
syms y1 y2; y_vars_sym = [y1; y2];
psi_raw = str2sym(strrep(spec.psi, '**', '^')); pmax = spec.pmax;
safe_set_sym = psi_raw / pmax;
target_set_sym = (y1 - spec.target(1))^2 + (y2 - spec.target(2))^2 - spec.RT^2;
SCALE = 4; KMAX = 1; LMIN = 0.1; DEG = 8;
t0 = tic; [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, DEG, DEG);
fprintf('__DESIGN__,lambda %.10g,delta %.10g,delta/lambda %.10g,solve time %g s\n', lam, del, del/lam, toc(t0));
fid = fopen('k1.txt', 'w');
fprintf(fid, '%s\n%s\n%.17g\n%.17g\n%s\n%.17g\n', char(k1y(1)), char(k1y(2)), lam, del, char(psi_raw), pmax); fclose(fid);
fid = fopen(fullfile('data', 'sos_time.txt'), 'w'); fprintf(fid, '%.3f\n', toc(t_total)); fclose(fid);
