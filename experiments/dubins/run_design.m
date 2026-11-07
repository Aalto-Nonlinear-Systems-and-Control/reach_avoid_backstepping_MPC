% k1 of Example 1 from the SOS program in scaled outputs y / SCALE, with lambda fixed (LMIN), |k1_i| <= KMAX and degree 8.
% Writes k1.txt: k1_1, k1_2, lambda, delta_s, alpha (string), pmax. Prints the solve time reported in Table I.
global LMIN KMAX SCALE
here = fileparts(mfilename('fullpath')); cd(here);
addpath(fullfile(here, '..', 'sos'));
syms y1 y2; y_vars_sym = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7)/0.5)^2 - 0.4;
alpha_sym = (5/2 - y1)^2 + (3/10)*y1*y2^2*(5/2 - y1);
SCALE = 2; KMAX = 1; LMIN = 0.06; d = 8;
[G1, G2] = meshgrid(linspace(-2.05, 2.05, 1641));
hf = matlabFunction(h_raw, 'Vars', [y1 y2]); Hg = hf(G1, G2);
af = matlabFunction(alpha_sym, 'Vars', [y1 y2]); Pg = af(G1, G2) .* Hg; pmax = max(Pg(Hg > 0));
safe_set_sym = alpha_sym * h_raw / pmax;
t0 = tic; [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, d, d);
fprintf('__DESIGN__,lambda %.10g,delta %.10g,delta/lambda %.10g,solve time %g s\n', lam, del, del/lam, toc(t0));
fid = fopen('k1.txt', 'w');
fprintf(fid, '%s\n%s\n%.17g\n%.17g\n%s\n%.17g\n', char(k1y(1)), char(k1y(2)), lam, del, char(alpha_sym), pmax); fclose(fid);
