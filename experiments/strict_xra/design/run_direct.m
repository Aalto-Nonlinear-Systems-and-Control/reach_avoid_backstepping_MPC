% Check: solving the k1 SOS directly at the final lambda (and |k1_i| <= 0.159) reproduces the time-scaled design.
global LMIN KMAX SCALE
here = fileparts(mfilename('fullpath')); cd(here);
addpath(fullfile(here, '..', '..', '..', 'matlab')); addpath(fullfile(here, '..'));
syms y1 y2; y_vars_sym = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7)/0.5)^2 - 0.4;
alpha = (5/2 - y1)^2 + (3/10)*y1*y2^2*(5/2 - y1);
[G1, G2] = meshgrid(linspace(-2.05, 2.05, 1641));
hf = matlabFunction(h_raw, 'Vars', [y1 y2]); Hg = hf(G1, G2);
af = matlabFunction(alpha, 'Vars', [y1 y2]); Pg = af(G1, G2) .* Hg; pmax = max(Pg(Hg > 0));
safe_set_sym = alpha * h_raw / pmax;
s = str2double(getenv("DESIGN_S")); if isnan(s), s = 0.159; end; SCALE = 2; KMAX = s; d = 8; LMIN = s * 0.060109801660967063;
t0 = tic; [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, d, d);
fprintf('__DIRECT__,%.10g,%.10g,%.10g,%g\n', LMIN, lam, del, toc(t0));
fid = fopen(sprintf("k1/W3_direct_s%g.txt", s), "w");
fprintf(fid, '%s\n%s\n%.17g\n%.17g\n%s\n%.17g\n', char(k1y(1)), char(k1y(2)), lam, del, char(alpha), pmax); fclose(fid);
