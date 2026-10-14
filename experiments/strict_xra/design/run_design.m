% Rigorous k1 for a list of weights and lambdas (scaled, |k1_i| <= KMAX, degree 8).
% Writes k1/<tag>_l<lambda>.txt: k1_1, k1_2, lambda, delta_s, alpha (string), pmax.
global LMIN KMAX SCALE
here = fileparts(mfilename('fullpath')); cd(here);
addpath(fullfile(here, '..', '..', '..', 'matlab')); addpath(fullfile(here, '..'));
syms y1 y2; y_vars_sym = [y1; y2];
h_raw = -(y1^4 + y2^4 - 16) * (y1^4 + y2^4 - 4);
target_set_sym = (y2 - 0)^2 + ((y1 + 1.7)/0.5)^2 - 0.4;
W = struct('tag', {'W0', 'W3'}, 'alpha', {(5/2 - y1)^2, (5/2 - y1)^2 + (3/10)*y1*y2^2*(5/2 - y1)});
lams = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2];
if exist('DESIGN_TAGS', 'var'), W = W(ismember({W.tag}, DESIGN_TAGS)); end
if exist('DESIGN_LAMS', 'var'), lams = DESIGN_LAMS; end
SCALE = 2; KMAX = 1; d = 8;
[G1, G2] = meshgrid(linspace(-2.05, 2.05, 1641));
hf = matlabFunction(h_raw, 'Vars', [y1 y2]); Hg = hf(G1, G2);
for w = 1:numel(W)
    af = matlabFunction(W(w).alpha, 'Vars', [y1 y2]); Pg = af(G1, G2) .* Hg; pmax = max(Pg(Hg > 0));
    safe_set_sym = W(w).alpha * h_raw / pmax;
    for LMIN = lams
        try
            t0 = tic; [k1y, lam, del] = solve_vanilla_k1_controller(y_vars_sym, safe_set_sym, target_set_sym, d, d);
            fprintf('__DESIGN__,%s,%g,%.10g,%.10g,%.10g,%g\n', W(w).tag, LMIN, lam, del, del/lam, toc(t0));
            fid = fopen(sprintf('k1/%s_l%g.txt', W(w).tag, LMIN), 'w');
            fprintf(fid, '%s\n%s\n%.17g\n%.17g\n%s\n%.17g\n', char(k1y(1)), char(k1y(2)), lam, del, char(W(w).alpha), pmax); fclose(fid);
        catch err
            fprintf('__ERR__,%s,%g,%s\n', W(w).tag, LMIN, strrep(err.message, newline, ' '));
        end
    end
end
