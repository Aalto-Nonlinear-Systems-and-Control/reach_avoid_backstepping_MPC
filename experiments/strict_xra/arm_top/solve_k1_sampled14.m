% Study: program (14) solved directly with the input bounds enforced at sampled states (top-view arm).
% lambda fixed; ECGBF condition as in (12) (SOS in scaled outputs y / SIGMA); optional |k1_i| <= K on {psi >= 0};
% sampled input bounds as linear constraints G * theta <= b on the coefficient vector theta of k1 in the monomial
% basis MONS of the scaled outputs (k1(y) = SIGMA * sum_k theta_k (y1/SIGMA)^p_k (y2/SIGMA)^q_k).
% Input: IN_MAT with fields psi_str, pmax, target, RT, lam, kmax (0: no bound), sigma, deg, mons, G, b, out_mat.
% Output: out_mat with theta, delta, solver info. Does not touch any existing design file.
here = fileparts(mfilename('fullpath')); cd(here);
addpath(fullfile(here, '..', '..', '..', 'matlab')); addpath('/Users/jianqiang/Downloads/SOSTOOLS-4.03');
D = load(IN_MAT); t0 = tic; D.deg = double(D.deg); D.mons = double(D.mons);
syms y1 y2; Y = [y1; y2]; sg = D.sigma;
psi = str2sym(strrep(D.psi_str, '**', '^')) / D.pmax; phi = (y1 - D.target(1))^2 + (y2 - D.target(2))^2 - D.RT^2;
psis = expand(subs(psi, Y, sg * Y)); phis = expand(subs(phi, Y, sg * Y));
pvar z1 z2; Z = [z1; z2]; psp = sym2pvar(psis, Y, Z); php = sym2pvar(phis, Y, Z);
nb = size(D.mons, 1); th = []; names = {};
for i = 1:2, for k = 1:nb, th = [th; dpvar(sprintf('th%d_%d', i, k))]; end, end
dpvar delta
prog = sosprogram(Z, [delta; th]);
k1a = 0; k1b = 0;
for k = 1:nb, mk = z1^D.mons(k, 1) * z2^D.mons(k, 2); k1a = k1a + th(k) * mk; k1b = k1b + th(nb + k) * mk; end
k1 = [k1a; k1b];
ms = monomials(Z, 0:D.deg);
[prog, s1] = sospolyvar(prog, ms); [prog, s2] = sospolyvar(prog, ms);
lie = diff(psp, z1) * k1(1) + diff(psp, z2) * k1(2);
prog = sosineq(prog, lie - D.lam * psp + delta - s1 * psp - s2 * php);
prog = sosineq(prog, s1); prog = sosineq(prog, s2); prog = sosineq(prog, delta);
if D.kmax > 0
    Kt = D.kmax / sg;
    for i = 1:2
        [prog, sa] = sospolyvar(prog, ms); [prog, sb] = sospolyvar(prog, ms);
        prog = sosineq(prog, sa); prog = sosineq(prog, sb);
        prog = sosineq(prog, Kt - k1(i) - sa * psp); prog = sosineq(prog, Kt + k1(i) - sb * psp);
    end
end
e = D.G * th;
for r = 1:numel(D.b), prog = sosineq(prog, D.b(r) - e(r)); end
fprintf('__BUILD__,%d constraints,%.1f s\n', numel(D.b), toc(t0));
prog = sossetobj(prog, delta);
opt.solver = 'mosek'; prog = sossolve(prog, opt); si = prog.solinfo.info;
theta = double(sosgetsol(prog, th)); dl = double(sosgetsol(prog, delta));
info = [si.pinf, si.dinf, si.numerr, si.feasratio, si.cpusec, si.iter];
fprintf('__SOL__,%g,%g,%g,%g,%g,%g,delta=%.6g,total %.1f s\n', info, dl, toc(t0));
save(D.out_mat, 'theta', 'dl', 'info');
