clear;
echo off;
close all;
tic

pvar x y;
dpvar lambda;
vars = [x, y];

ds = 8;
du = 4;
dm = 6;
epsilon = 1e-8;

% safe set h(x)>0
h = 1 - x ^ 2 - y ^ 2;
h = 1 - (x - 0.5) ^ 2 - (y - 0.5) ^ 2;
h =- (x .^ 4 + y .^ 4 - 16) .* (x .^ 4 + y .^ 4 - 8);

% target set g(x)<0
g = 2 * ((y - 0.1) / 2) ^ 2 + 3 * ((x + 0.4) / 3) ^ 4 + (3 * x + 0.3) ^ 2 * (4 * y + 0.2) ^ 2 - 0.01;
g = (y - 1.2) ^ 2 + ((x + 1.8) / 0.5) ^ 2 - 0.1;

alpha = 1e-3 * (-g + 300);
h = alpha .* h;

%% sos optimization for constructing the control input u
prog = sosprogram(vars, lambda);

monos_s = monomials(vars, 0:ds);
[prog, s0] = sospolyvar(prog, monos_s);
[prog, s1] = sospolyvar(prog, monos_s);
[prog, s2] = sospolyvar(prog, monos_s);

% constraints of the sos program
prog = sosineq(prog, lambda * h + s0 * h + s1 * g + s2);
prog = sosineq(prog, lambda - epsilon);
prog = sosineq(prog, s0);
prog = sosineq(prog, s1);
prog = sosineq(prog, s2);

% set the solver
solver_opt.solver = 'mosek';
% solve the sos program
prog = sossolve(prog, solver_opt);

% obtain the solutions
lambda_poly = sosgetsol(prog, lambda);
s0_poly = sosgetsol(prog, s0);
s1_poly = sosgetsol(prog, s1);
s2_poly = sosgetsol(prog, s2);

% print the obtianed polynomials
disp('The obtained polynomial lambda is:');
disp(poly2str(lambda_poly));
disp('The obtained polynomial s0 is:');
disp(poly2str(s0_poly));
disp('The obtained polynomial s1 is:');
disp(poly2str(s1_poly));
disp('The obtained polynomial s2 is:');
disp(poly2str(s2_poly));

%% construction of the control input u
% compute the jacobian of h with respect to vars
dh1 = diff(h, vars(1));
dh2 = diff(h, vars(2));

dh_sq = dh1 ^ 2 + dh2 ^ 2;
u_rhs_poly = lambda_poly * h + s0_poly * h + s1_poly * g + s2_poly;

%% convert the polynomials to symbolic expressions for plotting
syms x y;
dh1_sym = str2sym(poly2str(dh1));
dh2_sym = str2sym(poly2str(dh2));
dh_sq_sym = str2sym(poly2str(dh_sq));
u_rhs_sym = str2sym(poly2str(u_rhs_poly));
P = u_rhs_sym * dh1_sym / dh_sq_sym;
Q = u_rhs_sym * dh2_sym / dh_sq_sym;

%% plotting
x_min = -2.1;
x_max = 2.1;
y_min = -2.1;
y_max = 2.1;

hold on; axis equal;

plotPolyL(-g, vars, [x_min, x_max], [y_min, y_max], [0 0.4470 0.7410], [0 0]);
plotPolyL(h, vars, [x_min, x_max], [y_min, y_max], [0.8500 0.3250 0.0980], [0 0]);

% plotVectorFieldDirection(P,Q,[x_min,y_min],[x_max, y_max],70);
plotStreamlines(P, Q, [x_min, x_max], [y_min, y_max], 500);

%% functionals
function plotStreamlines(P, Q, xrange, yrange, density, numStart)
    % Plot streamlines of a 2D vector field
    % Input parameters:
    %   P      : x component expression (function handle/string/symbolic)
    %   Q      : y component expression
    %   xrange : x axis range, default [-5,5]
    %   yrange : y axis range, default [-5,5]
    %   density: grid density, default 30
    %   numStart: number of starting points (per direction), default 15

    % Handle input parameters
    if nargin < 6, numStart = 15; end
    if nargin < 5, density = 30; end
    if nargin < 4, yrange = [-5 5]; end
    if nargin < 3, xrange = [-5 5]; end

    % Convert input to function handle
    P = expr2fun(P);
    Q = expr2fun(Q);

    % Generate computation grid
    x = linspace(xrange(1), xrange(2), density);
    y = linspace(yrange(1), yrange(2), density);
    [X, Y] = meshgrid(x, y);

    % Compute vector field
    U = P(X, Y);
    V = Q(X, Y);

    % check if there are NaN values in U or V
    if any(isnan(U(:))) || any(isnan(V(:)))
        error('The vector field contains NaN values. Please check the input expressions P and Q.');
    end

    % Generate starting points
    startx = linspace(xrange(1) * 0.9, xrange(2) * 0.9, numStart);
    starty = linspace(yrange(1) * 0.9, yrange(2) * 0.9, numStart);
    [startX, startY] = meshgrid(startx, starty);

    % Compute streamlines
    s = streamslice(X, Y, U, V, 3);
    set(s, 'Color', [146/255 149/255 145/255])

    % Figure formatting
    axis([xrange yrange])
    xlabel('y1'), ylabel('y2')
    title('')
    grid on
    set(gca, 'Box', 'on')
    hold off
end

% Helper function: convert various inputs to function handle
function fh = expr2fun(expr)

    if isa(expr, 'function_handle')
        fh = expr;
    elseif ischar(expr)
        fh = str2func(['@(x,y)' vectorize(expr)]);
    elseif isa(expr, 'sym')
        vars = symvar(expr);

        if numel(vars) < 2
            vars = [sym('x'), sym('y')]; % Force x,y as variables
        end

        fh = matlabFunction(expr, 'Vars', vars(1:2));
    else
        error('Unsupported input type');
    end

end

function plotVectorFieldDirection(P, Q, xrange, yrange, density)
    % Plot 2D vector field (showing only direction)
    % Note: Normalize the vector field so all arrows have the same length

    % Handle input type for P
    if ischar(P)
        P = str2func(['@(x,y)' vectorize(P)]);
    elseif isa(P, 'sym')
        vars = symvar(P);

        if numel(vars) < 2
            vars = [sym('x'), sym('y')]; % Default variables x,y
        end

        P = matlabFunction(P, 'Vars', vars(1:2));
    end

    % Handle input type for Q
    if ischar(Q)
        Q = str2func(['@(x,y)' vectorize(Q)]);
    elseif isa(Q, 'sym')
        vars = symvar(Q);

        if numel(vars) < 2
            vars = [sym('x'), sym('y')]; % Default variables x,y
        end

        Q = matlabFunction(Q, 'Vars', vars(1:2));
    end

    % Generate grid
    x = linspace(xrange(1), xrange(2), density);
    y = linspace(yrange(1), yrange(2), density);
    [X, Y] = meshgrid(x, y);

    % Compute original vector components
    U = P(X, Y);
    V = Q(X, Y);

    stream2(X, Y, U, V, X, Y);

    axis tight;
    xlabel('x'); ylabel('y');
    title('Normalized vector field direction');
    grid on;
end

function s = poly2str(p)
    c = char(p);
    s = c{1};
end
