%clear;
echo off;
close all;
tic

pvar x y;
dpvar delta lambda;
vars = [x, y];

ds = 11;
du = 5;
factor = 1e3;

monos_ux1 = monomials([2 * (2 * rand() - 1) * x, 2 * (2 * rand() - 1) * y], 0:du);
monos_ux2 = monomials([2 * (2 * rand() - 1) * x, 2 * (2 * rand() - 1) * y], 0:du);

R = 3.5;
a = 2;
b = 1.5;

% safe set h(x)>0
%h = a * (R - y) ^ 2 - 1.5 * x - (x ^ 4 + y ^ 4 - R ^ 2) ^ 2;
h =- (x .^ 4 + y .^ 4 - 16) .* (x .^ 4 + y .^ 4 - 8);
% target set g(x)<0
g = (y - 1.2) ^ 2 + ((x + 1.8) / 0.5) ^ 2 - 0.1;
alpha = 1e-3 * (-g + 300);
h = alpha .* h;
%% sos optimization

prog = sosprogram(vars, [lambda, delta]);

[prog, u1] = sospolyvar(prog, monos_ux1);
[prog, u2] = sospolyvar(prog, monos_ux2);

Lhu = diff(h, vars(1)) * u1 + diff(h, vars(2)) * u2;

% auxiliary polynomials

monos_s = monomials(vars, 0:ds);

[prog, s0] = sospolyvar(prog, monos_s);
[prog, s1] = sospolyvar(prog, monos_s);

% constraints

prog = sosineq(prog, Lhu - lambda * h + delta - s0 * h - s1 * g);
prog = sosineq(prog, delta);
prog = sosineq(prog, lambda - factor * delta);
prog = sosineq(prog, s0);
prog = sosineq(prog, s1);

% set objective and solver
prog = sossetobj(prog, delta);
solver_opt.solver = 'mosek';

% solve
prog = sossolve(prog, solver_opt);

% obtain solutions
u1_poly = sosgetsol(prog, u1);
u2_poly = sosgetsol(prog, u2);
lambda_poly = sosgetsol(prog, lambda);
delta_poly = sosgetsol(prog, delta);

%% plotting

syms x y;
P = str2sym(poly2str(u1_poly));
Q = str2sym(poly2str(u2_poly));

x_min = -2.1;
x_max = 2.1;
y_min = -2.2;
y_max = 2;

hold on; axis equal;

plotPolyL(-g, vars, [x_min, x_max], [y_min, y_max], [0 0.4470 0.7410], [0 0]);
plotPolyL(h, vars, [x_min, x_max], [y_min, y_max], [0.8500 0.3250 0.0980], [0 0]);

plotStreamlines(P, Q, [x_min, x_max], [y_min, y_max], 100);

disp("u1_poly:----------------------------------------------------------")
disp(poly2str(u1_poly));
disp("u2_poly:----------------------------------------------------------")
disp(poly2str(u2_poly));
disp("lambda:----------------------------------------------------------")
disp(poly2str(lambda_poly));
disp("delta:----------------------------------------------------------")
disp(poly2str(delta_poly));
disp("delta/lambda:----------------------------------------------------------")
disp(poly2str(delta_poly / lambda_poly));

%% functional
function plotStreamlines(P, Q, xrange, yrange, density, numStart)
    % Plot streamlines of a 2D vector field
    % Input parameters:
    %   P      : x-component expression (function handle/string/symbolic)
    %   Q      : y-component expression
    %   xrange : x-axis range, default [-5,5]
    %   yrange : y-axis range, default [-5,5]
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

    % Generate starting points
    startx = linspace(xrange(1) * 0.9, xrange(2) * 0.9, numStart);
    starty = linspace(yrange(1) * 0.9, yrange(2) * 0.9, numStart);
    [startX, startY] = meshgrid(startx, starty);

    % Compute streamlines
    s = streamslice(X, Y, U, V, 3);
    set(s, 'Color', [146/255 149/255 145/255])

    % Figure decoration
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
            vars = [sym('x'), sym('y')]; % Force use x,y as variables
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
