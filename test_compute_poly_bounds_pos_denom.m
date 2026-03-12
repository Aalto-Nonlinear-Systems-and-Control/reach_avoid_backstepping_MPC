% Test script for computing polynomial bounds with positive denominator
clc;
clear;
close all;

syms x1 x2;
num_poly = 0.01 - x1 ^ 2 - x2 ^ 2; % general numerator polynomial (some part of it is positive, but not necessarily positive everywhere)
% den_poly = x1 ^ 2 + x2 ^ 2 + 1; % positive denominator everywhere
den_poly = x1 ^ 2 + x2 ^ 2 - 0.04;
superlevel_set_poly = 1 - x1 ^ 2 - x2 ^ 2; % superlevel set is the unit disk
ds_min = 6; % minimum degree for SOS multipliers
epsilon = 1e-4;

[lower_bound, upper_bound] = compute_poly_bounds(num_poly, den_poly, superlevel_set_poly, ds_min, epsilon);

% using sampling to verify the bounds
n_samples = 50000;
bounds = [-1, 1; -1, 1]; % sampling bounds for x1 and x2
x1_samples = bounds(1, 1) + (bounds(1, 2) - bounds(1, 1)) * rand(n_samples, 1);
x2_samples = bounds(2, 1) + (bounds(2, 2) - bounds(2, 1)) * rand(n_samples, 1);

% conver the super_level set polynomial to a function handle for efficient evaluation
superlevel_set_func = matlabFunction(superlevel_set_poly, 'Vars', [x1, x2]);
num_poly_func = matlabFunction(num_poly, 'Vars', [x1, x2]);
den_poly_func = matlabFunction(den_poly, 'Vars', [x1, x2]);

% get the indices of samples that are inside the superlevel set (feasible region)
feasible_indices = find(superlevel_set_func(x1_samples, x2_samples) >= 0 ...
    & abs(den_poly_func(x1_samples, x2_samples)) >= epsilon); % also exclude the samples where the denominator is close to zero
% evaluate the rational polynomial at the feasible samples
rational_poly_values = num_poly_func(x1_samples(feasible_indices), x2_samples(feasible_indices)) ...
    ./ den_poly_func(x1_samples(feasible_indices), x2_samples(feasible_indices));

% get the minimum and maximum values of the rational polynomial at the feasible samples
min_rational_poly = min(rational_poly_values);
max_rational_poly = max(rational_poly_values);

% disp the results
disp(['Computed lower bound: ', num2str(lower_bound)]);
disp(['Computed upper bound: ', num2str(upper_bound)]);
disp(['Minimum value of rational polynomial at feasible samples: ', num2str(min_rational_poly)]);
disp(['Maximum value of rational polynomial at feasible samples: ', num2str(max_rational_poly)]);
