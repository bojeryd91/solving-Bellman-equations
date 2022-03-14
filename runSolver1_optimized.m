function [X, curr_V, optimal_c, iter, err] = runSolver1_optimized(init_V)
%% Define economic parameters
beta = 0.95;
r = 1.0/beta - 1.0 - 0.01;
y = 1.0;
c_min = 0.0001;
N = length(init_V);

%% Initialize vectors
X = linspace(1.0, 4.0, N)';
curr_V = init_V;
optimal_c = zeros(N, 1);

%% The solver
iter = 0; iter_max = 1000; err = 1.0e10; err_tol = 1.0e-4;
while iter < iter_max && err > err_tol
    % For a particular x, compute the right-hand-side of the
    %   Bellman equation for each possible choice x' in X.
    Cs = X - (X' - y)./(1.0+r);
    u  = (Cs >= c_min).*log(Cs) + (Cs < c_min).*(-1.0e8);
    Bellman_rhs = u + beta.*curr_V';

    % Identify the maximum and maximizer
    [maxim, index_maximizer] = max(Bellman_rhs, [], 2);
    new_V = maxim;
    optimal_xp = X(index_maximizer);
    optimal_c = X - (optimal_xp - y)/(1+r);

    % Prepare for next iteration
    err    = max(abs(new_V - curr_V));
    curr_V = new_V;
    iter = iter + 1;
end
end

