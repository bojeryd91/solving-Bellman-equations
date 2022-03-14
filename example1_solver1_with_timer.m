%% Clean up
clear, clc, close all % Do NOT use "clear all"

tic()

%% Define economic parameters
beta = 0.95;
r = 1.0/beta - 1.0 - 0.01;
y = 1.0;
c_min = 0.0001;
u =@(c) (c >= c_min).*log(c) + (c < c_min).*(-1.0e8);

%% Define solver parameters
N = 100;

%% Initialize vectors
X = linspace(1.0, 4.0, N)';
curr_V = ones(N, 1);
new_V  = curr_V;
optimal_c = zeros(N, 1);

%% The solver
iter = 0; iter_max = 1000; err = 1.0e10; err_tol = 1.0e-4;
while iter < iter_max && err > err_tol
    for i_x = 1:N
        x = X(i_x);
        % For a particular x, compute the right-hand-side of the
        %   Bellman equation for each possible choice x' in X.
        Bellman_rhs = u(x - (X - y)./(1.0+r)) + beta.*curr_V;

        % Identify the maximum and maximizer
        [maxim, index_maximizer] = max(Bellman_rhs);
        new_V(i_x) = maxim;
        optimal_xp = X(index_maximizer);
        optimal_c(i_x) = x - (optimal_xp - y)/(1+r);
    end
    % Prepare for next iteration
    err    = max(abs(new_V - curr_V));
    curr_V = new_V;
    iter = iter + 1;
end

toc()

%% Output results
fprintf("Number iterations to convergence: %i\n", iter)
fprintf("Final maximum error: %8.5g\n", err)

plot(X, optimal_c); xlabel("Cash-on-hand"); ylabel("Consumption")
title("Optimal consumption by cash-on-hand")
figure()
plot(X, curr_V); xlabel("Cash-on-hand"); ylabel("Value function")
title("Value function by cash-on-hand")

%  Compute saving decision and plot it
optimal_s = X - optimal_c;
figure()
plot(X, optimal_s); xlabel("Cash-on-hand"); ylabel("Savings")
title("Optimal savings by cash-on-hand")
