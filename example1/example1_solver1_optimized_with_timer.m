%% Clean up
clear, clc, close all % Do NOT use "clear all"

tic()

N = 100; init_V = ones(N, 1);
[X, V, optimal_c, iter, err] = runSolver1_optimized(init_V);

toc()

%% Output results
fprintf("Number iterations to convergence: %i\n", iter)
fprintf("Final maximum error: %8.5g\n", err)

plot(X, optimal_c); xlabel("Cash-on-hand"); ylabel("Consumption")
title("Optimal consumption by cash-on-hand")
figure()
plot(X, V); xlabel("Cash-on-hand"); ylabel("Value function")
title("Value function by cash-on-hand")

%  Compute saving decision and plot it
optimal_s = X - optimal_c;
figure()
plot(X, optimal_s); xlabel("Cash-on-hand"); ylabel("Savings")
title("Optimal savings by cash-on-hand")
