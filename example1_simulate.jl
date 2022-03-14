## Solve model using one of the solvers and add decision rule to the workspace.
include("example1_solver3.jl")

## Create simulation command
function simulate_household(init_x, num_periods)
    c_dec = LinearInterpolation(X, optimal_c)
    Y = Array{Float64}(undef, num_periods, 3)
    Y[1, 1] = init_x
    # Start simulation
    for t in 1:num_periods
        x_t = Y[t, 1]
        # Compute today's optimal consumption
        Y[t, 2] = c_dec(x_t)
        # Compute saving decision
        Y[t, 3] = x_t - Y[t, 2]
        # Compute tomorrow's cash-on-hand
        if t < num_periods
            Y[t+1, 1] = Y[t, 3]*(1.0+r) + y
        end
    end
    return Y
end

## Simulate
sim1  = simulate_household(2.0, 10)
plot1 = plot(1:size(sim1)[1], sim1[:, 1], label="", title="Simulation")
display(plot1); savefig(plot1, "./figures/simulation_11")

sim2  = simulate_household(1.5, 20)
plot2 = plot(1:size(sim2)[1], sim2[:, 1], label="", title="Simulation")
sim_ss = simulate_household(1.5, 1000)
plot!([1, size(sim2)[1]], [sim_ss[end, 1], sim_ss[end, 1]], label="",
        linestyle=:dash)
display(plot2); savefig(plot2, "./figures/simulation_12")
