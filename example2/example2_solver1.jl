## Import packages
using Plots, Interpolations, Optim

## Import constants and functions
include("example2_parameters.jl")
include("solver21.jl")

## Run solver
(V, optimal_c, err, iter) = solver21(ones(N_S, N_ϵ), 0.05, 0.5)

## Output results
print(""); println("Number iterations to convergence: $iter")
println("Final maximum error: $err")

V_plot =
plot(  S, V[:, 1], label="ϵ_1", legend=:topleft,
        ylabel="Value function", xlabel="s")
[plot!(S, V[:, i_ϵ], label="ϵ_$i_ϵ") for i_ϵ in 2:N_ϵ]
display(V_plot)
display(plot(S, optimal_c[:, 1], label="", ylabel="optimal c", xlabel="s"))
