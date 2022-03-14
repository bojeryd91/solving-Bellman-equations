## Import packages
using Plots, Interpolations, Optim

## Define economic parameters
const β = 0.95
const r = 1.0/β - 1.0
const y = 1.0
const c_min = 0.0001
u(c) = c >= c_min ? log(c) : -1.0e8

## Define solver parameters and initialize vectors
const N_X = 100
const X   = range(1.0, 4.0, length=N_X)
curr_V    =  ones(N_X); new_V = copy(curr_V)
optimal_c = zeros(N_X)

## The solver
iter = 0; iter_max = 1000; err = 1.0e10; err_tol = 1.0e-4
while iter < iter_max && err > err_tol
    # Create an interpolation object V_interp
    V_temp = LinearInterpolation(X, curr_V, extrapolation_bc=
                                                Interpolations.Flat())
    V_interp(x_tp1) = x_tp1 >= y ? V_temp(x_tp1) : -1.0e8
    for (i_x, x) in enumerate(X)
        # Create the function to minimize
        Bellman_rhs(c) = u(c) + β*V_interp((x - c)*(1.0+r) + y)

        # For a particular x, optimize the whole right-hand-side expression
        #   of the Bellman equation w.r.t. c
        result = optimize(c -> -Bellman_rhs(c), c_min, x)

        # Identify the maximum and maximizer
        new_V[i_x] = -result.minimum
        global optimal_c[i_x] = result.minimizer
    end
    # Prepare for next iteration
    global err  = maximum(abs.(new_V .- curr_V))
    global curr_V = copy(new_V)
    global iter += 1
end

## Output results
print(""); println("Number iterations to convergence: $iter")
println("Final maximum error: $err")

plots = [
plot(X, optimal_c, label="", xlabel="Cash-on-hand", ylabel="Consumption",
        title="Optimal consumption by cash-on-hand")]
push!(plots,
plot(X, curr_V,    label="", xlabel="Cash-on-hand", ylabel="Value function",
        title="Value function by cash-on-hand"))

#  Compute saving decision and plot it
optimal_s = X .- optimal_c
push!(plots,
plot(X, optimal_s, label="", xlabel="Cash-on-hand", ylabel="Savings",
        title="Optimal savings by cash-on-hand"))

[display(plot) for plot in plots]
