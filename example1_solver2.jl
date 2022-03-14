## Import packages
using Plots, Interpolations

## Define economic parameters
const β = 0.95
const r = 1.0/β - 1.0
const y = 1.0
const c_min = 0.0001
u(c) = c >= c_min ? log(c) : -1.0e8

## Define solver parameters
const N_X = 100
const N_C = 800

## Initialize vectors
const X = range(1.0, 4.0, length=N_X)
const C = range(0.0, maximum(X), length=N_C)
curr_V  = ones(N_X); new_V = copy(curr_V)
optimal_c = zeros(N_X) # Why N_X?

## The solver
iter = 0; iter_max = 1000; err = 1.0e10; err_tol = 1.0e-4
while iter < iter_max && err > err_tol
    # Create an interpolation object V_interp
    V_temp = LinearInterpolation(X, curr_V, extrapolation_bc=Flat())
    V_interp(x_tp1) = x_tp1 >= y ? V_temp(x_tp1) : -1.0e8
    for (i_x, x) in enumerate(X)
        # For a particular x, compute the right-hand-side of the
        #   Bellman equation for each possible choice c ∈ C.
        Bellman_rhs = u.(C) .+ β.*V_interp.((x .- C).*(1.0+r) .+ y)

        # Identify the maximum and maximizer
        (maxim, index_maximizer) = findmax(Bellman_rhs)
        new_V[i_x] = maxim
        global optimal_c[i_x] = C[index_maximizer]
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
