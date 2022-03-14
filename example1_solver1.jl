## Import packages
using Plots

## Define economic parameters
const β = 0.95
const r = 1.0/β - 1.0
const y = 1.0
const c_min = 0.0001
u(c) = c >= c_min ? log(c) : -1.0e8

## Define solver parameters
const N = 100 # Note that I didn't write 100.0, so this will be a Int64

## Initialize vectors
const X = range(1.0, 4.0, length=N)
curr_V  = ones(N)
new_V   = copy(curr_V) # Use copy(), not new_V = curr_V
optimal_c = zeros(N)

## The solver
iter = 0; iter_max = 1000; err = 1.0e10; err_tol = 1.0e-4
while iter < iter_max && err > err_tol
    for (i_x, x) in enumerate(X)
        # For a particular x, compute the right-hand-side of the
        #   Bellman equation for each possible choice x′ ∈ X.
        Bellman_rhs = u.(x .- (X .- y)./(1.0+r)) .+ β.*curr_V

        # Identify the maximum and maximizer
        (maxim, index_maximizer) = findmax(Bellman_rhs)
        new_V[i_x]  = maxim
            optimal_x′ = X[index_maximizer] # Optimal decision for tomorrow's
            optimal_c[i_x] = x - (optimal_x′ - y)/(1+r) #        cash-on-hand
    end
    # Prepare for next iteration
    global err    = maximum(abs.(new_V .- curr_V))
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
