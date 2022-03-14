## Import packages
using Plots, Interpolations, Optim

## Define economic parameters
const β = 0.95
const r = 1.0/β - 1.0 - 0.01
const y = 1.0
const c_min = 0.0001
u(c) = c >= c_min ? log(c) : -1.0e8

## Define solver parameters and initialize vectors
const N_X = 100
const X   = range(1.0, 4.0, length=N_X)
const N_C = 800
const C = range(0.0, maximum(X), length=N_C)

## The solvers
function solver11(V_guess)
    # Create variables as in unwrapped code before
    curr_V = V_guess; new_V = copy(curr_V)
    optimal_c = copy(curr_V)

    # Below is exactly as in Solver 3 except that globals are removed
    iter = 0; iter_max = 1000; err = 1.0e10; err_tol = 1.0e-4
    while iter < iter_max && err > err_tol
        for (i_x, x) in enumerate(X)
            # For a particular x, compute the right-hand-side of the
            #   Bellman equation for each possible choice x′ ∈ X.
            Bellman_rhs = u.(x .- (X .- y)./(1.0+r)) .+ β.*curr_V

            # Identify the maximum and maximizer
            (maxim, index_maximizer) = findmax(Bellman_rhs)
            new_V[i_x]  = maxim
            optimal_x′ = X[index_maximizer]
            optimal_c[i_x] = x - (optimal_x′ - y)/(1+r)
        end
        # Prepare for next iteration
        err    = maximum(abs.(new_V .- curr_V))
        curr_V = copy(new_V)
        iter += 1
    end
    return (curr_V, optimal_c, iter, err)
end

function solver12(V_guess)
    # Create variables as in unwrapped code before
    curr_V = V_guess; new_V = copy(curr_V)
    optimal_c = copy(curr_V)

    # Below is exactly as in Solver 3 except that globals are removed,
    #   and the addition of Interpolations. before Flat()
    iter = 0; iter_max = 1000; err = 1.0e10; err_tol = 1.0e-4
    while iter < iter_max && err > err_tol
        # Create an interpolation object V_interp
        V_temp = LinearInterpolation(X, curr_V, extrapolation_bc=
                                                    Interpolations.Flat())
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
        err  = maximum(abs.(new_V .- curr_V))
        curr_V = copy(new_V)
        iter += 1
    end
    return (curr_V, optimal_c, iter, err)
end

function solver13(V_guess)
    # Create variables as in unwrapped code before
    curr_V = V_guess; new_V = copy(curr_V); optimal_c = copy(curr_V)

    # Below is exactly as in Solver 3 except that globals are removed
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
            optimal_c[i_x] = result.minimizer
        end
        # Prepare for next iteration
        err  = maximum(abs.(new_V .- curr_V))
        curr_V = copy(new_V)
        iter += 1
    end
    return (curr_V, optimal_c, iter, err)
end

## Run all solvers and save results
first_guess_V = ones(N_X)
(V_sol11, optimal_c_sol11, _, _) = solver11(first_guess_V)
(V_sol12, optimal_c_sol12, _, _) = solver12(first_guess_V)
(V_sol13, optimal_c_sol13, _, _) = solver13(first_guess_V)

## Compare computation time, iterations, and final error of each solver
iter = 1; print("")
for solver in (solver11, solver12, solver13)
    println("Using solver$iter:")
    (_, _, num_iter, final_err) = solver(first_guess_V)
    print(""); @time solver(first_guess_V)

    print(""); println("Number iterations to convergence: $num_iter")
    println("Final maximum error: $final_err")
    global iter += 1
end

## Output results
plot_c =
plot( X, optimal_c_sol11, label="Solver 1")
plot!(X, optimal_c_sol12, label="Solver 2",
                        color=:red,   linestyle=:dash)
plot!(X, optimal_c_sol13, label="Solver 3",
                        color=:green, linestyle=:dash)
plot!(xlabel="Cash-on-hand", ylabel="Consumption",
      title="Optimal consumption by cash-on-hand", legend=:bottomright)
display(plot_c)

plot_V =
plot( X, V_sol11, label="Solver 1")
plot!(X, V_sol12, label="Solver 2", color=:red,   linestyle=:dash)
plot!(X, V_sol13, label="Solver 3", color=:green, linestyle=:dash)
plot!(xlabel="Cash-on-hand", ylabel="Value function",
    title="Value function by cash-on-hand", legend=:bottomright)
display(plot_V)

#  Compute saving decision and plot it
optimal_s_sol_11 = X .- optimal_c_sol11
optimal_s_sol_12 = X .- optimal_c_sol12
optimal_s_sol_13 = X .- optimal_c_sol13
plot_s =
plot( X, optimal_s_sol_11, label="Solver 1")
plot!(X, optimal_s_sol_12, label="Solver 2",
                            color=:red,   linestyle=:dash)
plot!(X, optimal_s_sol_13, label="Solver 3",
                            color=:green, linestyle=:dash)
plot!(xlabel="Cash-on-hand", ylabel="Savings",
    title="Optimal savings by cash-on-hand", legend=:bottomright)
display(plot_s)
