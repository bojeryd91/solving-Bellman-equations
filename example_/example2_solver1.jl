## Import packages
using Plots, Interpolations, Optim

## Import constants and functions
include("example2_parameters.jl")

## The solver
function solver21(init_V, r, w)
    curr_V    = init_V; new_V = copy(curr_V)
    optimal_c = c_min.*ones(N_S, N_ϵ)
    optimal_h = 2.0.*ones(N_S, N_ϵ)

    iter = 1; iter_max = 500; err = 1.0e10; err_tol = 1.0e-4
    while iter < iter_max && err > err_tol
        print("Iteration $iter, error: $err");
        # Create an interpolation object V_interp
        V_interp = LinearInterpolation((S, ϵϵ), curr_V, extrapolation_bc=
                                                       Interpolations.Flat())
        for (i_tilde_s, tilde_s) in enumerate(S)
            for (i_ϵ, ϵ) in enumerate(ϵϵ)
                # Create the function to minimize
                R = function(c, h)
                    x = tilde_s + w*ϵ*h
                    s = x - c
                    if s < s_min
                        return -1.0e8
                    end
                    tilde_s_tp1 = s*(1.0+r)
                    return u(c, h) + β*Expect(tilde_s_tp1, ϵ, V_interp)
                end

                # For this (tilde_s, ϵ), maximize R(c, h)
                result  = optimize(arg -> -R(arg[1], arg[2]),
                                    [optimal_c[i_tilde_s, i_ϵ],
                                     optimal_h[i_tilde_s, i_ϵ]])
                result2 = optimize(arg -> -R(arg[1], arg[2]),
                                    [c_min+0.0001, 2.0])

                if -result.minimum < -result2.minimum
                    result = result2
                end

                # Identify the maximum and maximizer
                new_V[i_tilde_s, i_ϵ] = -result.minimum
                optimal_c[i_tilde_s, i_ϵ] = result.minimizer[1]
                optimal_h[i_tilde_s, i_ϵ] = result.minimizer[2]
            end
        end
        # Prepare for next iteration
        err  = maximum(abs.(new_V .- curr_V))
        curr_V = copy(new_V)
        iter += 1
        print("\e[2K\e[1G") # clear whole line and move cursor to column 1
    end
    print("Done!")
    return (new_V, optimal_c, optimal_h, err, iter)
end
## Run solver
(V, optimal_c, optimal_h, err, iter) = solver21(ones(N_S, N_ϵ), 0.08, 0.5)

## Output results
print(""); println("Number iterations to convergence: $iter")
println("Final maximum error: $err")

display(plot(S, V[:, 1], label=""))
display(plot(S, optimal_c[:, 1], label=""))
display(plot(S, optimal_h[:, 1], label=""))
