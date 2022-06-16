## The solver
function solver21(init_V, r, w; print_output=true)
    curr_V    = init_V; new_V = copy(curr_V)
    optimal_c = c_min.*ones(N_S, N_ϵ)

    iter = 1; iter_max = 500; err = 1.0e10; err_tol = 1.0e-6
    while iter < iter_max && err > err_tol
        if mod(iter, 10) == 1 && print_output == true
            print("\e[2K\e[1G") # clear whole line and move cursor to column 1
            print("Iteration $iter, error: $err")
        end
        # Create an interpolation object V_interp
        V_temp = LinearInterpolation((S, ϵϵ), curr_V, extrapolation_bc=
                                                       Interpolations.Flat())
        V_interp(s, ϵ) = s < s_min ? -1.0e8 : V_temp(s, ϵ)

        for (i_tilde_s, tilde_s) in enumerate(S)
            for (i_ϵ, ϵ) in enumerate(ϵϵ)
                x = tilde_s + exp(ϵ)*y_bar*w

                # Create the function to minimize
                R = function(c)
                    s = x - c
                    if s < s_min
                        return -1.0e8
                    end
                    tilde_s_tp1 = s*(1.0+r)
                    return u(c) + β*Expect(tilde_s_tp1, V_interp)
                end

                # For this (tilde_s, ϵ), maximize R(c)
                result = optimize(c -> -R(c), c_min, x - s_min)

                # Identify the maximum and maximizer
                    new_V[i_tilde_s, i_ϵ] = -result.minimum
                optimal_c[i_tilde_s, i_ϵ] =  result.minimizer
            end
        end
        # Prepare for next iteration
        err    = maximum(abs.(new_V .- curr_V))
        curr_V = copy(new_V)
        iter  += 1
    end
    # Print last iteration
    if print_output == true
        print("\e[2K\e[1G") # clear whole line and move cursor to column 1
        println("Iteration $iter, error: $err")
        println("Done!")
    end
    if iter == iter_max
        return # Return nothing if iter_max was reached
    else
        return (new_V, optimal_c, err, iter)
    end
end
