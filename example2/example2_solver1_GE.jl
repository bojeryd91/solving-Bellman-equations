## Import packages
using Plots, Interpolations, Optim, LinearAlgebra

## Import constants and functions
include("example2_parameters.jl")
include("solver21.jl")

function solveGE_solver21(r_0, w_0)
    V_0 = ones(N_S, N_ϵ); V_guess = V_0
    r_guess = r_0; w_guess = w_0
    iter = 1; iter_max = 500; err = 1.0e10; err_tol = 1.0e-6
    while iter < iter_max && err > err_tol
        print("\e[2K\e[1G") # clear whole line and move cursor to column 1
        print("Iteration $iter, error: $err")

        # For this guess of r and w, solve the household problem
        (V_guess, optimal_c, _, _) = solver21(V_guess, r_guess, w_guess,
                                                print_output=false)
        # Compute the optimal saving decision HERE FOR NOW, MOVE INTO solver21
        optimal_s = copy(optimal_c)
        for (i_s, s) in enumerate(S)
            for (i_ϵ, ϵ) in enumerate(ϵϵ)
                optimal_s[i_s, i_ϵ] = s + exp(ϵ)*y_bar*w_guess -
                                                    optimal_c[i_s, i_ϵ]
            end
        end

        # Construct the transition matrix G_sϵ
        G_sϵ = zeros(N_S*N_ϵ, N_S*N_ϵ)
        for (i_s, s) in enumerate(S)
            for (i_ϵ, ϵ) in enumerate(ϵϵ)
                optimal_s_t = optimal_s[i_s, i_ϵ]
                s_tp1 = optimal_s_t*(1.0+r_guess)
                i_s_tp1 = findlast(S .<= s_tp1)
                i_today_G    = (i_s-1)*N_ϵ + i_ϵ
                i_tomorrow_G = (i_s_tp1-1)*N_ϵ .+ (1:N_ϵ)
                G_sϵ[i_tomorrow_G, i_today_G] .=
                                G_sϵ[i_tomorrow_G, i_today_G] .+ Pϵ
            end
        end
        i_eig = findall(real.(eigvals(G_sϵ)).≈1.0)[1]
        distr_ss  = real.(eigvecs(G_sϵ)[:, i_eig])
        distr_ss .= sign(distr_ss[1]).*distr_ss
        if minimum(distr_ss) < 0.0
            println("problem")
        end
        distr_ss .= distr_ss./sum(distr_ss)
        fig =
        plot(S, distr_ss[1:N_ϵ:end], label="")
        display(fig)
        # sleep(4.0)

        # Aggregate savings across households
        Sˢ = sum(sum.([S.*distr_ss[i_ϵ:N_ϵ:end] for i_ϵ in eachindex(ϵϵ)]))

        # Compute demand for labor
        Lᶠ = (α*A/w_guess)^(1.0-α)
        # Make new guess for w and r
        w_guess = w_guess*(1.0 + log(Lᶠ/1.0)/2.0)
        r_guess = r_guess*(1.0-log(abs(Sˢ))/10.0)

        print(", r=$r_guess, w=$w_guess")
        sleep(4.0)
        iter += 1
        err   = 0.1
    end
end
solveGE_solver21(0.05, 0.5)
