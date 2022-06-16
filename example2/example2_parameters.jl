const β = 0.95
const A = 1.0
const α = 2/3
const c_min = 0.0001
const y_bar = 1.0
const s_min = -y_bar
u(c) = c >= c_min ? log(c) : -1.0e8

const ϵϵ  = [-0.3, -0.1, 0.0, 0.1, 0.3]
const N_ϵ = length(ϵϵ)
const Pϵ  = ones(N_ϵ)./N_ϵ

const N_S = 300
const S   = reduce(vcat, [s_min, s_min+0.001, s_min+0.01,
                            range(s_min+0.02, 2.0, length=N_S-3)])

function Expect(tilde_s, V_interp)
    sum(V_interp.((tilde_s, ), ϵϵ).*Pϵ)
end
