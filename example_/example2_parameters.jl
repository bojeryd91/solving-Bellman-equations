const β = 0.95
const ϕ = 0.8
const c_min = 0.0001
const a = 1.0
const s_min = -2.0
u(c, h) = c >= c_min && h >= 0.0 ? log(c) - a*h^2 : -1.0e8

const Pη = ones(9)./9.0
const ηη = range(-0.05, 0.05, length=length(Pη))

const N_S = 80; const N_ϵ = 40
const S   = range(s_min,  4.0, length=N_S)
const ϵϵ  = range( 0.75, 1.25, length=N_ϵ)

function Expect(tilde_s, ϵ, V_interp)
    possible_ϵ_tp1 = ϵ^ϕ*exp.(ηη)
    sum(V_interp.((tilde_s, ), possible_ϵ_tp1).*Pη)
end
