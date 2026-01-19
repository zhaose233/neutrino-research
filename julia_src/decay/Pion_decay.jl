# %% Prepare

using SymPy
using Plots

ENV["QT_QPA_PLATFORM"]="xcb"
gr()

# %% Symbols

function λ(x, y, z)
  x^2 + y^2 + z^2 - 2 * x * y - 2 * y * z - 2 * z * x
end

@syms m_π::(positive, real) m_l::(positive, real) m_N::(positive, real) p2::(positive, real)

eq_p = Eq(m_π , sqrt(m_l^2 + p2) + sqrt(m_N^2 + p2))

sol_p2 = solve(eq_p, p2)[1]

E_l = (m_π^2 - m_N^2 + m_l^2) / (2 * m_π)
E_N = (m_π^2 + m_N^2 - m_l^2) / (2 * m_π)

M2 = 2 * m_π^2 * E_l * E_N - m_π^2 * (E_l * E_N + p2)

LIPS  = sqrt(p2) / m_π

Γ = M2 * LIPS

Γ_pure = Γ(p2 => sol_p2) |> simplify

BR_eμ = Γ_pure(m_π => 139.0, m_l => 0.5, m_N => 0) / Γ_pure(m_π => 139.0, m_l => 105.0, m_N => 0)
println("e/μ Brand Ratio = $BR_eμ")

# %% BR HNL

m_HNL = range(10, 30, length = 1000)

BR_HNL = [Γ_pure(m_π => 139.0, m_l => 105.0, m_N => m) for m in m_HNL] ./ Γ_pure(m_π => 139.0, m_l => 105.0, m_N => 0) 

p1 = plot(m_HNL, BR_HNL, title = "π^± decay branch ratio for HNL/ν_μ ~ HNL mass", label = nothing)

mkpath("julia_results")
savefig(p1, "julia_results/HNL_BR.pdf")

# %% Energy distribution

@syms θ::(positive, real), β::(positive, real), E_π::(positive, real)

γ = 1/sqrt(1 - β^2)

sol_β = solve(Eq(E_π, m_π * γ), β)[1]


E_N_max = E_N * γ + sqrt(sol_p2) * γ * β

E_N_max(β => sol_β) |> simplify |> factor |> simplify

E_ratio = (m_N^2 - m_l^2 + m_π^2 + sqrt(λ(m_N^2, m_l^2, m_π^2))) / (2 * m_π^2)

ratio_l = [E_ratio(m_N => m, m_l => 105.0, m_π => 139.0) for m in m_HNL]

p1 = plot(m_HNL, ratio_l, title = "HNL max energy ~ mass", label = nothing)
