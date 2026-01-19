# %% Prepare

using SymPy
using Plots

ENV["QT_QPA_PLATFORM"]="xcb"
gr()

const M_MU = 105.658

# %% Symbols

function λ(x, y, z)
  x^2 + y^2 + z^2 - 2 * x * y - 2 * y * z - 2 * z * x
end

@syms m_μ::(positive, real), m_N::(positive, real), p1::(positive, real), p3::(positive, real), G_F::(positive, real), E_N::(positive, real)

M2 = 64 * G_F^2 * m_μ * p3 * ((m_μ^2 - m_N^2) / 2 - p3 * m_μ)

LIPS = 1 / (8 * (2 * PI)^3 * m_μ) * p1 / sqrt(p1^2 + m_N^2) 

dΓ = M2 * LIPS |> simplify

ineq = [Eq(- p1 + p3, m_μ - sqrt(p1^2 + m_N^2) - p3), Eq(m_μ - sqrt(p1^2 + m_N^2) - p3, p1 + p3)]

p3max = solve(ineq[1], p3)[1]
p3min = solve(ineq[2], p3)[1]

p1max = solve(Eq(p1, m_μ - sqrt(p1^2 + m_N^2)), p1)[1]

dΓdp1 = integrate(dΓ, (p3, p3min, p3max)) |> simplify

# Γ_ν = integrate(dΓdE(m_N => 0) |> simplify, (p1, 0, p1max(m_N => 0))) |> simplify
Γ_N = integrate(dΓdp1 |> simplify, (p1, 0, p1max)) |> simplify |> expand |> factor
Γ_N = replace(Γ_N, sqrt(m_N^4 + 2 * m_N^2 * m_μ^2 + m_μ^4) => (m_N^2 + m_μ^2)) |> simplify |> factor
Γ_ν = limit(Γ_N, (m_N => 0))

BR = Γ_N / Γ_ν |> simplify

dΓdE = (dΓdp1 * E_N / p1)(p1 => sqrt(E_N^2 - m_N^2))

# %% Calculate BR

m_HNL = range(10, 30, length=100)
BRs = [BR(m_μ => M_MU, m_N => m) for m in m_HNL]
plot(m_HNL, BRs)

# %% Calculate Spectrum
E_range_test = range(0, M_MU, length=100)
spec_N, spec_ν = let tmp_dΓ_N = dΓdE(m_μ => M_MU, m_N => 20, G_F => 1) |> simplify, 
                  tmp_dΓ_ν = dΓdE(m_μ => M_MU, m_N => 0, G_F => 1) |> simplify
  [let a = N(limit(tmp_dΓ_N, E_N => E_r))
    imag(a) == 0 && a > 0 ? a : 0.0
  end
  for E_r in E_range_test],
  [let a = N(limit(tmp_dΓ_ν, E_N => E_r))
    imag(a) == 0 && a > 0 ? a : 0.0
  end
  for E_r in E_range_test]
end

plot(E_range_test, spec_ratios)