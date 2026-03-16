# %% prepare

using LsqFit
using QuadGK
using Plots
using FastBroadcast
using DelimitedFiles
using LinearAlgebra
using DataInterpolations
using Base.Threads
using ProgressMeter
using LaTeXStrings

const M_MU = 105.6583715e-3
const M_PI = 139.57018e-3

const R = (M_MU / M_PI)^2

function λ_func(x, y, z)
  return x^2 + y^2 + z^2 - 2 * x * y - 2 * y * z - 2 * z * x
end

# 0, (1 - r)
function kernel_func_μ1(x)
  return 1/(1-R)
end

# r, 1
function kernel_func_μg(x)
  return (3 - 2 * R) / (9 * (1 - R)^2) * (9 * x^2 - 6 * log(x) - 4 * x^3 - 5)
end

# 0, r
function kernel_func_μh1(x)
  return (3 - 2 * R) / (9 * (1 - R)^2) * (9 * R^2 - 6 * log(R) - 4 * R^3 - 5)
end

# 0, r
function kernel_func_μh2(x)
  return (1 + 2 * R) * (R - x) / (9 * R^2) * (9 * (R + x) - 4 * (R^2 + R*x + x^2))
end

# r, 1
function kernel_func_eg(x)
  return (2 * (1 - x)) / (3 * (1 - R)^2) * (6 * (1 - x)^2 + R * (5 + 5 * x - 4 * x^2) + 6 * R * log(x))
end

# 0, r
function kernel_func_eh1(x)
  return 2 / (3 * (1 - R)^2) * ((1 - R) * (6 - 7 * R + 11 * R^2 - 4 * R^3) + 6 * R * log(R))
end

# 0, r
function kernel_func_eh2(x)
  return 2 * (R - x) / (3 * R^2) * (7 * R^2 - 4 * R^3 + 7 * x * R - 4 * x * R^2 - 2 * x^2 - 4 * x^2 * R)
end

const INT_0 = 1e-8

function get_ν_spec(E_val, π_spec_func)
  # predictions = similar(E_ν_array)

  function gen_int(kfs, E_val)
    return (x) -> begin
      E_π   = E_val / x
      J_val = π_spec_func(E_π)
      sum([ 2 / 3 * kf(x) * J_val / x for kf in kfs])
    end
  end

  prediction = 0.0
  int_0λ = gen_int([kernel_func_μ1], E_val)
  int_r1 = gen_int([kernel_func_μg, kernel_func_eg], E_val)
  int_0r = gen_int([kernel_func_μh1, kernel_func_μh2, kernel_func_eh1, kernel_func_eh2], E_val)
  try
    prediction = quadgk(int_0λ, INT_0, 1 - R)[1] + quadgk(int_0r, INT_0, R)[1] + quadgk(int_r1, R, 1)[1]
  catch e
    prediction = 0.0
  end
  return prediction
end

function gen_get_γ_spec(m_N; E_range=(0, Inf), U = 1.0)
  γ_N = (M_PI^2 - M_MU^2 + m_N^2) / (2 * m_N * M_PI)

  function Q_γ_R(E_val)
    2 * E_val / (γ_N^2 * m_N^2)
  end
  function Q_γ_L(E_val)
    2 * (γ_N * m_N - E_val) / (γ_N^2 * m_N^2)
  end

  function λ_func(x, y, z)
    return x^2 + y^2 + z^2 - 2 * x * y - 2 * y * z - 2 * z * x
  end
  temp_w = (M_MU^2 + m_N^2) * M_PI^2 - (M_MU^2 - m_N^2)^2
  w_R = temp_w - (M_MU^2 - m_N^2) * sqrt(λ_func(M_PI^2, M_MU^2, m_N^2))
  w_L = temp_w + (M_MU^2 - m_N^2) * sqrt(λ_func(M_PI^2, M_MU^2, m_N^2))
  Γ_N = (- (M_MU^2 - m_N^2)^2 + M_PI^2 * (M_MU^2 + m_N^2)) * sqrt(λ_func(M_PI^2, M_MU^2, m_N^2)) * U^2
  Γ_ν = M_MU^2 * (M_PI^2 - M_MU^2)^2
  BR_N = Γ_N / (Γ_N + Γ_ν)
  # println(w_R)
  # println(w_L)

  function Q_γ_mix(E_val)
    return (Q_γ_L(E_val) * w_L + Q_γ_R(E_val) * w_R) / (w_L + w_R)
  end

  function get_γ_spec_raw(E_val, π_spec_func)
    function gen_int_E_π(E_γ_in_π)
      return function int_E_π(E_π)
        π_spec_func(E_π) / (2 * E_π / M_PI * E_γ_in_π)
      end
    end

    prediction = 0.0
    prediction = 
      quadgk((E_γ_in_π) -> let low_b = max(M_PI * E_val / 2 / E_γ_in_π, E_range[1])
      if low_b > E_range[2]
        0.0
      else
        quadgk(gen_int_E_π(E_γ_in_π), low_b, E_range[2])[1] * Q_γ_mix(E_γ_in_π)
      end
      end,
      INT_0, γ_N * m_N)[1] * BR_N
  end
  return get_γ_spec_raw
end

function get_γ_spec_approx(HNL_spec; E_range=(0, Inf))
  function γ_spec(E)
    low_b = max(E, E_range[1])
    prediction = quadgk((E_N) -> HNL_spec(E_N) * 2 * (E_N - E) / E_N^2, low_b, E_range[2])[1]
    return prediction
  end

  return γ_spec
end

const s2w = 0.2397
function get_e_spec_approx(HNL_spec; E_range=(0, Inf))
  BR = (1 - 4 * s2w + 8 * s2w^2)/4
  function γ_spec(E)
    low_b = max(E, E_range[1])
    prediction = quadgk((E_N) -> HNL_spec(E_N) / (3 * E_N^4) * (11 * E_N^3 - 27 * E_N * E^2 + 16 * E^3) * BR, low_b, E_range[2])[1]
    return prediction
  end

  return γ_spec
end

function get_γ_from_e(e_spec)
  ϵ_0 = 6.3e-13
  m_e = 0.000510998951
  function γ_spec(E)
    E_e = m_e * sqrt(3 * E / 4 / ϵ_0)
    return e_spec(E_e) * E_e / E / 2
  end
end

# %% read file

lit_ν_spec = readdlm("datas/20220913_Evidence_for_neutrino_emission_from_the_nearby_active_galaxy_NGC_1068_data/resources/Fig4_SED/model_murase_et_al.txt", skipstart=1)

icecube_spec = readdlm("datas/20220913_Evidence_for_neutrino_emission_from_the_nearby_active_galaxy_NGC_1068_data/resources/Fig4_SED/ngc1068_spectrum_95.txt", skipstart=1)

low_γ_obs = readdlm("datas/20220913_Evidence_for_neutrino_emission_from_the_nearby_active_galaxy_NGC_1068_data/resources/Fig4_SED/gammaray_0.1_to_100_GeV.txt", skipstart=1)

high_γ_obs = readdlm("datas/20220913_Evidence_for_neutrino_emission_from_the_nearby_active_galaxy_NGC_1068_data/resources/Fig4_SED/gammaray_above_200_GeV.txt", skipstart=1)


# p_π = [1e-16, 0.48, 1000.0, 1500.0]

# @. fit_model(E, p) = get_ν_spec(E, (E_π) -> p[1] * (E_π)^(-p[2]) * exp(- E_π/p[3])) * E^2
# fit = curve_fit(fit_model, lit_ν_spec[:,1],  lit_ν_spec[:, 2] .* 10^-3, [1e-16, 0.48, 1500.0])
# p_π = [5e-17, 1000, 1/2, 2000.0, 1000.0]

# function π_spec_func(E)
#   p_π[1] * (E / p_π[2])^(- p_π[3]) * exp(- E / p_π[4]) 
# end

ν_spec_func = AkimaInterpolation(lit_ν_spec[:, 2] .* 1e-3 ./ (lit_ν_spec[:, 1] .^2), lit_ν_spec[:, 1])

E_ν_array = 10 .^ range(log10(minimum(lit_ν_spec[:,1]) + 1), log10(maximum(lit_ν_spec[:,1]) - 10), length=50)

const E_interp_min = minimum(lit_ν_spec[:, 1])
const E_interp_max = maximum(lit_ν_spec[:, 1])

# %% calculate

function π_spec_func(E)
  E_scaled = E * 0.4
  if E_scaled < E_interp_min || E_scaled > E_interp_max
    return 0.0
  end
  return 0.4^2 * ν_spec_func(E_scaled)
end

Q_predict = zeros(length(E_ν_array))
@showprogress @threads for i in eachindex(E_ν_array)
  E = E_ν_array[i]
  Q_predict[i] = get_ν_spec(E, π_spec_func)
end

Q_γ_predict_1, Q_func_1 = zeros(length(E_ν_array)),  gen_get_γ_spec(0.1e-3; E_range=(E_interp_min, E_interp_max))
@showprogress @threads for i in eachindex(E_ν_array)
  E = E_ν_array[i]
  Q_γ_predict_1[i] = Q_func_1(E, π_spec_func)
end

Q_γ_predict_2, Q_func_2 = zeros(length(E_ν_array)),  gen_get_γ_spec(30e-3; E_range=(E_interp_min, E_interp_max))
@showprogress @threads for i in eachindex(E_ν_array)
  E = E_ν_array[i]
  Q_γ_predict_2[i] = Q_func_2(E, π_spec_func)
end
# = gen_get_γ_spec(0.1)(E_ν_array, Ref(π_spec_func))

# Q_γ_predict_25 = @.. thread=true gen_get_γ_spec(25.0)(E_ν_array, Ref(π_spec_func))


plot(E_ν_array, Q_predict .* (E_ν_array .^ 2), xscale = :log10, yscale = :log10, ylims=(1e-14, 1e-9), xlims=(1,1e7), label="Neutrino confirm", color = :gray)
plot!(E_ν_array, π_spec_func.(E_ν_array) .* (E_ν_array .^ 2), label = "Pion", color = :black)
plot!(E_ν_array,  Q_γ_predict_1 .* (E_ν_array .^ 2), label="m_N = 0.1 MeV")
plot!(E_ν_array,  Q_γ_predict_2 .* (E_ν_array .^ 2), label="m_N = 30 MeV")
plot!(lit_ν_spec[:,1], lit_ν_spec[:, 2] .* 1e-3, label="Neutrino original", color = :red)
plot!(icecube_spec[:, 1], icecube_spec[:, 2], label="Neutrino IceCube", color = :blue)

# %% approx γ

HNL_spec = (E) -> 
  if E < E_interp_min || E > E_interp_max
    0.0
  else
    ν_spec_func(E) * 2 * 10^-1.5
  end

γ_spec = get_γ_spec_approx(HNL_spec; E_range=(E_interp_min, E_interp_max))

pgfplotsx()

p2 = plot(lit_ν_spec[:,1], lit_ν_spec[:, 2] .* 1e-3, label="Neutrino (Theory)", color = :black, xscale = :log10, yscale = :log10, ylims=(1e-14, 1e-9), xlims=(1e-1,1e5), linestyle = :dot, framestyle = :box, size=(400,300), xlabel=L"E \; \mathrm{(GeV)}", ylabel=L"E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})", legend=:topleft)
plot!(p2, icecube_spec[:, 1], icecube_spec[:, 2], label="Neutrino (IceCube)", color = :blue)
plot!(p2, (E) -> γ_spec(E) * E^2, label="HNL γ", color = :red)
scatter!(p2, low_γ_obs[:,1], low_γ_obs[:,2],
        label="γ-ray 0.1 to 100 GeV", color=:green)
scatter!(p2, high_γ_obs[:,1], high_γ_obs[:,2],
        label="γ-ray > 200 GeV")

# savefig(p2, "julia_plots/u1.5_gamma.pdf")

# %% constrain U

obs_γ_UL = [low_γ_obs[:, 1:2]; high_γ_obs[:,1:2]]

function get_U_UL(m_N)
  ρ_π =  (M_PI^2 * (M_MU^2 + m_N^2) - (M_MU^2 - m_N^2)^2) * sqrt(λ_func(M_PI^2, M_MU^2, m_N^2)) / (M_MU^2 * (M_PI^2 - M_MU^2)^2)
  r_N = m_N / M_MU
  ρ_μ = -r_N^8 + 8*r_N^6 - 24*r_N^4*log(r_N) - 8*r_N^2 + 1

  println("RHO PI: $ρ_π\nRHO MU: $ρ_μ")

  HNL_spec_approx = (E) ->
  if E < E_interp_min || E > E_interp_max
    0.0
  else
    ν_spec_func(E) * (ρ_π + ρ_μ)
  end
  γ_spec_approx = get_γ_spec_approx(HNL_spec_approx; E_range=(E_interp_min, E_interp_max))

  predic_Φ = [ γ_spec_approx(E) * E^2 for E in obs_γ_UL[:,1]]
  U_UL = minimum([obs_γ_UL[i, 2] / predic_Φ[i] for i in eachindex(predic_Φ) if predic_Φ[i] > 0])  

  return U_UL

end

m_N_range = range(1e-3, 3e-2; length=100)
U_UL_list = zeros(length(m_N_range))
@showprogress @threads for i in eachindex(U_UL_list)
  U_UL_list[i] = get_U_UL(m_N_range[i])
end
p3 = plot(m_N_range .* 1e3, U_UL_list, color = :black, framestyle = :box, size=(400,300), label=nothing, xlabel=L"m_N \; \mathrm{(MeV)}", ylabel=L"|U_\mu|^2")
# savefig(p3, "julia_plots/U_UL.pdf")


# %% approx e

HNL_spec_e = (E) -> 
  if E < E_interp_min || E > E_interp_max
    0.0
  else
    ν_spec_func(E) * 2 * 10^-1
  end

e_spec = get_e_spec_approx(HNL_spec_e; E_range=(E_interp_min, E_interp_max))

γ_N_spec = get_γ_from_e(e_spec)

pgfplotsx()

p4 = plot(lit_ν_spec[:,1], lit_ν_spec[:, 2] .* 1e-3, label="Neutrino (Theory)", color = :black, xscale = :log10, yscale = :log10, ylims=(1e-16, 1e-9), xlims=(1e-1,1e5), linestyle = :dot, framestyle = :box, size=(500,400), xlabel=L"E \; \mathrm{(GeV)}", ylabel=L"E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})", legend=:topleft)
plot!(p4, icecube_spec[:, 1], icecube_spec[:, 2], label="Neutrino (IceCube)", color = :blue)
plot!(p4, (E) -> e_spec(E) * E^2, label="HNL e", color = :red)
plot!(p4, (E) -> γ_N_spec(E) * E^2, label="HNL γ", color = :green)
scatter!(p4, low_γ_obs[:,1], low_γ_obs[:,2],
        label="γ-ray 0.1 to 100 GeV", color=:green)
scatter!(p4, high_γ_obs[:,1], high_γ_obs[:,2],
        label="γ-ray > 200 GeV")

# savefig(p2, "julia_plots/u1.5_gamma.pdf")