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
using PGFPlotsX
using ColorSchemes
using DataFrames
using CSV


const M_MU = 105.6583715
const M_PI = 139.57018

const R = (M_MU / M_PI)^2

function λ_func(x, y, z)
    return x^2 + y^2 + z^2 - 2 * x * y - 2 * y * z - 2 * z * x
end

function gen_get_γ_spec(m_N; E_range=(0, Inf), U=1.0)
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
        # prediction = quadgk((E_N) -> HNL_spec(E_N) / E_N, low_b, E_range[2])[1]
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

function get_ρ(m_N)
    ρ_π = (M_PI^2 * (M_MU^2 + m_N^2) - (M_MU^2 - m_N^2)^2) * sqrt(λ_func(M_PI^2, M_MU^2, m_N^2)) / (M_MU^2 * (M_PI^2 - M_MU^2)^2)
    r_N = m_N / M_MU
    ρ_μ = -r_N^8 + 8*r_N^6 - 24*r_N^4*log(r_N) - 8*r_N^2 + 1
    return ρ_π + ρ_μ
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

ν_spec_func = AkimaInterpolation(lit_ν_spec[:, 2] .* 1e-3 ./ (lit_ν_spec[:, 1] .^ 2), lit_ν_spec[:, 1])

E_ν_array = 10 .^ range(log10(minimum(lit_ν_spec[:, 1]) + 1), log10(maximum(lit_ν_spec[:, 1]) - 10), length=50)

const E_interp_min = minimum(lit_ν_spec[:, 1])
const E_interp_max = maximum(lit_ν_spec[:, 1])

# %% calculate

function gen_π_spec_func(A, Γ, E_cut, β)
    function π_spec(E)
        return A * E^(-Γ) * exp(-(E / E_cut)^β)
    end
    return π_spec
end

function get_ν_spec(π_spec)
    r = 0.573
    λ = 0.427

    # Eq (37)
    function g_ν_μ(x)
        return (3 - 2*r) / (9 * (1 - r)^2) * (9*x^2 - 6*log(x) - 4*x^3 - 5)
    end

    # Eq (38)
    function h_ν_μ_1(x)
        return (3 - 2*r) / (9 * (1 - r)^2) * (9*r^2 - 6*log(r) - 4*r^3 - 5)
    end

    # Eq (39)
    function h_ν_μ_2(x)
        return ((1 + 2*r) * (r - x)) / (9 * r^2) * (9*(r + x) - 4*(r^2 + r*x + x^2))
    end

    # Eq (41)
    function g_ν_e(x)
        # 將 (1 - x) 只乘在前兩項
        term12 = (1 - x) * (6 * (1 - x)^2 + r * (5 + 5*x - 4*x^2))
        term3 = 6 * r * log(x)
        return 2 / (3 * (1 - r)^2) * (term12 + term3)
    end

    # Eq (42)
    function h_ν_e_1(x)
        return 2 / (3 * (1 - r)^2) * ((1 - r) * (6 - 7*r + 11*r^2 - 4*r^3) + 6*r*log(r))
    end

    # Eq (43)
    function h_ν_e_2(x)
        return (2 * (r - x)) / (3 * r^2) * (7*r^2 - 4*r^3 + 7*x*r - 4*x*r^2 - 2*x^2 - 4*x^2*r)
    end

    function ν_spec(E_ν)
        ν_e_spec = quadgk((x) -> (h_ν_e_1(x) + h_ν_e_2(x)) * π_spec(E_ν / x) / x, 0, r)[1] +
            quadgk((x) -> g_ν_e(x) * π_spec(E_ν / x) / x, r, 1)[1]
        ν_μ_spec_1 = 1 / λ * quadgk((x) -> π_spec(E_ν / x) / x, 0, λ)[1]
        ν_μ_spec_2 = quadgk((x) -> (h_ν_μ_1(x) + h_ν_μ_2(x)) * π_spec(E_ν / x) / x, 0, r)[1] +
            quadgk((x) -> g_ν_μ(x) * π_spec(E_ν / x) / x, r, 1)[1]
        return ν_e_spec + ν_μ_spec_1 + ν_μ_spec_2
    end

    return ν_spec

end

function get_HNL_spec(π_spec, s; m_N=10)
    λ_v = λ_func(M_PI^2, m_N^2, M_MU^2)
    λ_min = (M_PI^2 - M_MU^2 + m_N^2 - sqrt(λ_v)) / (2 * M_PI^2)
    λ_max = (M_PI^2 - M_MU^2 + m_N^2 + sqrt(λ_v)) / (2 * M_PI^2)
    λ_diff = sqrt(λ_v) / M_PI^2

    Γ_p = ((M_MU^2 + m_N^2) * M_PI^2 - (M_MU^2 - m_N^2)^2)
    Γ_n = - (M_MU^2 - m_N^2) * sqrt(λ_v)
    ρ = Γ_p * sqrt(λ_v) / (M_MU^2 * (M_PI^2 - M_MU^2)^2)

    function HNL_spec(E_N)
        cosθ_r(x) = (x * (M_PI^2 - M_MU^2 + m_N^2) - 2 * m_N^2) / (x * sqrt((1 - (m_N/E_N)^2) * λ_v))
        r(x) = 1/2 + Γ_n / 2 / Γ_p * s * cosθ_r(x)
        return quadgk((x) -> π_spec(E_N / x) * ρ * r(x) * 1 / (λ_diff * x), λ_min, λ_max)[1]
    end

    return HNL_spec
end

function get_e_spec(HNL_spec, s)
    g(y) = 16 * y^3 - 27^y^2 + 11 - s * (32 * y^3 - 63 * y^2 + 36*y - 5)
    function e_spec(E)
        return 0.11 * quadgk((y) -> g(y) * HNL_spec(E / y) / (3 * y), 0, 1)[1]
    end
end

function get_μ_spec(π_spec)
    λ_n = M_MU^2 / M_PI^2

    function μ_spec(E)
        return quadgk((x) -> π_spec(E/x) / x, λ_n, 1)[1] / (1 - λ_n)
    end
end

function get_HNL_spec_μ(μ_spec, s; m_N=10)
    r_N = m_N / M_MU
    ρ = - r_N^8 + 8 * r_N^6 - 24 * r_N^4 * log(r_N) - 8 * r_N^2 + 1

    function HNL_spec(E_N)
        k(y) = -((M_MU^2*(y-1)*(M_MU^2*y-m_N^2)*(4*E_N*M_MU^4*s*y^4-4*sqrt(E_N-m_N)*sqrt(m_N+E_N)*M_MU^4*y^4-5*E_N*M_MU^4*s*y^3+E_N*m_N^2*M_MU^2*s*y^3+5*sqrt(E_N-m_N)*sqrt(m_N+E_N)*M_MU^4*y^3+5*sqrt(E_N-m_N)*m_N^2*sqrt(m_N+E_N)*M_MU^2*y^3-5*E_N*M_MU^4*s*y^2-8*E_N*m_N^2*M_MU^2*s*y^2+E_N*m_N^4*s*y^2+5*sqrt(E_N-m_N)*sqrt(m_N+E_N)*M_MU^4*y^2-22*sqrt(E_N-m_N)*m_N^2*sqrt(m_N+E_N)*M_MU^2*y^2+5*sqrt(E_N-m_N)*m_N^4*sqrt(m_N+E_N)*y^2+19*E_N*m_N^2*M_MU^2*s*y+E_N*m_N^4*s*y+5*sqrt(E_N-m_N)*m_N^2*sqrt(m_N+E_N)*M_MU^2*y+5*sqrt(E_N-m_N)*m_N^4*sqrt(m_N+E_N)*y-8*E_N*m_N^4*s-4*sqrt(E_N-m_N)*m_N^4*sqrt(m_N+E_N)))/(3*sqrt(E_N-m_N)*sqrt(m_N+E_N)*(24*m_N^4*M_MU^4*log(M_MU)+M_MU^8-8*m_N^2*M_MU^6-24*m_N^4*log(m_N)*M_MU^4+8*m_N^6*M_MU^2-m_N^8)*y^3))
        return quadgk((y) -> μ_spec(E_N / y) * k(y) * ρ / (2 * y), m_N^2 / M_MU^2, 1)[1]
    end

    return HNL_spec
end

# %% fit PI

ν_nums = ν_spec_func.(E_ν_array)
par = ()
min_delta = Inf
for A in 14.2:0.01:14.4
    for Γ in 0.10:0.01:0.13
        for E_cut in 0.7:0.01:0.9
            for β in 0.2:0.01:0.4
                π_spec_GeV = gen_π_spec_func(10^-A, Γ, 10^E_cut, β)
                ν_spec_pre = get_ν_spec((E) -> π_spec_GeV(E * 1e-3))

                ν_nums_pre = ν_spec_pre.(E_ν_array .* 1e3)

                delta = sum([(log(a/b))^2 for (a, b, E) in zip(ν_nums, ν_nums_pre, E_ν_array)])

                if delta < min_delta
                    par = (A=10^-A, Γ, E_cut=10^E_cut, β)
                    min_delta = delta
                end
            end
        end
    end
end

print(par)

open("./julia_results/pion_spec_par.txt", "w") do f
    write(f, string(par))
end

p_fit = plot(lit_ν_spec[:, 1], lit_ν_spec[:, 2] .* 1e-3, label="Neutrino (Theory)", color=:black, xscale=:log10, yscale=:log10, ylims=(1e-14, 1e-10), xlims=(1e-1, 1e5), linestyle=:dot, framestyle=:box, xlabel=L"E \; \mathrm{(GeV)}", ylabel=L"E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})", legend=:outerright)

# plot!(p_fit, E_ν_array, ((E) -> ν_spec_func(E) * E^2).(E_ν_array), label="int_ν")

π_spec_GeV = gen_π_spec_func(par.A, par.Γ, par.E_cut, par.β)

plot!(p_fit, (E) -> π_spec_GeV(E) * E^2, label="π_prox")

ν_spec_pre = get_ν_spec((E) -> π_spec_GeV(E * 1e-3))
plot!(p_fit, (E) -> ν_spec_pre(E * 1e3) * E^2, label="ν_prox")

# %% run HNL spec PI

π_spec_GeV = gen_π_spec_func(par.A, par.Γ, par.E_cut, par.β)

p_HNL_π_test = plot((E) -> π_spec_GeV(E) * E^2, label="π_prox", xscale=:log10,
    yscale=:log10, ylims=(1e-14, 1e-10), xlims=(1e-1, 1e5),
    linestyle=:dot, framestyle=:box, xlabel=L"E \; \mathrm{(GeV)}",
    ylabel=L"E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})",
    legend=:outerright)

ν_spec_pre = get_ν_spec((E) -> π_spec_GeV(E * 1e-3))
plot!(p_HNL_π_test, (E) -> ν_spec_pre(E * 1e3) * E^2, label="ν_prox")

# plot!(p_HNL_π_test ,lit_ν_spec[:,1], lit_ν_spec[:, 2] .* 1e-3, label="Neutrino (Theory)", color = :black)

HNL_spec_π_1 = get_HNL_spec((E) -> π_spec_GeV(E * 1e-3), 1; m_N=10)
# plot!(p_HNL_π_test, (E) -> HNL_spec_π_1(E * 1e3) * E^2, label = "HNL_1_prox")

HNL_spec_π_m1 = get_HNL_spec((E) -> π_spec_GeV(E * 1e-3), -1; m_N=10)
# plot!(p_HNL_π_test, (E) -> HNL_spec_π_m1(E * 1e3) * E^2, label = "HNL_-1_prox")

μ_spec = get_μ_spec((E) -> π_spec_GeV(E * 1e-3))
plot!(p_HNL_π_test, (E) -> μ_spec(E * 1e3) * E^2, label="μ spec")

HNL_spec_μ_1 = get_HNL_spec_μ(μ_spec, 1)
# plot!(p_HNL_π_test, (E) -> HNL_spec_μ_1(E * 1e3) * E^2, label = "HNL μ spec")

HNL_spec_μ_m1 = get_HNL_spec_μ(μ_spec, -1)
# plot!(p_HNL_π_test, (E) -> HNL_spec_μ_m1(E * 1e3) * E^2, label = "HNL μ spec")


# %% approx e by HNL

e_spec_1 = get_e_spec(HNL_spec_π_1, 1)

plot!(p_HNL_π_test, (E) -> e_spec_1(E * 1e3) * E^2, label="e (s=1)")

e_spec_m1 = get_e_spec(HNL_spec_π_m1, -1)

plot!(p_HNL_π_test, (E) -> e_spec_m1(E * 1e3) * E^2, label="e (s=-1)")
# %% approx γ

HNL_spec = (E) ->
    if E < E_interp_min || E > E_interp_max
        0.0
    else
        ν_spec_func(E) * 2 * 10^-2
    end

γ_spec = get_γ_spec_approx(HNL_spec; E_range=(E_interp_min, E_interp_max))

pgfplotsx()

p2 = plot(lit_ν_spec[:, 1], lit_ν_spec[:, 2] .* 1e-3, label="Neutrino (Theory)", color=:black, xscale=:log10, yscale=:log10, ylims=(1e-14, 1e-10), xlims=(1e-1, 1e5), linestyle=:dot, framestyle=:box, size=(400, 300), xlabel=L"E \; \mathrm{(GeV)}", ylabel=L"E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})", legend=:outerright)
plot!(p2, icecube_spec[:, 1], icecube_spec[:, 2], label="Neutrino (IceCube)", color=:blue)
plot!(p2, (E) -> γ_spec(E) * E^2, label="HNL γ", color=:red)
scatter!(p2, low_γ_obs[:, 1], low_γ_obs[:, 2],
    label="γ-ray 0.1 to 100 GeV", color=:green)
scatter!(p2, high_γ_obs[:, 1], high_γ_obs[:, 2],
    label="γ-ray > 200 GeV", color=:yellow)

savefig(p2, "julia_plots/u1.5_gamma.pdf")

# %% Get Plots
palette = ColorSchemes.okabe_ito

plt = @pgf LogLogAxis(
    {
    width = "8.5cm",
    height = "6.5cm",
    grid = "major",
    grid_style = {dashed, gray!20},
    axis_lines = "box",
    axis_line_style = {black, line_width = "1pt"},

    # 軸標籤與範圍
    xlabel = L"E \; \mathrm{(GeV)}",
    ylabel = L"$E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})$",
    xmin = 0.1,
    xmax = 1e6,
    ymin = 1e-14,
    ymax = 1e-10,

    # 刻度樣式 (出版級細節)
    tick_pos = "left",
    # xtick = "{0.0, 0.1, 0.2, 0.3, 0.4, 0.5}",
    minor_x_tick_num = 1,
    minor_y_tick_num = 1,
    axis_line_shift = "0pt",

    # 圖例樣式
    legend_style = {
        at = "(0.03, 0.97)",
        anchor = "north west",
        fill = "white",
        fill_opacity = 0.85,
        text_opacity = 1.0,
        draw = "gray!50",
        font = L"\small",
        nodes = "{scale=0.7, transform shape}"
    },
    legend_cell_align = "left"
}
)

E_array = 10 .^ (-1:0.05:6)
π_v = π_spec_GeV.(E_array) .* (E_array .^ 2)
plot_options = @pgf {
    color = palette[1],
    line_width = "1.0pt",
    no_marks
}
@pgf push!(plt, Plot(plot_options, Coordinates(E_array, π_v)))
@pgf push!(plt, LegendEntry(L"$\pi^\pm$ fit"))

ν_v = ν_spec_pre.(E_array .* 1e3) .* (E_array .^ 2)
plot_options = @pgf {
    color = palette[2],
    line_width = "1.0pt",
    no_marks
}
@pgf push!(plt, Plot(plot_options, Coordinates(E_array, ν_v)))
@pgf push!(plt, LegendEntry(L"$\nu$ fit"))

plot_options = @pgf {
    color = "black",
    line_width = "1.0pt",
    no_marks,
    densely_dashed
}
@pgf push!(plt, Plot(plot_options, Coordinates(lit_ν_spec[:, 1], lit_ν_spec[:, 2] .* 1e-3)))
@pgf push!(plt, LegendEntry(L"$\nu$ Murase (2020)"))

pgfsave("julia_plots/pion_spec_show.pdf", plt)

# %% Get Plot 2 Prepare

palette = ColorSchemes.okabe_ito

E_array = 10 .^ (-1:0.05:6)
π_v = π_spec_GeV.(E_array) .* (E_array .^ 2)

μ_v = μ_spec.(E_array .* 1e3) .* (E_array .^ 2)

N_π_m1_v = HNL_spec_π_m1.(E_array .* 1e3) .* (E_array .^ 2)

N_π_1_v = HNL_spec_π_1.(E_array .* 1e3) .* (E_array .^ 2)

N_μ_m1_v = HNL_spec_μ_m1.(E_array .* 1e3) .* (E_array .^ 2)

N_μ_1_v = HNL_spec_μ_1.(E_array .* 1e3) .* (E_array .^ 2)

df_s1 = DataFrame(E=E_array, Pion=π_v, Muon=μ_v)

CSV.write("julia_results/E_Pion_Muon.csv", df_s1)

# %% Get Plot 2

plt2 = @pgf LogLogAxis(
    {
    width = "8.5cm",
    height = "6.5cm",
    grid = "major",
    grid_style = {dashed, gray!20},
    axis_lines = "box",
    axis_line_style = {black, line_width = "1pt"},

    # 軸標籤與範圍
    xlabel = L"E \; \mathrm{(GeV)}",
    ylabel = L"$E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})$",
    xmin = 0.1,
    xmax = 1e6,
    ymin = 1e-14,
    ymax = 1e-10,

    # 刻度樣式 (出版級細節)
    tick_pos = "left",
    # xtick = "{0.0, 0.1, 0.2, 0.3, 0.4, 0.5}",
    minor_x_tick_num = 1,
    minor_y_tick_num = 1,
    axis_line_shift = "0pt",

    # 圖例樣式
    legend_style = {
        at = "(0.03, 0.97)",
        anchor = "north west",
        fill = "white",
        fill_opacity = 0.85,
        text_opacity = 1.0,
        draw = "gray!50",
        font = L"\small",
        nodes = "{scale=0.7, transform shape}"
    },
    legend_cell_align = "left"
}
)

plot_options = @pgf {
    color = "gray",
    line_width = "0.8pt"
}
@pgf push!(plt2, Plot(plot_options, Coordinates(E_array, π_v)))
@pgf push!(plt2, LegendEntry(L"$\pi^\pm$"))

plot_options = @pgf {
    color = "gray",
    line_width = "0.8pt",
    densely_dotted
}
@pgf push!(plt2, Plot(plot_options, Coordinates(E_array, μ_v)))
@pgf push!(plt2, LegendEntry(L"$\mu^\pm$"))

plot_options = @pgf {
    color = palette[1],
    line_width = "0.8pt"
}
@pgf push!(plt2, Plot(plot_options, Coordinates(E_array, N_π_m1_v)))
@pgf push!(plt2, LegendEntry(L"$N^{(\pi)},\,s=-1$"))

plot_options = @pgf {
    color = palette[1],
    line_width = "0.8pt",
    densely_dashed
}
@pgf push!(plt2, Plot(plot_options, Coordinates(E_array, N_π_1_v)))
@pgf push!(plt2, LegendEntry(L"$N^{(\pi)},\,s=+1$"))

plot_options = @pgf {
    color = palette[2],
    line_width = "0.8pt"
}
@pgf push!(plt2, Plot(plot_options, Coordinates(E_array, N_μ_m1_v)))
@pgf push!(plt2, LegendEntry(L"$N^{(\mu)},\,s=-1$"))

plot_options = @pgf {
    color = palette[2],
    line_width = "0.8pt",
    densely_dashed
}
@pgf push!(plt2, Plot(plot_options, Coordinates(E_array, N_μ_1_v)))
@pgf push!(plt2, LegendEntry(L"$N^{(\mu)},\,s=+1$"))

plot_options = @pgf {
    color = palette[3],
    line_width = "0.8pt"
}
e_v = CSV.read("julia_results/e_spec.csv", DataFrame)[:, 2]
@pgf push!(plt2, Plot(plot_options, Coordinates(E_array, e_v)))
@pgf push!(plt2, LegendEntry(L"overall $e^\pm$"))

pgfsave("julia_plots/HNL_spec_sep_show.pdf", plt2)

# %% Get Plot Electron Spectrum Showcase Prepare

N_m1_all_v = N_π_m1_v .+ N_μ_m1_v

N_1_all_v = N_π_1_v .+ N_μ_1_v

e_spec_π_1 = get_e_spec(HNL_spec_π_1, 1)

e_spec_π_m1 = get_e_spec(HNL_spec_π_m1, -1)

e_spec_μ_1 = get_e_spec(HNL_spec_μ_1, 1)

e_spec_μ_m1 = get_e_spec(HNL_spec_μ_m1, -1)

e_π_m1_v = e_spec_π_m1.(E_array .* 1e3) .* (E_array .^ 2)

e_π_1_v = e_spec_π_1.(E_array .* 1e3) .* (E_array .^ 2)

e_μ_m1_v = e_spec_μ_m1.(E_array .* 1e3) .* (E_array .^ 2)

e_μ_1_v = e_spec_μ_1.(E_array .* 1e3) .* (E_array .^ 2)

e_all_v = e_π_m1_v .+ e_μ_m1_v .+ e_π_1_v .+ e_μ_1_v

df = DataFrame(E=E_array, e_Phi=e_all_v)

CSV.write("julia_results/e_spec.csv", df)

# %% Get Plot Electron Spectrum Showcase

plt3 = @pgf LogLogAxis(
    {
    width = "8.5cm",
    height = "6.5cm",
    grid = "major",
    grid_style = {dashed, gray!20},
    axis_lines = "box",
    axis_line_style = {black, line_width = "1pt"},

    # 軸標籤與範圍
    xlabel = L"E \; \mathrm{(GeV)}",
    ylabel = L"$E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})$",
    xmin = 0.1,
    xmax = 1e6,
    ymin = 1e-14,
    ymax = 1e-10,

    # 刻度樣式 (出版級細節)
    tick_pos = "left",
    # xtick = "{0.0, 0.1, 0.2, 0.3, 0.4, 0.5}",
    minor_x_tick_num = 1,
    minor_y_tick_num = 1,
    axis_line_shift = "0pt",

    # 圖例樣式
    legend_style = {
        at = "(0.03, 0.97)",
        anchor = "north west",
        fill = "white",
        fill_opacity = 0.85,
        text_opacity = 1.0,
        draw = "gray!50",
        font = L"\small",
        nodes = "{scale=0.7, transform shape}"
    },
    legend_cell_align = "left"
}
)

palette = ColorSchemes.okabe_ito

plot_options = @pgf {
    color = "gray",
    line_width = "0.8pt"
}
@pgf push!(plt3, Plot(plot_options, Coordinates(E_array, N_m1_all_v)))
@pgf push!(plt3, LegendEntry(L"$N,\,s=-1$"))

plot_options = @pgf {
    color = "gray",
    line_width = "0.8pt",
    densely_dashed
}
@pgf push!(plt3, Plot(plot_options, Coordinates(E_array, N_1_all_v)))
@pgf push!(plt3, LegendEntry(L"$N,\,s=+1$"))

plot_options = @pgf {
    color = palette[1],
    line_width = "0.8pt"
}
@pgf push!(plt3, Plot(plot_options, Coordinates(E_array, e_all_v)))
@pgf push!(plt3, LegendEntry(L"$e$"))

pgfsave("julia_plots/HNL_to_e_show.pdf", plt3)

# %% different m_N
E_array = 10 .^ (-1:0.05:6)

π_spec_norm = (E) -> π_spec_GeV(E * 1e-3)

μ_spec = get_μ_spec(π_spec_norm)

m_N_list = range(10.0, 30.0, 21)

@showprogress @threads for m_N in m_N_list

    l_HNL_spec_π_1 = get_HNL_spec(π_spec_norm, 1; m_N=m_N)
    l_HNL_spec_π_m1 = get_HNL_spec(π_spec_norm, -1; m_N=m_N)

    l_HNL_spec_μ_1 = get_HNL_spec_μ(μ_spec, 1; m_N=m_N)
    l_HNL_spec_μ_m1 = get_HNL_spec_μ(μ_spec, -1; m_N=m_N)

    l_HNL_all_1 = (E) -> l_HNL_spec_π_1(E) + l_HNL_spec_μ_1(E)
    l_HNL_all_m1 = (E) -> l_HNL_spec_π_m1(E) + l_HNL_spec_μ_m1(E)

    l_e_spec_1 = get_e_spec(l_HNL_all_1, 1)
    l_e_spec_m1 = get_e_spec(l_HNL_all_m1, -1)

    l_e_spec_all = (E) -> l_e_spec_1(E) + l_e_spec_m1(E)

    e_spec_v = l_e_spec_all.(E_array .* 1e3)

    df = DataFrame(E=E_array, e_spec=e_spec_v)

    CSV.write("julia_results/e_spec_mN_" * string(m_N) * ".csv", df)

end

# %% Get Plot E

df1 = CSV.read("julia_results/E_Pion_Muon.csv", DataFrame)
# df2 = CSV.read("py_code/gamma_spec_30.0.txt", DataFrame; header=false)

dfs = [CSV.read("py_code/gamma_spec_$i.txt", DataFrame; header=false) for i in 10.0:1.0:30.0]

E_γs = dfs[1][:, 1]
γ_ϕs_max = zeros(length(E_γs))
γ_ϕs_min = zeros(length(E_γs))

@showprogress @threads for i in range(1, length(E_γs))
    tmp_vec = [df[i, 2] for df in dfs]
    max_phi = maximum(tmp_vec)
    min_phi = minimum(tmp_vec)
    γ_ϕs_max[i] = max_phi
    γ_ϕs_min[i] = min_phi
end

pe = plot(df1[:, 1], df1[:, 2], label=L"$e^\pm$", color=:black, xscale=:log10, yscale=:log10, ylims=(1e-19, 1e-10), xlims=(1e-6, 1e5), linestyle=:dot, framestyle=:box, size=(400, 300), xlabel=L"E \; \mathrm{(GeV)}", ylabel=L"E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})", legend=:outerright)

plot!(pe, E_γs, γ_ϕs_max .* (E_γs .^ 2), label=L"$\gamma$")
plot!(pe, E_γs, γ_ϕs_min .* (E_γs .^ 2), label=L"$\gamma$")

# %% Plot gamma show case

push!(PGFPlotsX.CUSTOM_PREAMBLE, raw"\usepgfplotslibrary{fillbetween}")
plt4 = @pgf LogLogAxis(
    {
    width = "8.5cm",
    height = "6.5cm",
    grid = "major",
    grid_style = {dashed, gray!20},
    axis_lines = "box",
    axis_line_style = {black, line_width = "1pt"},

    # 軸標籤與範圍
    xlabel = L"E \; \mathrm{(GeV)}",
    ylabel = L"$E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})$",
    xmin = 0.001,
    xmax = 1e6,
    ymin = 1e-19,
    ymax = 1e-10,

    # 刻度樣式 (出版級細節)
    tick_pos = "left",
    # xtick = "{0.0, 0.1, 0.2, 0.3, 0.4, 0.5}",
    minor_x_tick_num = 1,
    minor_y_tick_num = 1,
    axis_line_shift = "0pt",

    # 圖例樣式
    legend_style = {
        at = "(1.03, 0.97)",
        anchor = "north west",
        fill = "white",
        fill_opacity = 0.85,
        text_opacity = 1.0,
        draw = "gray!50",
        font = L"\small",
        nodes = "{scale=0.7, transform shape}"
    },
    legend_cell_align = "left"
}
)
palette = ColorSchemes.okabe_ito

plot_options = @pgf {
    "name path=f",
    color = "red",
    line_width = "0.5pt"
}
@pgf push!(plt4, Plot(plot_options, Coordinates(E_γs, γ_ϕs_max .* (E_γs .^ 2))))
@pgf push!(plt4, LegendEntry(L"$\gamma ,\times |U_\mu|^{-2}$  HNL"))

plot_options = @pgf {
    "forget plot",
    "name path=g",
    color = "red",
    line_width = "0.5pt"
}
@pgf push!(plt4, Plot(plot_options, Coordinates(E_γs, γ_ϕs_min .* (E_γs .^ 2))))

plot_options = @pgf {
    "forget plot",
    fill = "red",
    draw = "none",
    opacity = 0.4
}
@pgf push!(plt4, Plot(plot_options, raw"fill between [of=f and g]"))

plot_options = @pgf {
    color = "gray",
    line_width = "0.8pt"
}
@pgf push!(plt4, Plot(plot_options, Coordinates(lit_ν_spec[:, 1], lit_ν_spec[:, 2] .* 1e-3)))
@pgf push!(plt4, LegendEntry(L"$\nu_\mu$ Murase (2020)"))

plot_options = @pgf {
    only_marks,
    mark = "*",
    mark_size = "0.8pt",
    color = palette[3],
    "error bars/x dir" = "both",
    "error bars/x explicit" = true,
    "error bars/y dir" = "both",
    "error bars/y explicit" = true,
    "error bars/error bar style" = {line_width = "0.8pt"}
}
@pgf push!(plt4, Plot(plot_options, Coordinates(
    low_γ_obs[1:6, 1], low_γ_obs[1:6, 2];
    xerrorminus=low_γ_obs[1:6, 3],
    xerrorplus=low_γ_obs[1:6, 4],
    yerrorminus=low_γ_obs[1:6, 5],
    yerrorplus=low_γ_obs[1:6, 6]
)))
@pgf push!(plt4, LegendEntry(L"$\gamma$ 4FGL-DR2"))
plot_options = @pgf {
    "forget plot",
    only_marks,
    mark_size = "0.8pt",
    color = palette[3],
    "error bars/x dir" = "both",
    "error bars/x explicit" = true,
    # "error bars/y dir" = "both",
    # "error bars/y explicit" = true,
    "error bars/error bar style" = {line_width = "0.8pt"},
    # mark="triangle*",
    # "every mark/.append style"="{rotate=180}"
}
@pgf push!(plt4, Plot(plot_options, Coordinates(
    low_γ_obs[7:7, 1], low_γ_obs[7:7, 2];
    xerrorminus=low_γ_obs[7:7, 3],
    xerrorplus=low_γ_obs[7:7, 4],
    yerrorminus=low_γ_obs[7:7, 5],
    yerrorplus=low_γ_obs[7:7, 6]
)))

plot_options = @pgf {
    only_marks,
    mark_size = "0.8pt",
    color = palette[1],
    "error bars/x dir" = "both",
    "error bars/x explicit" = true,
    # "error bars/y dir" = "both",
    # "error bars/y explicit" = true,
    "error bars/error bar style" = {line_width = "0.8pt"},
    # mark="triangle*",
    # "every mark/.append style"="{rotate=180}"
}
@pgf push!(plt4, Plot(plot_options, Coordinates(
    high_γ_obs[:, 1], high_γ_obs[:, 2];
    xerrorminus=high_γ_obs[:, 3],
    xerrorplus=high_γ_obs[:, 4]
)))
@pgf push!(plt4, LegendEntry(L"$\gamma$ MAGIC"))

pgfsave("julia_plots/gamma_show_compare.pdf", plt4)

# %% constrain U

obs_γ_UL = [low_γ_obs[:, 1:2]; high_γ_obs[:, 1:2]]

function get_U_UL(m_N)
    ρ_π = (M_PI^2 * (M_MU^2 + m_N^2) - (M_MU^2 - m_N^2)^2) * sqrt(λ_func(M_PI^2, M_MU^2, m_N^2)) / (M_MU^2 * (M_PI^2 - M_MU^2)^2)
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

    predic_Φ = [γ_spec_approx(E) * E^2 for E in obs_γ_UL[:, 1]]
    U_UL = minimum([obs_γ_UL[i, 2] / predic_Φ[i] for i in eachindex(predic_Φ) if predic_Φ[i] > 0])

    return U_UL

end

m_N_range = range(1e-3, 3e-2; length=100)
U_UL_list = zeros(length(m_N_range))
@showprogress @threads for i in eachindex(U_UL_list)
    U_UL_list[i] = get_U_UL(m_N_range[i])
end
p3 = plot(m_N_range .* 1e3, U_UL_list, color=:black, framestyle=:box, size=(400, 300), label=nothing, xlabel=L"m_N \; \mathrm{(MeV)}", ylabel=L"|U_\mu|^2")
savefig(p3, "julia_plots/U_UL.pdf")


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

p4 = plot(lit_ν_spec[:, 1], lit_ν_spec[:, 2] .* 1e-3, label="Neutrino (Theory)", color=:black, xscale=:log10, yscale=:log10, ylims=(1e-16, 1e-10), xlims=(1e-1, 1e5), linestyle=:dot, framestyle=:box, size=(400, 300), xlabel=L"E \; \mathrm{(GeV)}", ylabel=L"E^2 \Phi \; (\mathrm{TeV} \;  \mathrm{cm}^{-2} \;  \mathrm{s}^{-1})", legend=:outerright)
plot!(p4, icecube_spec[:, 1], icecube_spec[:, 2], label="Neutrino (IceCube)", color=:blue)
plot!(p4, (E) -> e_spec(E) * E^2, label="HNL e", color=:green)
plot!(p4, (E) -> γ_N_spec(E) * E^2, label="HNL γ", color=:red)
scatter!(p4, low_γ_obs[:, 1], low_γ_obs[:, 2],
    label="γ-ray 0.1 to 100 GeV", color=:green)
scatter!(p4, high_γ_obs[:, 1], high_γ_obs[:, 2],
    label="γ-ray > 200 GeV", color=:yellow)

savefig(p4, "julia_plots/u_e_gamma.pdf")