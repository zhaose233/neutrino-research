# %% Prepare
using PGFPlotsX
using LaTeXStrings
using ColorSchemes
using Plots

push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{amsmath, amsfonts}")

# %% Formula Muon Decay Spectrum

M_MU=105.6583755

function gen_dΓdE(m_N, s; m_μ=M_MU)
    function dΓdE(x)
        E = x * m_μ
        if E < m_N || E > (m_μ^2 + m_N^2) / (2 * m_μ)
            return 0
        else
            p1 = sqrt(E^2 - m_N^2)
            return (8*(3*m_N^2*m_μ^2*s-3*E^2*m_μ^2*s-4*E*m_N^2*m_μ*s+4*E^3*m_μ*s+m_N^4*s-E^2*m_N^2*s-4*m_μ*p1^3+3*E*m_μ^2*p1-6*m_N^2*m_μ*p1+3*E*m_N^2*p1))/m_μ^4
        end
    end

    return ((min=m_N/m_μ, max=(m_μ^2 + m_N^2) / (2 * m_μ^2)), dΓdE)
end

# %% Plot
f01 = gen_dΓdE(0, -1)

f11 = gen_dΓdE(15, 1)
f12 = gen_dΓdE(15, -1)

f31 = gen_dΓdE(30, 1)
f32 = gen_dΓdE(30, -1)


rg = 0.0:0.0001:0.55
p1=plot(rg, f01.(rg))
plot!(p1,rg, f11.(rg))
plot!(p1,rg, f12.(rg))

plot!(p1,rg, f31.(rg))
plot!(p1,rg, f32.(rg))
savefig(p1, "julia_plots/mu_HNL_spec.png")

# %% PGF

palette = ColorSchemes.okabe_ito

plt = @pgf Axis(
    {
        width = "8.5cm",
        height = "6.5cm",
        grid = "major",
        grid_style = {dashed, gray!20},

        axis_lines = "box",
        axis_line_style = {black, line_width = "1pt"},
        
        # 軸標籤與範圍
        xlabel = L"$x = E / m_\mu$",
        ylabel = L"$\mathrm{d}\Gamma / \mathrm{d}x \quad [\text{a.u.}]$",
        xmin = 0.0,
        xmax = 0.56,
        ymin = 0.0,
        ymax = 4.2,
        
        # 刻度樣式 (出版級細節)
        tick_pos = "left",
        xtick = "{0.0, 0.1, 0.2, 0.3, 0.4, 0.5}",
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

curves = [
    (m_N=0,  s=-1, label=L"$m_N = 0\text{ MeV},\, s = -1$", color=palette[1]),
    (m_N=15, s=-1, label=L"$m_N = 15\text{ MeV},\, s = -1$", color=palette[2]),
    (m_N=15, s=1,  label=L"$m_N = 15\text{ MeV},\, s = +1$", color=palette[2]),
    (m_N=30, s=-1, label=L"$m_N = 30\text{ MeV},\, s = -1$", color=palette[3]),
    (m_N=30, s=1,  label=L"$m_N = 30\text{ MeV},\, s = +1$", color=palette[3])
]

# 添加曲線
for c in curves
    x_tuple, f = gen_dΓdE(c.m_N, c.s)
    rg = vcat(x_tuple.min:0.001:x_tuple.max, x_tuple.max + 1e-6)
    ys = f.(rg)

    plot_options = @pgf {
        color = c.color,
        line_width = "1.2pt",
        no_marks
    }
    if c.s == 1
        plot_options["densely dashed"] = true
    end

    @pgf push!(plt, Plot(plot_options, Coordinates(rg, ys)))
    @pgf push!(plt, LegendEntry(c.label))
end

pgfsave("julia_plots/muon_decay_spectrum.pdf", plt)
pgfsave("julia_plots/muon_decay_spectrum.tex", plt)
