# %% Prepare
using PGFPlotsX
using LaTeXStrings
using ColorSchemes
using Plots

push!(PGFPlotsX.CUSTOM_PREAMBLE, "\\usepackage{amsmath, amsfonts}")

# %% Formula

M_MU=105.6583755
M_PI=139.57039
STEP=0.001

function λ(x, y, z)
    x^2 + y^2 + z^2 - 2 * x * y - 2 * y * z - 2 * z * x
end

function gen_cosθ(m_N; m_π=M_PI, m_μ=M_MU)
    x_max = (m_π^2 - m_μ^2 + m_N^2 + sqrt(λ(m_π^2, m_μ^2, m_N^2))) / (2 * m_π^2)
    x_min = (m_π^2 - m_μ^2 + m_N^2 - sqrt(λ(m_π^2, m_μ^2, m_N^2))) / (2 * m_π^2)
        
    function cosθ_r(x)
        if x > x_max || x < x_min
            return 0
        end
        return (x * (m_π^2 - m_μ^2 + m_N^2) - 2 * m_N^2) / x / sqrt(λ(m_π^2, m_μ^2, m_N^2))
    end

    return (x_min, x_max), cosθ_r
end

# %% Plot

rgys = []

for i in [10, 20, 30]
    (xmin, xmax), f = gen_cosθ(i)
    rg=xmin:STEP:xmax
    push!(rgys, (xmin:STEP:xmax, f.(rg)))
end

p1 = plot(rgys[1][1], rgys[1][2])
plot!(p1, rgys[2][1], rgys[2][2])
plot!(p1, rgys[3][1], rgys[3][2])


# %% PGF

palette = ColorSchemes.okabe_ito
curves = [(color=palette[1], m=10.0, label=L"$m_N = 10 \text{MeV}$"),
          (color=palette[2], m=20.0, label=L"$m_N = 20 \text{MeV}$"),
          (color=palette[3], m=30.0, label=L"$m_N = 30 \text{MeV}$")]

plt = @pgf Axis(
    {
        width = "8.5cm",
        height = "6.5cm",
        grid = "major",
        grid_style = {dashed, gray!20},

        axis_lines = "box",
        axis_line_style = {black, line_width = "1pt"},
        
        # 軸標籤與範圍
        xlabel = L"$E_N / E_\pi$",
        ylabel = L"$\cos\theta_r$",
        xmin = 0.0,
        xmax = 0.45,
        ymin = -1.05,
        ymax = 1.05,
        
        # 刻度樣式 (出版級細節)
        tick_pos = "left",
        xtick = "{0.0, 0.1, 0.2, 0.3, 0.4, 0.5}",
        minor_x_tick_num = 1,
        minor_y_tick_num = 1,
        axis_line_shift = "0pt",
        
        # 圖例樣式
        legend_style = {
            at = "(0.97, 0.03)",
            anchor = "south east",
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

# 添加曲線
for c in curves
    ((xmin, xmax), f) = gen_cosθ(c.m)
    rg = xmin:STEP:xmax
    ys = f.(rg)

    plot_options = @pgf {
        color = c.color,
        line_width = "1.2pt",
        no_marks
    }

    dash_line_1 = (x = [xmin, xmin], y = [ys[1], 0])
    dash_line_2 = (x = [xmax, xmax], y = [ys[end], 0])

    @pgf push!(plt, Plot(plot_options, Coordinates(rg, ys)))
    @pgf push!(plt, LegendEntry(c.label))
end

for c in curves
    ((xmin, xmax), f) = gen_cosθ(c.m)
    rg = xmin:STEP:xmax
    ys = f.(rg)

    plot_options = @pgf {
        color = "gray",
        line_width = "0.5pt",
        no_marks,
        densely_dashed
    }

    dash_line_1 = (x = [xmin, xmin], y = [ys[1], 0])
    dash_line_2 = (x = [xmax, xmax], y = [ys[end], 0])

    @pgf push!(plt, Plot(plot_options, Coordinates(dash_line_1.x, dash_line_1.y)))
    @pgf push!(plt, Plot(plot_options, Coordinates(dash_line_2.x, dash_line_2.y)))
end

plot_options = @pgf {
    color = "black",
    line_width = "0.5pt",
    no_marks
}
@pgf push!(plt, Plot(plot_options, Coordinates([0, 0.45], [0, 0])))

pgfsave("julia_plots/pion_decay_lr.pdf", plt)
pgfsave("julia_plots/pion_decay_lr.tex", plt)

# %% PGF_P_H

function gen_Γs(s; m_μ=M_MU, m_π=M_PI)
    function Γs(m_N)
        NOM = (m_μ^2 + m_N^2) * m_π^2 - (m_μ^2 - m_N^2)^2
        return (NOM - s * (m_μ^2 -m_N^2) * sqrt(λ(m_π^2, m_μ^2, m_N^2))) / (2 * NOM)
    end

    return Γs
end

rg = 0.0:0.01:(M_PI - M_MU)

Γn = gen_Γs(-1)
Γp = gen_Γs(+1)

ysn = Γn.(rg)
ysp = Γp.(rg)

# plot(rg, ysn)

plt2 = @pgf Axis(
    {
        width = "8.5cm",
        height = "6.5cm",
        grid = "major",
        grid_style = {dashed, gray!20},

        axis_lines = "box",
        axis_line_style = {black, line_width = "1pt"},
        
        # 軸標籤與範圍
        xlabel = L"$m_N \quad \text{[MeV]}$",
        ylabel = L"$BR^{(s)}$",
        xmin = 0.0,
        xmax = rg[end],
        ymin = 0.0,
        ymax = 1.0,
        
        # 刻度樣式 (出版級細節)
        tick_pos = "left",
        xtick = "{0, 5, 10, 15, 20, 25, 30}",
        minor_x_tick_num = 1,
        minor_y_tick_num = 1,
        axis_line_shift = "0pt",
        
        # 圖例樣式
        legend_style = {
            at = "(0.97, 0.03)",
            anchor = "south east",
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

plot_options_n = @pgf {
    color = palette[1],
    line_width = "1.2pt",
    no_marks
}

@pgf push!(plt2, Plot(plot_options_n, Coordinates(rg, ysn)))
@pgf push!(plt2, LegendEntry(L"$s = - 1$"))

plot_options_p = @pgf {
    color = palette[2],
    line_width = "1.2pt",
    no_marks
}
@pgf push!(plt2, Plot(plot_options_p, Coordinates(rg, ysp)))
@pgf push!(plt2, LegendEntry(L"$s = + 1$"))

pgfsave("julia_plots/pion_decay_hel_br.pdf", plt2)