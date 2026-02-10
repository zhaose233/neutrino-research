using Plots
using Printf

# ==========================================
# 1. 物理常数与基本单位
# ==========================================
const c_cgs = 2.99792458e10
const G_cgs = 6.67430e-8
const k_B   = 1.380649e-16
const m_p   = 1.6726e-24
const M_sun = 1.989e33
const Mpc_cgs = 3.086e24
const eV_to_erg = 1.60218e-12
const sigma_T = 6.65e-25

# ==========================================
# 2. NGC 1068 环境参数计算器
# ==========================================
function get_ngc1068_params()
    # 观测数据
    L_X = 10.0^43       # erg/s
    d_Mpc = 12.7        # Mpc
    
    # 标度关系 (来自论文 Table S1 & S2)
    # M_BH ~ L_X^0.746
    M_BH = 2.0e7 * M_sun * (L_X / 1.16e43)^0.746
    R_s = 2 * G_cgs * M_BH / c_cgs^2
    
    # 模型设定 R = 30 Rs
    R_norm = 30.0
    R_cor = R_norm * R_s
    
    # 估算 Coronal Density (n_p) 和 Magnetic Field (B)
    # 简化：使用论文中 implicit 的数值
    # 对于 Lx=10^43, B ~ 1000-2000 G (Table S2 says B~1.3 kG for beta=1)
    # n_p ~ 10^9 - 10^10 cm^-3
    # 这里我们反推 n_p 以保持光学深度一致
    log_L_X = log10(L_X)
    tau_T = 0.46 - 0.05 * (log_L_X - 44.0) # 约 0.5
    H = R_cor / sqrt(3)
    n_p = tau_T / (sigma_T * H)
    
    # 质子热能 (Virial Temperature)
    T_p = (m_p * c_cgs^2) / (6 * k_B * R_norm) # ~ 100 MeV
    
    # Infall velocity (Escape velocity proxy)
    V_K = c_cgs / sqrt(2 * R_norm)
    alpha = 0.1
    V_fall = alpha * V_K
    
    return L_X, d_Mpc, R_cor, n_p, T_p, V_fall, H
end

# ==========================================
# 3. 核心物理：计算归一化和截止能量
# ==========================================
function calc_spectrum_thin_line(E_nu_GeV_array)
    # --- 1. 获取环境参数 ---
    L_X, d_Mpc, R_cor, n_p, T_p, V_fall, H = get_ngc1068_params()
    
    # --- 2. 设定特定曲线的参数 (Thin Line in Fig 4) ---
    eta = 70.0              # 湍流参数 (影响截止能量)
    P_CR_ratio = 0.30       # 宇宙线压强比 (影响归一化, 30%)
    
    # --- 3. 计算截止能量 E_cut (依赖于 eta) ---
    # 物理逻辑：加速时间 t_acc = 冷却时间 t_cool
    # t_acc ∝ eta. t_cool (Bethe-Heitler) 是常数.
    # 因此 E_max ∝ 1/eta.
    # 参考点：当 eta=10 时，E_p,cl ≈ 220 TeV (Table S2 for Lx=10^43)
    # 所以当 eta=70 时：
    E_p_cl_ref = 220.0 # TeV (at eta=10)
    E_p_cut_TeV = E_p_cl_ref * (10.0 / eta) # 反比缩放
    
    # 中微子截止能量 ≈ 质子能量 / 20
    E_nu_cut_GeV = (E_p_cut_TeV * 1000.0) / 20.0
    
    # --- 4. 计算总光度归一化 (Luminosity) ---
    # 热压力 P_th
    P_th = n_p * k_B * T_p
    # 宇宙线压力 P_CR
    P_CR = P_CR_ratio * P_th
    
    # 宇宙线总注入光度 L_CR
    # L_CR ≈ Energy / t_esc ≈ (3 * P_CR * Vol) / (R/V_fall)
    Vol = 2 * pi * R_cor^2 * H
    t_esc = R_cor / V_fall
    L_CR = (3 * P_CR * Vol) / t_esc
    
    # --- 5. 效率因子 (f_pp vs f_BH) ---
    # 需要计算 f_pp 和 f_BH
    # 使用简化近似，基于 Lx=10^43
    # f_pp ~ 0.5 (from previous calc)
    # f_BH ~ 30.0 (from previous calc, dominated by MeV gamma rays)
    f_pp = 0.6 # 近似值，基于 Table S1
    f_BH = 30.0 
    
    # 中微子产生分支比
    # Efficiency = f_pp / (1 + f_BH + f_pp)
    # 还要乘以 1/8 (p -> pi -> nu 能量分配)
    f_eff = f_pp / (1.0 + f_BH + f_pp)
    L_nu_total = (1.0/8.0) * f_eff * L_CR
    
    # --- 6. 构建谱形 (Energy Dependence) ---
    # Fig 4 的 Thin Line 特征：
    # 低能区：像 E^-2 上升
    # 高能区：在 1-3 TeV 处指数截断
    # 公式：Flux(E) = A * E^(-Gamma) * exp(-E/E_cut)
    # 我们这里算的是 E^2 * Flux，所以应该是 E^(2-Gamma) * exp(...)
    
    Gamma = 2.0 # 标准费米加速谱指数
    
    fluxes = zeros(length(E_nu_GeV_array))
    d_cm = d_Mpc * Mpc_cgs
    
    # 归一化常数 A 的计算
    # L_nu_total = Integral[ Flux(E) * 4pi d^2 * dE ]
    # 对于 E^-2 谱，E^2 dN/dE 是常数 C0
    # Integral (C0/E^2 * E) dE = Integral (C0/E) dE = C0 * ln(Emax/Emin)
    # 所以 C0 ≈ L_nu / (4pi d^2 * ln(Emax/Emin))
    # Bolometric correction factor ln(Emax/Emin) ~ ln(1e6/1e0) ~ 14
    bolo_factor = 14.0 
    
    peak_flux_approx = L_nu_total / (4 * pi * d_cm^2 * bolo_factor)
    
    for (i, E_val) in enumerate(E_nu_GeV_array)
        # 谱形因子
        shape = (E_val / E_nu_cut_GeV)^(2 - Gamma) * exp(-E_val / E_nu_cut_GeV)
        
        # 低能修正 (Optional): 
        # 考虑到 pion production threshold，在 < 100 GeV 处急剧下降
        # 论文 Fig 4 在 100 GeV 以下断崖式下跌
        threshold_suppression = 1.0
        if E_val < 200.0
             # 简单的经验函数模拟阈值效应
             threshold_suppression = (E_val / 200.0)^2 
        end
        
        fluxes[i] = peak_flux_approx * shape * threshold_suppression
    end
    
    return fluxes, E_nu_cut_GeV, L_nu_total
end

# ==========================================
# 4. 绘图与执行
# ==========================================

E_grid = 10.0 .^ range(0, 6, length=200) # 1 GeV - 1 PeV
flux_thin, E_cut_val, L_tot = calc_spectrum_thin_line(E_grid)

# 转换单位: L_nu 是 erg/s -> Flux 是 erg cm^-2 s^-1
# 图中 Y 轴是 GeV cm^-2 s^-1
# 1 erg = 624.15 GeV
flux_thin_GeV = flux_thin .* 624.15e9

plot(E_grid, flux_thin_GeV, 
    xscale=:log10, yscale=:log10,
    lw=2, color=:black, linestyle=:solid, # Thin line in paper is thin, usually solid
    label="Model: \$\\eta=70, P_{CR}/P_{th}=0.3\$",
    xlabel="Neutrino Energy [GeV]",
    ylabel="E F_E [GeV cm^{-2} s^{-1}]",
    title="Reproduction of Fig 4 (Thin Black Line)",
    xlim=(1e0, 1e6), ylim=(1e-11, 1e-6),
    grid=true, legend=:topright)

# 标记特征点
vline!([E_cut_val], ls=:dash, color=:blue, label="Cutoff ~ $(round(E_cut_val/1000, digits=1)) TeV")
# 标记 IceCube 数据范围 (近似)
plot!([1e3, 1e5], [2e-8, 2e-8], color=:gray, alpha=0.5, label="IceCube Approx Level")

println("=== Calculation Results ===")
@printf("Neutrino Cutoff Energy: %.2f TeV\n", E_cut_val/1000)
@printf("Total Neutrino Luminosity: %.2e erg/s\n", L_tot)