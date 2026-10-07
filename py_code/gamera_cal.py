import sys

sys.path.append('../GAMERA/lib')
import gappa as gp
import numpy as np

for i in range(10, 31):

    E_raw, Phi_raw = np.loadtxt('../julia_results/e_spec_mN_' + str(i * 1.0) + '.csv', skiprows=1, delimiter=',', unpack=True)

    d_L = 14400 * gp.kpc_to_cm

    e = E_raw * gp.GeV_to_erg
    phi = Phi_raw / gp.GeV_to_erg * (d_L**2 * 4 * np.pi) * 1000

    fr = gp.Radiation()

    # parameter -----------------

    c_speed = 2.9979e10

    b_field = 0.139 # Gauss

    n_density = 2.3e5

    grain_temp = 1500
    L_bol = 9.55e44
    r_sub = 0.4 * gp.pc_to_cm
    U_IR = L_bol / (4 * np.pi * r_sub**2 * c_speed)
    T_dust = 1500

    print(U_IR)

    # --------------------------

    fr.SetBField(b_field)
    fr.AddThermalTargetPhotons(T_dust, U_IR)
    fr.SetAmbientDensity(n_density)
    fr.SetDistance(14.4e6)

    p = list(zip(e, phi))
    fr.SetElectrons(p)

    e_out = np.logspace(-6,6,300) * gp.GeV_to_erg

    fr.CalculateDifferentialPhotonSpectrum(e_out)

    gamma_spec = np.array(fr.GetTotalSpectrum())
    gamma_spec[:, 1] = gamma_spec[:, 1] / 1000.0 * gp.GeV_to_erg
    gamma_spec[:, 0] = gamma_spec[:, 0] / gp.GeV_to_erg

    np.savetxt('gamma_spec_' + str(i * 1.0) + '.txt', gamma_spec, delimiter=",")