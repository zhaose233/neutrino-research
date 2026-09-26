import sys

sys.path.append('../GAMERA/lib')
import gappa as gp
import numpy as np

E_raw, Phi_raw = np.loadtxt('../julia_results/e_spec.csv', skiprows=1, delimiter=',', unpack=True)

d_L = 14400 * gp.kpc_to_cm

e = E_raw * gp.GeV_to_erg
phi = Phi_raw / gp.GeV_to_erg * (d_L**2 * 4 * np.pi) * 1000

fp = gp.Particles()

fp.SetCustomInjectionSpectrum(list(zip(e, phi)))

fp.SetBField(10e-3)
fp.AddThermalTargetPhotons(500, 1e4 * gp.eV_to_erg)

fp.SetAge(1e4)

fp.CalculateElectronSpectrum()

sed_s = np.array(fp.GetParticleSED())

np.savetxt('out.csv', sed_s)