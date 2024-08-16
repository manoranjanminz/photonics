import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use(['science','no-latex', 'grid'])


# Optical Bistability in MRR
c = 299792458
n_0 = 2.03  # Refractive index
n_g = 2.03  # Group index
tau_linear = (1 / 1e9)  # Optical linear absorption rate in the cavity
n_2 = 2.4 * 1e-15 * 1e-7   # Kerr coefficient
V_kerr = 1.35 * 1e-16  # Cavity volume
dn_by_dT = 2.6 * 1e-5  # Thermo-optic coefficient
R = 120  # Thermal resistance
Q_in = 31400
Q_v = 21300
Q_l = 12700
lambda_0 = 1550e-9

# Initial values
delta_lambda_therm = 0
delta_lambda_kerr = 0
omega_0 = 2 * np.pi * c / (lambda_0 + delta_lambda_therm + delta_lambda_kerr)
P_in = np.array([0.001, 0.01, 0.1, 1, 5]) # in mW
tau_in = Q_in / omega_0
tau_v = Q_v / omega_0

omega = np.linspace(2 * np.pi * 193.62e12, 2 * np.pi * 193.20e12, 2000)
# omega_rev = np.linspace(2 * np.pi * 193.20e12, 2 * np.pi * 193.62e12, 2000)
photon_energy = np.zeros((len(P_in), len(omega)))

fixed_omega = 2 * np.pi * c / (lambda_0 + 0.1e-9)
delta_lambda_therm_fo = 0
delta_lambda_kerr_fo = 0
omega_0_fo = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_fo + delta_lambda_kerr_fo)
tau_in_fo = tau_in
tau_v_fo = tau_v
P_in_fo = np.linspace(1,15,200)
photon_energy_fixed_omega = np.zeros((len(P_in_fo), 1))

delta_lambda_therm_fo_rev = 0
delta_lambda_kerr_fo_rev = 0
omega_0_fo_rev = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_fo + delta_lambda_kerr_fo)
tau_in_fo_rev = tau_in
tau_v_fo_rev = tau_v
P_in_fo_rev = np.linspace(15,1,200)
photon_energy_fixed_omega_rev = np.zeros((len(P_in_fo_rev), 1))

n_itr = 20 # No. of iterations for a fixed omega

photon_energy_evo = np.zeros((len(P_in), len(omega), n_itr))
x_n_itr = np.linspace(1,n_itr,n_itr)
# Nonlinear CMT model
for i in range(len(P_in)):
    for j in range(len(omega)):
        for z in range(n_itr):
            photon_energy[i, j] = ((1 / tau_in) * (P_in[i] / 2)) / ((omega[j] - omega_0) ** 2 + (((1 / tau_in) + (1 / tau_v) + (1 / tau_linear)) ** 2) / 4)
            photon_energy_evo[i,j,z] = photon_energy[i,j]

            omega_0 = 2 * np.pi * c / (lambda_0 + delta_lambda_therm + delta_lambda_kerr)
            tau_in = Q_in / omega_0
            tau_v = Q_v / omega_0
            delta_lambda_therm = (lambda_0 / n_0) * dn_by_dT * photon_energy[i, j] * (1 / tau_linear) * R
            delta_lambda_kerr = (n_2 * c * photon_energy[i, j] * lambda_0) / (n_0 * n_g * V_kerr)

# Cavity Energy vs Input Power for a fixed omega
for i in range(len(P_in_fo)):
    for z in range(n_itr):
        photon_energy_fixed_omega[i] = ((1 / tau_in_fo) * (P_in_fo[i] / 2)) / ((fixed_omega - omega_0_fo) ** 2 + (((1 / tau_in_fo) + (1 / tau_v_fo) + (1 / tau_linear)) ** 2) / 4)
        
        omega_0_fo = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_fo + delta_lambda_kerr_fo)
        tau_in_fo = Q_in / omega_0_fo
        tau_v_fo = Q_v / omega_0_fo
        delta_lambda_therm_fo = (lambda_0 / n_0) * dn_by_dT * photon_energy_fixed_omega[i] * (1 / tau_linear) * R
        delta_lambda_kerr_fo = (n_2 * c * photon_energy_fixed_omega[i] * lambda_0) / (n_0 * n_g * V_kerr)


for i in range(len(P_in_fo_rev)):
    for z in range(n_itr):
        photon_energy_fixed_omega_rev[i] = ((1 / tau_in_fo_rev) * (P_in_fo_rev[i] / 2)) / ((fixed_omega - omega_0_fo_rev) ** 2 + (((1 / tau_in_fo_rev) + (1 / tau_v_fo_rev) + (1 / tau_linear)) ** 2) / 4)
        
        omega_0_fo_rev = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_fo_rev + delta_lambda_kerr_fo_rev)
        tau_in_fo_rev = Q_in / omega_0_fo_rev
        tau_v_fo_rev = Q_v / omega_0_fo_rev
        delta_lambda_therm_fo_rev = (lambda_0 / n_0) * dn_by_dT * photon_energy_fixed_omega_rev[i] * (1 / tau_linear) * R
        delta_lambda_kerr_fo_rev = (n_2 * c * photon_energy_fixed_omega_rev[i] * lambda_0) / (n_0 * n_g * V_kerr)



# Plotting
plt.figure(dpi=400)
for i in range(len(P_in)):
    plt.plot(c / (omega * 1e-6 / (2 * np.pi)), 1e-3 * photon_energy[i, :], label=f'$P_i$= {P_in[i]:0.3f} mW', linewidth=1.5)
plt.title('Spectral Response (All-Pass)')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Cavity Photon Energy [J]')
plt.xlim([min(c / (omega * 1e-6 / (2 * np.pi))), max(c / (omega * 1e-6 / (2 * np.pi)))])
plt.legend(loc='upper right', fontsize=6)
plt.show()

plt.figure(dpi=400)
for i in range(len(P_in)):
    plt.plot(x_n_itr, 1e-3 * photon_energy_evo[i,round(len(omega)/2),:], label=f'$P_i$= {P_in[i]:0.3f} mW', linewidth=1.5)
plt.title('$\lambda$ = '+ np.array2string(2*np.pi*1e9*c/omega[round(len(omega)/2)], formatter={'float_kind':lambda x: "%.3f" % x}) + ' nm')
plt.xlabel('Number of iterations')
plt.ylabel('Cavity Photon Energy [J]')
plt.xlim([min(x_n_itr), max(x_n_itr)])
plt.legend(loc='lower right', fontsize=6)
plt.show()

delta_omega_loc = len(omega)/32
print('round(len(omega)/2 + len(omega/4)):: ',round(len(omega)/2 + delta_omega_loc))
print('round(len(omega)/2 - len(omega/4)):: ',round(len(omega)/2 - delta_omega_loc))


plt.figure(dpi=400)
plt.plot(x_n_itr, 1e-3 * photon_energy_evo[4,round(len(omega)/2),:], label=f'$\lambda$= {2*np.pi*1e9*c/omega[round(len(omega)/2)]:0.3f} nm', linewidth=1.5)
plt.plot(x_n_itr, 1e-3 * photon_energy_evo[4,round(len(omega)/2 + delta_omega_loc),:], label=f'$\lambda$= {2*np.pi*1e9*c/omega[round(len(omega)/2 + delta_omega_loc)]:0.3f} nm', linewidth=1.5)
plt.plot(x_n_itr, 1e-3 * photon_energy_evo[4,round(len(omega)/2 - delta_omega_loc),:], label=f'$\lambda$= {2*np.pi*1e9*c/omega[round(len(omega)/2 - delta_omega_loc)]:0.3f} nm', linewidth=1.5)
plt.title('$P_i$ = '+ np.array2string(P_in[4], formatter={'float_kind':lambda x: "%.2f" % x}) + ' mW')
plt.xlabel('Number of iterations')
plt.ylabel('Cavity Photon Energy [J]')
plt.xlim([min(x_n_itr), max(x_n_itr)])
plt.legend(loc='upper right', fontsize=5)
plt.show()

plt.figure(dpi=400)
for i in range(len(P_in)):
    plt.plot(c / (omega * 1e-6 / (2 * np.pi)), photon_energy[i, :] / P_in[i], label=f'$P_i$= {P_in[i]:0.3f} mW', linewidth=1.5)
plt.title('Spectral Response (All-Pass)')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Normalized Cavity Photon Energy')
plt.xlim([min(c / (omega * 1e-6 / (2 * np.pi))), max(c / (omega * 1e-6 / (2 * np.pi)))])
plt.legend(loc='upper right', fontsize=6)
plt.show()

plt.figure(dpi=400)
plt.plot(P_in_fo, photon_energy_fixed_omega, linewidth=1.5)
plt.plot(P_in_fo_rev, photon_energy_fixed_omega_rev, linewidth=1.5)
plt.title('Cavity Energy vs Input Power')
plt.xlabel('P$_i$ [mW]')
plt.ylabel('Cavity Photon Energy')
plt.xlim([min(P_in_fo), max(P_in_fo)])
# plt.legend(loc='upper right', fontsize=6)
plt.show()

