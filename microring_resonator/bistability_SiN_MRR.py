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
P_in = np.array([0.001, 0.01, 0.1, 1, 5, 10, 16])
tau_in = Q_in / omega_0
tau_v = Q_v / omega_0

omega = np.linspace(2 * np.pi * 193.57e12, 2 * np.pi * 193.25e12, 300)
omega_rev = np.linspace(2 * np.pi * 193.25e12, 2 * np.pi * 193.57e12, 300)
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

n_itr = 1 # No. of iterations for a fixed omega

# Nonlinear CMT model
for i in range(len(P_in)):
    for j in range(len(omega)):
        for z in range(n_itr):
            photon_energy[i, j] = ((1 / tau_in) * (P_in[i] / 2)) / ((omega[j] - omega_0) ** 2 + (((1 / tau_in) + (1 / tau_v) + (1 / tau_linear)) ** 2) / 4)
            
            omega_0 = 2 * np.pi * c / (lambda_0 + delta_lambda_therm + delta_lambda_kerr)
            tau_in = Q_in / omega_0
            tau_v = Q_v / omega_0
            delta_lambda_therm = (lambda_0 / n_0) * dn_by_dT * photon_energy[i, j] * (1 / tau_linear) * R
            delta_lambda_kerr = (n_2 * c * photon_energy[i, j] * lambda_0) / (n_0 * n_g * V_kerr)

# Cavity Energy vs Input Power
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


# Nonlinear CMT model (Transmission) [MM]
delta_lambda_therm = 0
delta_lambda_kerr = 0
omega_0 = 2 * np.pi * c / (lambda_0 + delta_lambda_therm + delta_lambda_kerr)
tau_in = Q_in / omega_0
tau_v = Q_v / omega_0

delta_lambda_therm_fo = 0
delta_lambda_kerr_fo = 0
omega_0_fo = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_fo + delta_lambda_kerr_fo)
tau_in_fo = tau_in
tau_v_fo = tau_v
transmission_fixed_omega = np.zeros((len(P_in_fo), 1))
output_power_fixed_omega = np.zeros((len(P_in_fo), 1))

delta_lambda_therm_fo_rev = 0
delta_lambda_kerr_fo_rev = 0
omega_0_fo_rev = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_fo_rev + delta_lambda_kerr_fo_rev)
tau_in_fo_rev = tau_in
tau_v_fo_rev = tau_v
transmission_fixed_omega_rev = np.zeros((len(P_in_fo), 1))
output_power_fixed_omega_rev = np.zeros((len(P_in_fo), 1))

transmission = np.zeros((len(P_in), len(omega)))
transmission_or = np.zeros((len(P_in), len(omega)))
photon_energy_or = np.zeros((len(P_in), len(omega)))

radius = 50e-6
L=2*np.pi*radius

for i in range(len(P_in)):
    for j in range(len(omega)):
        for z in range(n_itr):
            a = np.exp(-((omega_0*n_g)/(Q_v*c))*L/2)
            k = np.sqrt((omega_0*L*n_g)/(Q_in*c))
            r = np.sqrt(1-k**2)
            gamma = ((1 / tau_in) + (1 / tau_v) + (1 / tau_linear))/2
            transmission[i, j] = 1-(((1-a**2)*k**2)/((1-r*a)**2))*(1/(1+((omega[j]-omega_0)/(gamma))**2))
            
            omega_0 = 2 * np.pi * c / (lambda_0 + delta_lambda_therm + delta_lambda_kerr)
            tau_in = Q_in / omega_0
            tau_v = Q_v / omega_0
            delta_lambda_therm = (lambda_0 / n_0) * dn_by_dT * photon_energy[i, j] * (1 / tau_linear) * R
            delta_lambda_kerr = (n_2 * c * photon_energy[i, j] * lambda_0) / (n_0 * n_g * V_kerr)

delta_lambda_therm = 0
delta_lambda_kerr = 0
omega_0 = 2 * np.pi * c / (lambda_0 + delta_lambda_therm + delta_lambda_kerr)
tau_in = Q_in / omega_0
tau_v = Q_v / omega_0

for i in range(len(P_in)):
    for j in range(len(omega_rev)):
        for z in range(n_itr):
            a_or = np.exp(-((omega_0*n_g)/(Q_v*c))*L/2)
            k_or = np.sqrt((omega_0*L*n_g)/(Q_in*c))
            r_or = np.sqrt(1-k_or**2)
            gamma_or = ((1 / tau_in) + (1 / tau_v) + (1 / tau_linear))/2
            photon_energy_or[i, j] = ((1 / tau_in) * (P_in[i] / 2)) / ((omega_rev[j] - omega_0) ** 2 + (((1 / tau_in) + (1 / tau_v) + (1 / tau_linear)) ** 2) / 4)
            transmission_or[i, j] = 1-(((1-a_or**2)*k_or**2)/((1-r_or*a_or)**2))*(1/(1+((omega_rev[j]-omega_0)/(gamma_or))**2))
            
            omega_0 = 2 * np.pi * c / (lambda_0 + delta_lambda_therm + delta_lambda_kerr)
            tau_in = Q_in / omega_0
            tau_v = Q_v / omega_0
            delta_lambda_therm = (lambda_0 / n_0) * dn_by_dT * photon_energy_or[i, j] * (1 / tau_linear) * R
            delta_lambda_kerr = (n_2 * c * photon_energy_or[i, j] * lambda_0) / (n_0 * n_g * V_kerr)


for i in range(len(P_in_fo)):
    for z in range(n_itr):
        a_fo = np.exp(-((omega_0_fo*n_g)/(Q_v*c))*L/2)
        k_fo = np.sqrt((omega_0_fo*L*n_g)/(Q_in*c))
        r_fo = np.sqrt(1-k_fo**2)
        gamma_fo = (1 / tau_in_fo) + (1 / tau_v_fo) + (1 / tau_linear)
        transmission_fixed_omega[i] = 1-(((1-a_fo**2)*k_fo**2)/((1-r_fo*a_fo)**2))*(1/(1+((fixed_omega-omega_0_fo)/(gamma_fo))**2))
        output_power_fixed_omega[i] = P_in_fo[i]*transmission_fixed_omega[i]

        omega_0_fo = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_fo + delta_lambda_kerr_fo)
        tau_in_fo = Q_in / omega_0_fo
        tau_v_fo = Q_v / omega_0_fo
        delta_lambda_therm_fo = (lambda_0 / n_0) * dn_by_dT * photon_energy_fixed_omega[i] * (1 / tau_linear) * R
        delta_lambda_kerr_fo = (n_2 * c * photon_energy_fixed_omega[i] * lambda_0) / (n_0 * n_g * V_kerr)

photon_energy_fixed_omega_rev = np.zeros((len(P_in_fo_rev), 1))
for i in range(len(P_in_fo)):
    for z in range(n_itr):
        a_fo_rev = np.exp(-((omega_0_fo_rev*n_g)/(Q_v*c))*L/2)
        k_fo_rev = np.sqrt((omega_0_fo_rev*L*n_g)/(Q_in*c))
        r_fo_rev = np.sqrt(1-k_fo_rev**2)
        gamma_fo_rev = (1 / tau_in_fo_rev) + (1 / tau_v_fo_rev) + (1 / tau_linear)
        transmission_fixed_omega_rev[i] = 1-(((1-a_fo_rev**2)*k_fo_rev**2)/((1-r_fo_rev*a_fo_rev)**2))*(1/(1+((fixed_omega-omega_0_fo_rev)/(gamma_fo_rev))**2))
        output_power_fixed_omega_rev[i] = P_in_fo_rev[i]*transmission_fixed_omega_rev[i]
        photon_energy_fixed_omega_rev[i] = ((1 / tau_in_fo_rev) * (P_in_fo_rev[i] / 2)) / ((fixed_omega - omega_0_fo_rev) ** 2 + (((1 / tau_in_fo_rev) + (1 / tau_v_fo_rev) + (1 / tau_linear)) ** 2) / 4)

        omega_0_fo_rev = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_fo_rev + delta_lambda_kerr_fo_rev)
        tau_in_fo_rev = Q_in / omega_0_fo_rev
        tau_v_fo_rev = Q_v / omega_0_fo_rev
        delta_lambda_therm_fo_rev = (lambda_0 / n_0) * dn_by_dT * photon_energy_fixed_omega_rev[i] * (1 / tau_linear) * R
        delta_lambda_kerr_fo_rev = (n_2 * c * photon_energy_fixed_omega_rev[i] * lambda_0) / (n_0 * n_g * V_kerr)



# Nonlinear CMT model (Transmission) [Wang etal]
delta_lambda_therm_wang = 0
delta_lambda_kerr_wang = 0
omega_0_wang = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_wang + delta_lambda_kerr_wang)
tau_in_wang = Q_in / omega_0_wang
tau_v_wang = Q_v / omega_0_wang

delta_lambda_therm_wang_fo = 0
delta_lambda_kerr_wang_fo = 0
omega_0_wang_fo = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_wang_fo + delta_lambda_kerr_wang_fo)
tau_in_wang_fo = tau_in_wang
tau_v_wang_fo = tau_v_wang
transmission_wang_fixed_omega = np.zeros((len(P_in_fo), 1))
output_power_wang_fixed_omega = np.zeros((len(P_in_fo), 1))

delta_lambda_therm_wang_fo_rev = 0
delta_lambda_kerr_wang_fo_rev = 0
omega_0_wang_fo_rev = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_wang_fo_rev + delta_lambda_kerr_wang_fo_rev)
tau_in_wang_fo_rev = tau_in_wang
tau_v_wang_fo_rev = tau_v_wang
transmission_wang_fixed_omega_rev = np.zeros((len(P_in_fo), 1))
output_power_wang_fixed_omega_rev = np.zeros((len(P_in_fo), 1))

transmission_wang = np.zeros((len(P_in), len(omega)))
transmission_wang_or = np.zeros((len(P_in), len(omega_rev)))
# radius = 50e-6
# L=2*np.pi*radius
for i in range(len(P_in)):
    for j in range(len(omega_rev)):
        # a = np.exp(-((omega_0*n_g)/(Q_v*c))*L/2)
        # k = np.sqrt((omega_0*L*n_g)/(Q_in*c))
        # r = np.sqrt(1-k**2)
        # gamma = (1 / tau_in) + (1 / tau_v) + (1 / tau_linear)
        for z in range(n_itr):
            transmission_wang[i, j] = (((1/tau_in_wang)-(1/tau_v_wang)-(1/tau_linear))**2+(omega[j] - omega_0_wang)**2)/(((1/tau_in_wang)+(1/tau_v_wang)+(1/tau_linear))**2+(omega[j] - omega_0_wang)**2)
            
            omega_0_wang = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_wang + delta_lambda_kerr_wang)
            tau_in_wang = 2*Q_in / (omega_0_wang)
            tau_v_wang = 2*Q_v / (omega_0_wang)
            delta_lambda_therm_wang = (lambda_0 / n_0) * dn_by_dT * photon_energy[i, j] * (1 / tau_linear) * R
            delta_lambda_kerr_wang = (n_2 * c * photon_energy[i, j] * lambda_0) / (n_0 * n_g * V_kerr)

delta_lambda_therm_wang = 0
delta_lambda_kerr_wang = 0
omega_0_wang = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_wang + delta_lambda_kerr_wang)
tau_in_wang = Q_in / omega_0_wang
tau_v_wang = Q_v / omega_0_wang

for i in range(len(P_in)):
    for j in range(len(omega_rev)):
        # a = np.exp(-((omega_0*n_g)/(Q_v*c))*L/2)
        # k = np.sqrt((omega_0*L*n_g)/(Q_in*c))
        # r = np.sqrt(1-k**2)
        # gamma = (1 / tau_in) + (1 / tau_v) + (1 / tau_linear)
        for z in range(n_itr):
            transmission_wang_or[i, j] = (((1/tau_in_wang)-(1/tau_v_wang)-(1/tau_linear))**2+(omega_rev[j] - omega_0_wang)**2)/(((1/tau_in_wang)+(1/tau_v_wang)+(1/tau_linear))**2+(omega_rev[j] - omega_0_wang)**2)
            
            omega_0_wang = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_wang + delta_lambda_kerr_wang)
            tau_in_wang = 2*Q_in / (omega_0_wang)
            tau_v_wang = 2*Q_v / (omega_0_wang)
            delta_lambda_therm_wang = (lambda_0 / n_0) * dn_by_dT * photon_energy[i, j] * (1 / tau_linear) * R
            delta_lambda_kerr_wang = (n_2 * c * photon_energy[i, j] * lambda_0) / (n_0 * n_g * V_kerr)

delta_lambda_therm_wang_fo = 0
delta_lambda_kerr_wang_fo = 0
omega_0_wang_fo = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_wang_fo + delta_lambda_kerr_wang_fo)
tau_in_wang_fo = Q_in/omega_0_wang_fo
tau_v_wang_fo = Q_v/omega_0_wang_fo

for i in range(len(P_in_fo)):
    for z in range(n_itr):
        a_fo = np.exp(-((omega_0_wang_fo*n_g)/(Q_v*c))*L/2)
        k_fo = np.sqrt((omega_0_wang_fo*L*n_g)/(Q_in*c))
        r_fo = np.sqrt(1-k_fo**2)
        gamma_fo = (1 / tau_in_wang_fo) + (1 / tau_v_wang_fo) + (1 / tau_linear)
        transmission_wang_fixed_omega[i] = 1-(((1-a_fo**2)*k_fo**2)/((1-r_fo*a_fo)**2))*(1/(1+((fixed_omega-omega_0_wang_fo)/(gamma_fo))**2))
        output_power_wang_fixed_omega[i] = P_in_fo[i]*transmission_wang_fixed_omega[i]

        omega_0_wang_fo = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_fo + delta_lambda_kerr_fo)
        tau_in_wang_fo = Q_in / omega_0_wang_fo
        tau_v_wang_fo = Q_v / omega_0_wang_fo
        delta_lambda_wang_therm_fo = (lambda_0 / n_0) * dn_by_dT * photon_energy_fixed_omega[i] * (1 / tau_linear) * R
        delta_lambda_wang_kerr_fo = (n_2 * c * photon_energy_fixed_omega[i] * lambda_0) / (n_0 * n_g * V_kerr)


delta_lambda_therm_fo_rev = 0
delta_lambda_kerr_fo_rev = 0
omega_0_fo_rev = 2 * np.pi * c / (lambda_0 + delta_lambda_therm_fo_rev + delta_lambda_kerr_fo_rev)
tau_in_fo_rev = Q_in/omega_0_fo_rev
tau_v_fo_rev = Q_v/omega_0_fo_rev

for i in range(len(P_in_fo_rev)):
    for z in range(n_itr):
        a_fo_rev = np.exp(-((omega_0_fo_rev*n_g)/(Q_v*c))*L/2)
        k_fo_rev = np.sqrt((omega_0_fo_rev*L*n_g)/(Q_in*c))
        r_fo_rev = np.sqrt(1-k_fo_rev**2)
        gamma_fo_rev = (1 / tau_in_fo_rev) + (1 / tau_v_fo_rev) + (1 / tau_linear)
        transmission_wang_fixed_omega_rev[i] = 1-(((1-a_fo_rev**2)*k_fo_rev**2)/((1-r_fo_rev*a_fo_rev)**2))*(1/(1+((fixed_omega-omega_0_fo_rev)/(gamma_fo_rev))**2))
        output_power_wang_fixed_omega_rev[i] = P_in_fo_rev[i]*transmission_wang_fixed_omega_rev[i]

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
    plt.plot(c / (omega * 1e-6 / (2 * np.pi)), photon_energy[i, :] / P_in[i], label=f'$P_i$= {P_in[i]:0.3f} mW', linewidth=1.5)
plt.title('Spectral Response (All-Pass)')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Normalized Cavity Photon Energy')
plt.xlim([min(c / (omega * 1e-6 / (2 * np.pi))), max(c / (omega * 1e-6 / (2 * np.pi)))])
plt.legend(loc='upper right', fontsize=6)
plt.show()

plt.figure(dpi=400)
for i in range(len(P_in)):
    plt.plot(c / (omega * 1e-6 / (2 * np.pi)), transmission[i, :], label=f'$P_i$= {P_in[i]:0.3f} mW', linewidth=1.5)
plt.title('Spectral Response (All-Pass) [MM]')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Transmission')
plt.xlim([min(c / (omega * 1e-6 / (2 * np.pi))), max(c / (omega * 1e-6 / (2 * np.pi)))])
plt.legend(loc='lower right', fontsize=6)
plt.show()

plt.figure(dpi=400)
for i in range(len(P_in)):
    plt.plot(c / (omega * 1e-6 / (2 * np.pi)), transmission_wang[i, :], label=f'$P_i$= {P_in[i]:0.3f} mW', linewidth=1.5)
plt.title('Spectral Response (All-Pass) [Wang]')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Transmission')
plt.xlim([min(c / (omega * 1e-6 / (2 * np.pi))), max(c / (omega * 1e-6 / (2 * np.pi)))])
plt.legend(loc='lower right', fontsize=6)
plt.show()

plt.figure(dpi=400)
for i in range(len(P_in)):
    plt.plot(c / (omega_rev * 1e-6 / (2 * np.pi)), transmission_or[i, :], label=f'$P_i$= {P_in[i]:0.3f} mW', linewidth=1.5)
plt.title('Spectral Response (All-Pass) [omega_rev, MM]')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Transmission')
plt.xlim([min(c / (omega_rev * 1e-6 / (2 * np.pi))), max(c / (omega_rev * 1e-6 / (2 * np.pi)))])
plt.legend(loc='lower right', fontsize=6)
plt.show()

plt.figure(dpi=400)
for i in range(len(P_in)):
    plt.plot(c / (omega_rev * 1e-6 / (2 * np.pi)), transmission_wang_or[i, :], label=f'$P_i$= {P_in[i]:0.3f} mW', linewidth=1.5)
plt.title('Spectral Response (All-Pass) [omega_rev, Wang]')
plt.xlabel('Wavelength [$\mu$m]')
plt.ylabel('Transmission')
plt.xlim([min(c / (omega_rev * 1e-6 / (2 * np.pi))), max(c / (omega_rev * 1e-6 / (2 * np.pi)))])
plt.legend(loc='lower right', fontsize=6)
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


plt.figure(dpi=400)
plt.plot(P_in_fo, output_power_fixed_omega, linewidth=1.5)
plt.plot(P_in_fo_rev, output_power_fixed_omega_rev, linewidth=1.5)
plt.title('Output Power vs Input Power [MM]')
plt.xlabel('P$_i$ [mW]')
plt.ylabel('Output Power [mW]')
plt.xlim([min(P_in_fo), max(P_in_fo)])
# plt.legend(loc='upper right', fontsize=6)
plt.show()

plt.figure(dpi=400)
plt.plot(P_in_fo, output_power_wang_fixed_omega, linewidth=1.5)
plt.plot(P_in_fo_rev, output_power_wang_fixed_omega_rev, linewidth=1.5)
plt.title('Output Power vs Input Power [Wang]')
plt.xlabel('P$_i$ [mW]')
plt.ylabel('Output Power [mW]')
plt.xlim([min(P_in_fo), max(P_in_fo)])
# plt.legend(loc='upper right', fontsize=6)
plt.show()

