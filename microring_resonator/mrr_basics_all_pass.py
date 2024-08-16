## MRR simulation (All-pass filter)

import numpy as np
import matplotlib.pyplot as plt 
import scienceplots

# plt.style.use('seaborn-v0_8-muted')
#plt.style.use(['science','no-latex', 'grid'])


alpha_dB = 10.56           # propgation loss in dB/m
alpha = alpha_dB/4.343  # propagation loss in 1/m
radius = 113e-6         # Radius of the MRR
L = 2*np.pi*radius      # Circumference of the ring
Lc = 0                  # Coupling length
#k = np.sqrt([0.0043, 0.0025])         # field coupling 
k = np.sqrt([0.0109, 0.0067, 0.004, 0.0025, 0.0015])# #linspace(0.01,0.05,5) #
power_coupling = k**2
r = np.sqrt(1-power_coupling)
a = np.sqrt(np.exp(-alpha*(L)))
ng = 2.0873             # Group index

lamda = np.linspace(start=1.54965, stop=1.54969,num=1000001) # linspace(249666,249672,1001)# # simulation wavelength range in um
#neff = 1.820350376124900  - 0.173507629941756*(lamda-1.55)  + 0.001013914079182*(lamda-1.55).^2# #1.75um Weff
neff = 1.822441325575143  - 0.171279836013422*(lamda-1.55) +  0.001158534242567*((lamda-1.55)**2)# #1.8um Weff
lamda = lamda*1e-6
beta = 2*np.pi*np.divide(neff, lamda)

T = np.zeros((len(k), len(lamda)))
Px = np.zeros((len(k), len(lamda)))
T_dB = np.zeros((len(k), len(lamda)))
Px_dB = np.zeros((len(k), len(lamda)))
FSR = np.zeros((len(k), len(lamda)))
FWHM = np.zeros((len(k), len(lamda)))
Q = np.zeros((len(k), len(lamda)))
Finesse = np.zeros((len(k), len(lamda)))


for i in range(len(k)):
    T[i,:] = (r[i]**2 + a**2 - 2*r[i]*a*np.cos(beta*(L))) / (1 + (r[i]*a)**2 - 2*r[i]*a*np.cos(beta*(L)))# # transmission
    Px[i,:] = (k[i]*a)**2 /(1 + (r[i]*a)**2 - 2*r[i]*a*np.cos(beta*(L)))# # power enhancement factor 
    T_dB[i,:] = 10*np.log10(T[i,:])#
    Px_dB[i,:] = 10*np.log10(Px[i,:])#

    FSR[i,:] = lamda**2 / (ng*(Lc+L))#
    FWHM[i,:] = ((1-r[i]*a)*lamda**2)/(np.pi*ng*(Lc+L)*np.sqrt(r[i]*a))#
    Q[i,:] = lamda / FWHM[i,:]#
    Finesse[i,:] = FSR[i,:] / FWHM[i,:]#

plt.figure(dpi=400)
for i in range(len(k)):
    plt.plot(lamda*1e6,T[i,:],label = f"$k^2={k[i]**2 * 100:.2f}$%",linewidth=1.5)#
    plt.title("Spectral Response (All-Pass)")
    plt.xlabel("Wavelength [$\mu$m]")#
    plt.ylabel("Transmission")#
    plt.xlim([min(lamda*1e6),max(lamda*1e6)])#
    #pbaspect([3 1 1])#
    #plt.hold(True)#
plt.legend(loc="upper right",fontsize="7.5")#

plt.figure(dpi=400)#
for i in range(len(k)):
    plt.plot(lamda*1e6,T_dB[i,:],label = f"$k^2={k[i]**2 * 100:.2f}$%",linewidth=1.5)#
    plt.title("Spectral Response (All-Pass)")
    plt.xlabel("Wavelength [$\mu$m]")#
    plt.ylabel("Transmission [dB]")#
    plt.xlim([min(lamda*1e6),max(lamda*1e6)])#
    #plt.hold(True)#
plt.legend(loc="upper right",fontsize="7.5")#

plt.figure(dpi=400)#
for i in range(len(k)):
    plt.plot(lamda*1e6,Px[i,:],label = f"$k^2={k[i]**2 * 100:.2f}$%",linewidth=1.5)#
    plt.title("Power Enhancement (All-Pass)")
    plt.xlabel("Wavelength [$\mu$m]")#
    plt.ylabel("Power Enhancement")#
    plt.xlim([min(lamda*1e6),max(lamda*1e6)])#
    #plt.hold(True)#
plt.legend(loc="upper right",fontsize="7.5")#

plt.figure(dpi=400)#
for i in range(len(k)):
    plt.plot(lamda*1e6,Px_dB[i,:],label = f"$k^2={k[i]**2 * 100:.2f}$%",linewidth=1.5)#
    plt.title("Power Enhancement (All-Pass)")
    plt.xlabel("Wavelength [$\mu$m]")#
    plt.ylabel("Power Enhancement [dB]")#
    plt.xlim([min(lamda*1e6),max(lamda*1e6)])#
    #plt.hold(True)#
plt.legend(loc="upper right",fontsize="7.5")#


plt.figure(dpi=400)#
for i in range(len(k)):
    plt.plot(lamda*1e6,FSR[i,:]*1e9,label = f"$k^2={k[i]**2 * 100:.2f}$%",linewidth=1.5)#
    plt.title("Free Spectral Range (All-Pass)")
    plt.xlabel("Wavelength [$\mu$m]")#
    plt.ylabel("FSR [nm]")#
    plt.xlim([min(lamda*1e6),max(lamda*1e6)])#
    #pbaspect([3 1 1])#
plt.legend(loc="upper right",fontsize="7.5")


plt.figure(dpi=400)#
for i in range(len(k)):
    plt.plot(lamda*1e6,FWHM[i,:]*1e12,label = f"$k^2={k[i]**2 * 100:.2f}$%",linewidth=1.5)#
    plt.title("Full Width at Half Maximum (All-Pass)")
    plt.xlabel("Wavelength [$\mu$m]")#
    plt.ylabel("FWHM [pm]")#
    plt.xlim([min(lamda*1e6),max(lamda*1e6)])#
    #pbaspect([3 1 1])#
    #hold on#
plt.legend(loc="upper right",fontsize="7.5")#


plt.figure(dpi=400)#
for i in range(len(k)):
    plt.plot(lamda*1e6,Q[i,:],label = f"$k^2={k[i]**2 * 100:.2f}$%",linewidth=1.5)#
    plt.title("Quality Factor (All-Pass)")
    plt.xlabel("Wavelength [$\mu$m]")#
    plt.ylabel("Q factor")#
    plt.xlim([min(lamda*1e6),max(lamda*1e6)])#
    
plt.legend(loc="upper right",fontsize="7.5")#


plt.figure(dpi=400)#
for i in range(len(k)):
    plt.plot(lamda*1e6,Finesse[i,:],label = f"$k^2={k[i]**2 * 100:.2f}$%",linewidth=1.5)#
    plt.title("Finesse (All-Pass)")
    plt.xlabel("Wavelength [$\mu$m]")#
    plt.ylabel("Finesse")#
    plt.xlim([min(lamda*1e6),max(lamda*1e6)])#
    
plt.legend(loc="upper right",fontsize="7.5")#

plt.figure(dpi=400)
plt.plot(lamda*1e6,neff,linewidth=1.5)
plt.title("Effective index")
plt.xlabel("Wavelength [$\mu$m]")
plt.ylabel("Effective index")
plt.xlim([min(lamda*1e6),max(lamda*1e6)])
#plt.legend(fontsize="7.5")

