import numpy as np
import matplotlib.pyplot as plt
"""Produces plot of local energy and variance of
last 4 files in Bulk1 directory results and minimum
values of local energy with associated alpha and N_mc"""

N_files = 4
alphas = np.zeros(N_files)
EL_min = np.zeros_like(alphas)

baseFilename = "../../Bulk1/datPartI1NP2D3MC"
F = 19  #Fontsize
fig1 = plt.figure(figsize=(10.5,8)) #Figure for local energy

plt.style.use("seaborn-notebook")
for i in range(N_files):
    filename = baseFilename + str(i+4)
    data = np.loadtxt(filename, unpack=True)
    plt.plot(data[0], data[5], label=r"$N_{mc} = 10^%d$" % np.log10(data[-1][0]))
    minIndex = np.argmin(data[5])
    alphas[i] = data[0][minIndex]
    EL_min[i] = data[5][minIndex]
plt.xlabel(r"$\alpha$", fontsize=F+5)
plt.ylabel(r"$ \langle E_L \rangle $", fontsize=F+5)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.axhline(3.0,0 ,1.9, linestyle="--", color="black", label=r"$\mathcal{E}=3.0$")
plt.axis([0.4, 1.6, 2.95, 3.2])
plt.legend(fontsize=F)

fig2 = plt.figure(figsize=(10.5,8)) #Figure for local energy

plt.style.use("seaborn-notebook")
for i in range(N_files):
    filename = baseFilename + str(i+4)
    data = np.loadtxt(filename, unpack=True)
    plt.plot(data[0], data[7], label=r"$N_{mc} = 10^%d$" % np.log10(data[-1][0]))
plt.xlabel(r"$\alpha$", fontsize=F+5)
plt.ylabel(r"$ \sigma^2 $", fontsize=F+5)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.axhline(0.0,0 ,1.9, linestyle="--", color="black", label=r"$\sigma^2=0$")
plt.axis([0.1, 1.9, -0.1, 2.0])
plt.legend(fontsize=F)
plt.show()

for i in range(4):
    print "alpha=%.8f" % alphas[i], "log10(N_mc)=", i+4, "EL_min= %.8f" % EL_min[i]

fig1.savefig("../../../project/figsPartI/IStabilityThermalizedE.png")
fig2.savefig("../../../project/figsPartI/IStabilityThermalizedSigma.png")
