import numpy as np
import matplotlib.pyplot as plt
"""Produces plot of local energy and variance of
first 5 files in Bulk2 directory results"""

N_files = 1
baseFilename = "../../Bulk2/datPartI2NP2D3MC"
F = 19
fig1 = plt.figure(figsize=(11.5,8))

plt.style.use("seaborn-notebook") #Figure for local energy
for i in range(N_files):
    filename = baseFilename + str(i+6)
    data = np.loadtxt(filename, unpack=True)
    plt.plot(data[0], data[5], label=r"$N_{mc} = 10^%d$" % np.log10(data[-1][0]))
#    plt.plot(data[0], data[6], label=r"Closed-Form  $E_L$")
plt.xlabel(r"$\alpha$", fontsize=F+5)
plt.ylabel(r"$\langle E_L \rangle$", fontsize=F+5)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.axhline(np.min(data[6]),0.0 ,1.9, linestyle="--", label=r"$Closed-Form \;(E_L)_{min}=%g$"%np.min(data[6]))
plt.axhline(np.min(data[5]),0.0 ,1.9, linestyle="--", color="red", label=r"$Numerical \;(E_L)_{min}=%g$"%np.min(data[5]))
optimalAlpha = data[0][np.argmin( data[5] )]
plt.axvline(optimalAlpha , 0.0, 1.9, linestyle="--", color="red", label=r"Optimal $\alpha =%g$"% optimalAlpha)
plt.axis([0.2, 1.9, 3, 5])
plt.legend(fontsize=F)

fig2 = plt.figure(figsize=(10.5,8)) #Figure for variance
plt.style.use("seaborn-notebook")

for i in range(N_files):
    filename = baseFilename + str(i+6)
    data = np.loadtxt(filename, unpack=True)
    plt.plot(data[0], data[7], label=r"$N_{mc} = 10^%d$" % np.log10(data[-1][0]))
plt.axhline(np.min(data[7]),0.0 ,1.9, linestyle="--", color="red", label=r"$\sigma^2_{min}=%g$"%np.min(data[7]))
plt.xlabel(r"$\alpha$", fontsize=F+5)
plt.ylabel(r" $\sigma^2 $", fontsize=F+5)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.axhline(0.0,0 ,1.9, linestyle="--", color="black", label=r"$\sigma^2=0$")
plt.axis([0.1, 1.9, -0.5, 3])
plt.legend(fontsize=F)

plt.show()

fig1.savefig("../../../project/figsPartI/IStabilityPertE.png")
fig2.savefig("../../../project/figsPartI/IStabilityPertSigma.png")
