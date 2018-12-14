import numpy as np
import matplotlib.pyplot as plt
"""Produces plot of local energy and variance of
first 5 files in Bulk1 directory results"""

N_files = 5
baseFilename = "../../Bulk1/datPartI1NP2D3MC"
F = 19
fig1 = plt.figure(figsize=(11.5,8))

plt.style.use("seaborn-notebook") #Figure for local energy
for i in range(N_files):
    filename = baseFilename + str(i+1)
    data = np.loadtxt(filename, unpack=True)
    plt.plot(data[0], data[5], label=r"$N_{mc} = 10^%d$" % np.log10(data[-1][0]))
plt.xlabel(r"$\alpha$", fontsize=F+5)
plt.ylabel(r"$\langle E_L \rangle$", fontsize=F+5)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.axhline(3.0,0 ,1.9, linestyle="--", color="red", label=r"$\mathcal{E}_U=3.0$")
plt.axis([0.2, 1.8, -10, 10])
plt.legend(fontsize=F)

fig2 = plt.figure(figsize=(10.5,8)) #Figure for variance
plt.style.use("seaborn-notebook")

for i in range(N_files):
    filename = baseFilename + str(i+1)
    data = np.loadtxt(filename, unpack=True)
    plt.plot(data[0], data[7], label=r"$N_{mc} = 10^%d$" % np.log10(data[-1][0]))
plt.xlabel(r"$\alpha$", fontsize=F+5)
plt.ylabel(r" $\sigma^2 $", fontsize=F+5)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.axhline(0.0,0 ,1.9, linestyle="--", color="black", label=r"$\sigma^2=0$")
plt.axis([0.3, 1.7, -0.5, 8])
plt.legend(fontsize=F)

plt.show()

fig1.savefig("../../../project/figsPartI/IStabilityE.png")
fig2.savefig("../../../project/figsPartI/IStabilitySigma.png")
