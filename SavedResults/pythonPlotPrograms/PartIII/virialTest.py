import numpy as np
import matplotlib.pyplot as plt

fileU = "../../Bulk4/datOmegasT2U"
fileP = "../../Bulk4/datOmegasT2P"

datU = np.loadtxt(fileU, unpack=True)
datP = np.loadtxt(fileP, unpack=True)

T_expectU = datU[3]
T_expectP = datP[3]
U_expectU = datU[4]
U_expectP = datP[4]
omegas = datU[2]
ratioU = T_expectU/U_expectU
ratioP = T_expectP/U_expectP
F = 19

fig1 = plt.figure(figsize=(11.5,8))
fitdist = 5
plt.style.use("seaborn-notebook") #Figure for local energy
plt.scatter(omegas, ratioU, label="Unperturbed")
plt.scatter(omegas, ratioP, label="Perturbed")
plt.axhline(0.5, linestyle="--", color="red")
omega_fit = omegas[fitdist:]
pP = np.polyfit(omega_fit, ratioP[fitdist:], 1)
pU = np.polyfit(omega_fit, ratioU[fitdist:], 1)
plt.plot(omega_fit, pP[0]*omega_fit + pP[1],"--",color="orange", label="LSQ-Fit")
plt.plot(omega_fit, pU[0]*omega_fit + pU[1],"--",color="lightblue", label="LSQ-Fit")
print "Slope of Perturbed", pP[0]
print "Slope of Unperturbed", pU[0]
plt.xlabel(r"$\omega$", fontsize=F+5)
plt.ylabel(r"$\langle T \rangle/\langle U \rangle$", fontsize=F+5, rotation=0)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.axis("equal")
plt.legend(fontsize=F)

r12U = datU[-3]
r12P = datP[-3]
fig2 = plt.figure(figsize=(11.5,8))
fitdist = 5
plt.style.use("seaborn-notebook") #Figure for local energy
plt.plot(omegas, r12U, label="Unperturbed")
plt.plot(omegas, r12P, label="Perturbed")
plt.xlabel(r"$\omega$", fontsize=F+5)
plt.ylabel(r"$r_{12}$", fontsize=F+5)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.legend(fontsize=F)

plt.show()
