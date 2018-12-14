import numpy as np
import matplotlib.pyplot as plt

alphaT1U = 1.0
alphaT1P = 0.9

fileT1U = "../../Bulk2/datOmegasT1P"
fileT1P = "../../Bulk2/datOmegasT1U"

omegas = np.array([0.01, 0.5, 1.0])

datT1U = np.loadtxt(fileT1U, unpack=True)
datT1P = np.loadtxt(fileT1P, unpack=True)
meandistT1U = datT1U[-3]
meandistT1P = datT1P[-3]

F = 19  #Fontsize
fig1 = plt.figure(figsize=(10.5,8)) #Figure for local energy

ydist = [16, 3, 2]

plt.style.use("seaborn-notebook")
plt.plot(omegas, meandistT1U, label=r"Unperturbed $\alpha=%g$"%alphaT1U )
plt.plot(omegas, meandistT1P, label=r"Perturbed $\alpha=%g$"%alphaT1P )
plt.scatter(omegas, meandistT1U)
plt.scatter(omegas, meandistT1P)

plt.text(omegas[0]+0.05, ydist[0]+0.95, "(%g, %g)" % (omegas[0], meandistT1U[0]), fontsize=F-2)
plt.text(omegas[0]+0.05, ydist[0]-0.05, "(%g, %g)" % (omegas[0], meandistT1P[0]), fontsize=F-2)
plt.text(omegas[1], ydist[1]+1, "(%g, %g)" % (omegas[1], meandistT1U[1]), fontsize=F-2)
plt.text(omegas[1], ydist[1], "(%g, %g)" % (omegas[1], meandistT1P[1]), fontsize=F-2)
plt.text(omegas[2]-0.15, ydist[2]+1.4, "(%g, %g)" % (omegas[2], meandistT1U[2]), fontsize=F-2)
plt.text(omegas[2]-0.15, ydist[2]+0.4, "(%g, %g)" % (omegas[2], meandistT1P[2]), fontsize=F-2)

plt.xlabel(r"$\omega$", fontsize=F+5)
plt.ylabel(r"$ r_{12} $", fontsize=F+5)
plt.tick_params(axis='both', which='major', labelsize=F)
plt.legend(fontsize=F)
plt.show()

fig1.savefig("../../../project/figsPartI/IOmegas.png")
