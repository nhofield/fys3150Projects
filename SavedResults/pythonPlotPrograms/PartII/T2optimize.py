import numpy as np
import matplotlib.pyplot as plt
"""Program for Part II, reads produced file, alpha and beta are
interchanged for finding optimal parameters"""

file = "dataNP2D3MC6"   #Specify file here

data = np.loadtxt(file, unpack=True)
par = "alpha"
EL = data[5]
varParam = data[0]  #0 for alpha, 1 for beta
F = 19

fig1 = plt.figure(figsize=(12.5,8))
plt.style.use("seaborn-notebook") #Figure for local energy
plt.plot(varParam, EL, label=r"$N_{mc}=10^6, \beta = 0.277$"  )
fitdegree = 10
p = np.polyfit(varParam, EL, fitdegree)
fit = 0
for i in range(fitdegree+1):
    fit = fit + p[i]*varParam**(fitdegree-i)

plt.plot(varParam, fit, "--", label="%d-Degree LSQ-Fit" % fitdegree)
plt.axhline(np.min(fit), 0.0, 3.0, linestyle="--", label=r"$\langle E_L \rangle_{min}=%g$"%np.min(fit))
plt.axvline(varParam[np.argmin(fit)], 0.0, 3.0, linestyle="--", label=r"$Optimal\; \%s=%g$" % (par, varParam[np.argmin(fit)]))

plt.xlabel(r"$\%s$" % par, fontsize=F+5)
plt.ylabel(r"$\langle E_L \rangle$", fontsize=F+5)
"""
plt.axhline(np.min(EL), 0.0, 3.0, linestyle="--", label=r"$\langle E_L \rangle_{min}=%g$"%np.min(EL))
plt.axvline(varParam[np.argmin(EL)], 0.0, 3.0, linestyle="--", label=r"$Optimal\; \%s=%g$" % (par, varParam[np.argmin(EL)]))
"""
plt.tick_params(axis='both', which='major', labelsize=F)
plt.legend(fontsize=F)

plt.show()
