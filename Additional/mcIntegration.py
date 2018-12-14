import numpy as np
import matplotlib.pyplot as plt

#1, 3
f = lambda x: x**2 + 3
F_analytical = 18. - 10./3
a = 1
b = 3
M = 10000
def mcIntegrate(N, function, x_0, x_1):
    x_vals = np.random.uniform(a, b, N)
    width = x_1 - x_0
    return width*np.sum(f(x_vals))/N

N_ = range(M)
F = 19
plt.figure(figsize=(9, 8.5))
plt.subplot(2, 1, 1)
mcIntegrate = np.vectorize(mcIntegrate)
Integrations = mcIntegrate(N_, f, a, b)
plt.plot(N_ , Integrations, "r", label="MC-Integration")
plt.axhline(F_analytical, linestyle="--",label="Exact")
plt.xlabel("x", fontsize=F)
plt.ylabel("I", fontsize=F)
plt.tick_params(axis='both', which='major', labelsize=F)

plt.grid(True)
plt.legend(fontsize=F)

plt.subplot(2, 1, 2)
plt.plot(N_ , abs(Integrations-F_analytical)/F_analytical, "r", label="Error")
plt.xlabel("x", fontsize=F)
plt.ylabel("Error", fontsize=F)
plt.tick_params(axis='both', which='major', labelsize=F)

plt.grid(True)
plt.legend(fontsize=F)
plt.show()
