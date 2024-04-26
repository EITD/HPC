import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sizes = np.array([16, 25, 36, 49, 64, 100, 144, 196, 225, 256]) 
times = np.array([1.160940, 0.744931, 0.518897, 0.383246, 0.293978, 0.190490, 0.133844, 0.105423, 0.096787, 0.088070])

def performance_model(p, a, b):
    return a/p + b/np.sqrt(p)

params, params_covariance = curve_fit(performance_model, sizes, times)

plt.figure(figsize=(10, 7))
plt.plot(sizes, times, marker='o', linestyle='-', label='Measured Data')
plt.plot(sizes, performance_model(sizes, *params), color='red', label='Fitted function: a/p + b/sqrt(p)')
plt.xlabel('Number of Processes')
plt.ylabel('Execution Time (seconds)')
plt.title('Performance Modeling of Distributed Matrix Multiplication')
plt.xticks(sizes)
plt.legend()
plt.savefig('model.png')

print("Fitted parameters: a =", params[0], ", b =", params[1])
