import matplotlib.pyplot as plt
import numpy as np

sizes = np.array([16, 25, 36, 49, 64, 100, 144, 196, 225, 256])
times = np.array([1.160940, 0.744931, 0.518897, 0.383246, 0.293978, 0.190490, 0.133844, 0.105423, 0.096787, 0.088070])

plt.figure(figsize=(10, 7))
plt.plot(sizes, times, marker='o', linestyle='-')

plt.xlabel('Number of Processes')
plt.ylabel('Execution Time (seconds)')

plt.title('Fox Benchmark')
plt.xticks(sizes)
plt.savefig('fox.png')
