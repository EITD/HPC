import matplotlib.pyplot as plt
import numpy as np

sizes = np.array([8, 16, 32, 64, 128, 256])
times = np.array([2.790698, 1.401615, 1.079528, 0.542800, 0.319607, 0.211083])

plt.figure(figsize=(10, 7))
plt.plot(sizes, times, marker='o', linestyle='-')

plt.xlabel('Number of Processes')
plt.ylabel('Execution Time (seconds)')

plt.title('Pi Benchmark')
plt.xticks(sizes)
plt.savefig('pi.png')
