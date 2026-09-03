import matplotlib.pyplot as plt
import numpy as np

# ./test_collide3_f32 --timings2 > timings2.csv

T = np.zeros([2000, 2]);
with open('../timings2.csv', 'r') as fid:
    print(fid.readline()) # consume header
    for ii, line in enumerate(fid):
        parts = line.split(',')
        T[ii, 0] = float(parts[0])
        T[ii, 1] = float(parts[1])
    T = T[0:ii, :]

fig, ax = plt.subplots(nrows=1, ncols=1)
t_ns = T[:, 1]*1000000
ax.scatter(T[:, 0], t_ns/T[:,0], 10, 'k')
ax.plot(T[:, 0], t_ns/T[:,0], 'k')
ax.set_xlabel('number of points')
ax.set_ylabel('time per point [ns]')
ax.set_xscale('log')
plt.show()
