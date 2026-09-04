import os
import glob

import matplotlib.pyplot as plt
import numpy as np

"""
# Generate tables with
stress -t 20 -c 1
./test_collide3_f32 --timings --timeout 2000 --backend bf > timings_bf.csv
./test_collide3_f32 --timings --timeout 2000 --backend lowmem > timings_lowmem.csv
./test_collide3_f32 --timings --timeout 2000 --backend spatial > timings_spatial.csv
"""

files = glob.glob('timings_*.csv')

names = [x.replace('timings_', '') for x in files]
names = [x.replace('.csv', '') for x in names]

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[12, 5])

for fname in files:
    T = np.genfromtxt(fname, dtype=None,  delimiter=',', names=True)

    #ax.scatter(T['N'], T['t_ms'], 10, 'k')
    x = T['N']
    y = T['t_ms']
    ax[0].plot(x, y)
    x = T['N']
    y = 1e6*T['t_ms']/T['N']
    ax[1].plot(x, y);

for a in [0, 1]:
    ax[a].legend(names)
    ax[a].set_xlabel('number of points')
    ax[a].set_xscale('log')
    ax[a].set_yscale('log')

ax[0].set_ylabel('total time [ms]')
ax[1].set_ylabel('time per point [ns]')
ax[1].set_ylim([50, 500])
if not os.path.exists('test/timings.png'):
    plt.savefig('test/timings.png')

plt.show()
