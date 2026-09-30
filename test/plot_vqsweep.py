import os

import matplotlib.pyplot as plt
import numpy as np

"""
./test_collide3_f32 --vqsweep --npoint 100000 > vqsweep_100000.csv
./test_collide3_f32 --vqsweep --npoint 1000000 > vqsweep_1000000.csv

"""

T5 = np.genfromtxt('vqsweep_100000.csv', dtype=None,  delimiter=',', names=True)
T6 = np.genfromtxt('vqsweep_1000000.csv', dtype=None,  delimiter=',', names=True)

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[11, 5])
for (ii, T, ax1) in zip([0, 1], [T5, T6], ax):
    ax1.plot(T['Vq'], T['t_total_ms'], 'k')
    ax1.set_xlabel('Volume Quotient [-]')
    if ii == 0:
        ax1.set_ylabel('t_total [ms]', color='k')
    ax1.set_ylim([0, 1.1*np.max(T['t_total_ms'])])
    ax2 = ax1.twinx()
    y = T['collisions']/(T['t_total_ms']/1000)
    ax2.plot(T['Vq'], y, 'g')
    if ii == 1:
        ax2.set_ylabel('Collisions per second', color='g')
    ax2.set_ylim([0, 1.01*np.max(y)])

ax[0].set_title('N=100,000')
ax[1].set_title('N=1,000,000')

if not os.path.exists('test/vqsweep.png'):
    print(f'Writing to test/vqsweep.png')
    plt.savefig('test/vqsweep.png')

plt.show()
