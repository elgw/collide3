import sys
import time
import numpy as np
import math
from scipy.spatial import KDTree


def test_for_N(Vq, N):
    t_construct = 0
    t_scan = 0
    t_total = 0
    n_tries = 0
    while t_total < 1000.0 or n_tries < 3:
        X = np.random.random((N, 3))*2 -1
        X = X.astype(np.float64)

        t0 = time.time()
        tree = KDTree(X)
        t1 = time.time()
        t_construct += (t1-t0)*1000

        radius = ((Vq*6)/(math.pi*n))**(1/3)
        Vb = n*4/3*math.pi*radius**3
        assert(abs(Vq - Vb/8 ) < 1e-6)
        detection_distance = 2*radius
        t2 = time.time()
        #for ii in range(0, X.shape[0]):
        #    indices = tree.query_ball_point(X[ii, :], radius)
        results = tree.query(X, k=12, p=2, distance_upper_bound=detection_distance)
        t3 = time.time()
        t_scan += (t3-t2)*1000
        t_total += ( (t3-t2) + (t1-t0) )*1000
        n_tries += 1

    t_construct = t_construct/n_tries
    t_scan = t_scan/n_tries
    t_total = t_total/n_tries
    print(f'| KDTree | {N:,} | {t_construct:,.3f} | {t_scan:,.3f} | {t_total:,.3f} |')

if __name__ == '__main__':
    N = 2*8388608

    if len(sys.argv) > 1:
        N = int(sys.argv[1])


    print('| method   |       N | t_construct [ms] | t_scan [ms] | t_total [ms] |')
    print('|----------|--------:|-----------------:|------------:|-------------:|')
    n = 1024
    Vq = 0.1
    while n <= N:
        test_for_N(Vq, n)
        n*=2
