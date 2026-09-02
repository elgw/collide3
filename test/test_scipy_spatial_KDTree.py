import sys
import time
import numpy as np
from scipy.spatial import KDTree


def test_for_N(N):
    t_construct = 0
    t_scan = 0
    t_total = 0
    n_tries = 0
    while t_total < 1000.0:
        X = np.random.random((N, 3))*2 -1
        X = X.astype(np.float64)

        t0 = time.time()
        tree = KDTree(X)
        t1 = time.time()
        t_construct += (t1-t0)*1000

        radius = 2 / N**(1/3)
        t2 = time.time()
        #for ii in range(0, X.shape[0]):
        #    indices = tree.query_ball_point(X[ii, :], radius)
        results = tree.query(X, k=12, p=2, distance_upper_bound=radius)
        t3 = time.time()
        t_scan += (t3-t2)*1000
        t_total += ( (t3-t2) + (t1-t0) )*1000
        n_tries += 1

    t_construct = t_construct/n_tries
    t_scan = t_scan/n_tries
    t_total = t_total/n_tries
    print(f'| KDTree | {N:,} | {t_construct:,.2f} | {t_scan:,.2f} | {t_total:,.2f} |')

if __name__ == '__main__':
    N = 2*8388608

    if len(sys.argv) > 1:
        N = int(sys.argv[1])


    print('| method   |       N | t_construct [ms] | t_scan [ms] | t_total [ms] |')
    print('|----------|--------:|-----------------:|------------:|-------------:|')
    n = 1024
    while n <= N:
        test_for_N(n)
        n*=2
