# collide3 v1.0.2

<img src="test/screenshot.png" width="100%" />

A collisions detector for points/spheres in $`[-1,1]^3`$ with equal
radii. Used in [chromflock](https://www.github.com/elgw/chromflock).

The implementation combines spatial hashing/spatial partitioning (each
bin corresponds to a 3D box) with [counting
sort](https://en.wikipedia.org/wiki/Counting_sort) to generate a hash
table with a load factor of 100%. Good performance can only be
expected when the points are somewhat uniformly distributed.

The "demo" image above was generated with
``` shell
make ; ./test_collide3_f32 --chimera --npoint 10000; chimerax test.cmm
```

## API and Usage

See `collide3.h` for the latest version and usage notes. In essence:
the interface asks for a list of points and a callback function to
be used on each collision. See also the file `test/collide3_test.c`
for examples.

``` c
int
collide3_f32(const float * points,
         uint32_t n_point,
         float collision_distance,
         collide3_info * info,
         collide3_cb cb, void * cb_data);
```

## Performance indicators

This test use randomly distributed points in $`[-1,1]^3`$.  The bead
radius $`r`$ was set so that the volume quotient
$`V_{beads}/V_{domain}`$ is $`0.1`$, the collision distance
(between the beads centers) as $`2r`$. Smaller beads will of course
result in faster run-times and vice versa.

For 32-bit floating points:

| method   |          N | t_construct [ms] | t_scan [ms] | t_total [ms] |  mem [kb] |
|----------|-----------:|-----------------:|------------:|-------------:|----------:|
| collide3 |      1,024 |            0.008 |       0.068 |        0.076 |        17 |
| collide3 |      2,048 |            0.016 |       0.136 |        0.153 |        34 |
| collide3 |      4,096 |            0.032 |       0.270 |        0.303 |        68 |
| collide3 |      8,192 |            0.105 |       0.515 |        0.620 |       136 |
| collide3 |     16,384 |            0.217 |       0.975 |        1.192 |       273 |
| collide3 |     32,768 |            0.447 |       1.934 |        2.382 |       548 |
| collide3 |     65,536 |            0.602 |       4.061 |        4.663 |     1,097 |
| collide3 |    131,072 |            1.428 |       8.688 |       10.116 |     2,195 |
| collide3 |    262,144 |            3.104 |      15.503 |       18.606 |     4,397 |
| collide3 |    524,288 |            8.699 |      31.375 |       40.074 |     8,804 |
| collide3 |  1,048,576 |           22.380 |      62.443 |       84.823 |    17,599 |
| collide3 |  2,097,152 |           54.930 |     127.131 |      182.061 |    35,175 |
| collide3 |  4,194,304 |          121.150 |     253.942 |      375.092 |    70,431 |
| collide3 |  8,388,608 |          311.449 |     513.884 |      825.332 |   140,790 |
| collide3 | 16,777,216 |          902.530 |   1,108.272 |    2,010.802 |   281,667 |
| collide3 | 33,554,432 |        2,413.126 |   2,303.453 |    4,716.579 |   563,450 |
| collide3 | 67,108,864 |        4,982.988 |   4,331.374 |    9,314.363 | 1,126,990 |

By default `collide3` use brute force (`be_brute_force`)
for small problems and switches to `be_spatial` when there are 60
points or more. That is based on the following timings:

<img src="test/timings.png" width="100%" />


<details><summary>Results 64-bit floating points</summary>

| method   |          N | t_construct [ms] | t_scan [ms] | t_total [ms] |  mem [kb] |
|----------|-----------:|-----------------:|------------:|-------------:|----------:|
| collide3 |      1,024 |            0.006 |       0.058 |        0.064 |        33 |
| collide3 |      2,048 |            0.013 |       0.115 |        0.127 |        67 |
| collide3 |      4,096 |            0.065 |       0.221 |        0.287 |       134 |
| collide3 |      8,192 |            0.141 |       0.422 |        0.563 |       267 |
| collide3 |     16,384 |            0.127 |       1.057 |        1.184 |       535 |
| collide3 |     32,768 |            0.292 |       2.249 |        2.541 |     1,072 |
| collide3 |     65,536 |            0.654 |       4.159 |        4.812 |     2,146 |
| collide3 |    131,072 |            1.530 |       8.292 |        9.822 |     4,292 |
| collide3 |    262,144 |            4.816 |      16.712 |       21.529 |     8,591 |
| collide3 |    524,288 |           14.741 |      34.657 |       49.398 |    17,193 |
| collide3 |  1,048,576 |           38.532 |      71.617 |      110.150 |    34,376 |
| collide3 |  2,097,152 |           75.470 |     123.444 |      198.914 |    68,730 |
| collide3 |  4,194,304 |          158.516 |     236.271 |      394.787 |   137,540 |
| collide3 |  8,388,608 |          392.250 |     481.762 |      874.012 |   275,008 |
| collide3 | 16,777,216 |        1,051.769 |   1,091.246 |    2,143.015 |   550,103 |
| collide3 | 33,554,432 |        2,349.243 |   2,227.254 |    4,576.497 | 1,100,321 |
| collide3 | 67,108,864 |        5,461.413 |   4,066.265 |    9,527.677 | 2,200,732 |

</details>

## Comparisons


A [k-d tree](https://en.wikipedia.org/wiki/K-d_tree) can of course
also be used for the same problem (having much broader applicability).

Just to have some reference point I've run
[`scipy.spatial.KDTree`](https://docs.scipy.org/doc/scipy/reference/spatial.html)
with similar input as above (see `test/test_scipy_spatial_KDTree.py`),
which boils down to:

``` Python
tree = KDTree(X)
results = tree.query(X, k=2, p=2, distance_upper_bound=detection_distance)
```

Results:

| method |         N | t_construct [ms] | t_scan [ms] | t_total [ms] | mem [kb] |
|--------|----------:|-----------------:|------------:|-------------:|---------:|
| KDTree |     1,024 |            0.224 |       0.590 |        0.814 |          |
| KDTree |     2,048 |            0.428 |       1.135 |        1.564 |          |
| KDTree |     4,096 |            0.896 |       2.502 |        3.398 |          |
| KDTree |     8,192 |            1.734 |       4.830 |        6.564 |          |
| KDTree |    16,384 |            3.666 |      10.389 |       14.055 |          |
| KDTree |    32,768 |            8.218 |      23.239 |       31.457 |          |
| KDTree |    65,536 |           20.785 |      53.978 |       74.762 |          |
| KDTree |   131,072 |           42.853 |     110.646 |      153.500 |          |
| KDTree |   262,144 |           89.353 |     264.876 |      354.229 |          |
| KDTree |   524,288 |          183.615 |     839.618 |    1,023.233 |          |
| KDTree | 1,048,576 |          450.215 |   2,031.953 |    2,482.168 |          |
| KDTree | 2,097,152 |        1,013.275 |   4,512.824 |    5,526.099 |          |
| KDTree | 4,194,304 |        2,178.211 |   9,514.155 |   11,692.366 |          |
| KDTree | 8,388,608 |        4,995.838 |  20,391.409 |   25,387.246 |          |


Changing the data type of the input array from `np.float64` to
`np.float32` did not seem alter the timings.

## Notes/ Relevant links / See also

- The tests were performed with an AMD Ryzen 7 3700X 8-Core Processor using GCC 13.3.0.

- The visualization above was made with [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/).

- [`scipy.spatial.KDTree`](https://docs.scipy.org/doc/scipy/reference/spatial.html)

- I'd be happy to list your alternative here.
