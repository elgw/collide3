#include <assert.h>
#include <locale.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>

#include "../collide3.h"

typedef uint32_t u32;
typedef int32_t i32;
typedef float f32;

#ifdef COLLIDE3_f64
typedef double fxx;
#define fxx(X) X##_f64
#else
typedef float fxx;
#define fxx(X) X##_f32
#endif


static fxx inbox(void)
{
    return 2*(0.5 -((fxx) rand() / (fxx) RAND_MAX));
}

static fxx eudist2(const fxx * X, const fxx * Y)
{
    return pow(X[0]-Y[0], 2) + pow(X[1]-Y[1], 2) + pow(X[2]-Y[2], 2);
}

typedef struct {
    u32 n;
    i32 * C;
} cbdata;

static fxx * random_points(u32 n)
{
    fxx * X = malloc(3*n*sizeof(fxx));
    assert(X != NULL);
    for(u32 kk = 0; kk < 3*n; kk++) {
        X[kk] = inbox();
    }
    return X;
}

void cb(u32 u, u32 v, __attribute__ ((unused)) double d2, void * _data)
{
    cbdata * data = (cbdata * ) _data;
    i32 * C = data->C;
    u32 n = data->n;
    if(u < v)
    {
        C[u + n*v]++;
    } else {
        C[v + n*u]++;
    }
    return;
}

void test_vs_brute_force(u32 n)
{
    if(n < 1){
        return;
    }
    fxx radius = 2.0/cbrt(n);
    fxx radius2 = radius*radius;

    fxx * X = random_points(n);

    // Find contacts by brute force
    i32 * C0 = calloc(n*n, sizeof(i32));
    assert(C0 != NULL);
    i32 ncont = 0;
    for(u32 kk = 0; kk < n; kk++) {
        for(u32 ll = kk+1; ll < n; ll++) {
            if(eudist2(X+3*kk, X+3*ll) < radius2) {
                C0[kk + n*ll] = 1;
                ncont++;
            }
        }
    }
    printf(" ncont = %u", ncont);
    fflush(stdout);

    // Find contacts by collide3 and note in C
    i32 * C = calloc(n*n, sizeof(i32));
    cbdata data;
    data.C = C;
    data.n = n;
    fxx(collide3)(X, n, radius,
                  NULL,
                  cb, &data);

    // Check if the contacts were found
    for(u32 kk = 0; kk < n; kk++) {
        for(u32 ll = kk+1; ll < n; ll++) {
            if(C0[kk + n*ll] != C[kk + n*ll]) {
                fxx dist = 0;
                if(C0[kk + n*ll] == 1)
                {
                    dist = sqrt(eudist2(X+3*kk, X+3*ll));
                    printf("\n");
                    printf("Missing contact between %u and %u\n", kk, ll);
                    printf("[%f, %f, %f], [%f, %f, %f], d=%f (r=%f)\n",
                           X[3*kk], X[3*kk+1], X[3*kk+2],
                           X[3*ll], X[3*ll+1], X[3*ll+2],
                           dist,
                           radius);
                } else {
                    dist = sqrt(eudist2(X+3*kk, X+3*ll));
                    printf("\n");
                    printf("False contact between %u and %u\n", kk, ll);
                    printf("[%f, %f, %f], [%f, %f, %f], d=%f (r=%f)\n",
                           X[3*kk], X[3*kk+1], X[3*kk+2],
                           X[3*ll], X[3*ll+1], X[3*ll+2],
                           dist,
                           radius);
                }
                // Only consider a failure if the distance
                // is not too close to the radius
                if( fabs(dist-radius) < radius / 1e-5){
                    continue;
                }
                printf("Test failed\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    free(C);
    free(C0);
    free(X);
}

void timings(void)
{
    printf("--- Will perform some timing experiments\n");
    printf("\n");
    printf("| method |    N | t_construct [ms] | t_scan [ms] | t_total [ms] |  mem [kb] |\n");
    printf("|  ----  | ---: |     ---:         | ---:        | ---:         | ---:     |\n");
    for(i32 n = 1024; n < 2*10000000; n*=2)
    {
        fxx radius = 2.0/cbrt(n);
        fxx ntries = 0;
        double t_total = 0;
        double t_scan = 0;
        double t_construct = 0;
        double mem_b = 0;
        while(t_total < 1000.0)
        {
            fxx * X = random_points(n);
            collide3_info info = {};
            fxx(collide3)(X, n, radius,
                          &info,
                          NULL, NULL);
            free(X);
            t_total += info.t_create_ms + info.t_scan_ms;
            t_construct += info.t_create_ms;
            t_scan += info.t_scan_ms;
            mem_b = info.mem_alloc;
            ntries++;
        }
        printf("| collide3 | %'u | %'.2f |  %'.2f | %'.2f | %'.2f |\n",
               n,
               t_construct / ntries,
               t_scan / ntries,
               t_total / ntries,
               mem_b / 1000.0);
    }
    printf("\n");
    printf("--done\n\n");
}

static void
example1(void)
{
    printf("-- Example1: just counting the number of collisions\n");
    u32 n = 4;
    fxx X[3*4] = {
        0, 0, 0,
        .1, .1, .1,
        .5, .5, .5,
        .7, .7, .7};
    collide3_info info = {};
    for(fxx radius = 0; radius < 1; radius+=0.1)
    {
        fxx(collide3)(X, n, radius,
                      &info,
                      NULL, NULL);
        printf("%lu collisions using detection radius=%.1f\n",
               info.n_collisions, radius);
    }
    printf("---done\n\n");
}


static void
validate(u32 ntest, u32 max_size)
{
    printf("--- Will perform %u tests comparing to brute force, n < %u\n", ntest, max_size);
    for(u32 kk = 0; kk < ntest; kk++)
    {
        u32 n = rand() % max_size;
        printf("\r %6u/%6u test_vs_brute_force, N = %u ",
               kk+1, ntest, n);
        fflush(stdout);

        test_vs_brute_force(n);
    }
    printf("\n");
    printf("--- done\n\n");
    return;
}

int main(int argc, char ** argv)
{
    setlocale(LC_NUMERIC, "");
    example1();

    validate(1000, 1000);

    timings();

    validate(1000000, 1234);

    printf("\n");

    return EXIT_SUCCESS;
}
