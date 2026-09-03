#include <assert.h>
#include <locale.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

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

static const char* geqs(fxx a, fxx b)
{
    if(a < b) {
        return "<-";
    }
    if(a > b) {
        return "->";
    }
    return "==";
}

static const char *
backend_s(collide3_backend be)
{
    switch(be) {
    case be_auto:
        return "auto";
    case be_spatial:
        return "spatial";
    case be_brute_force:
        return "bf";
    case be_spatial_lowmem:
        return "lowmem";
    }
    return "???";
}

static double
timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}


static fxx inbox(void)
{
    return 2*(0.5 -((fxx) rand() / (fxx) RAND_MAX));
}

static fxx eudist2(const fxx * X, const fxx * Y)
{
    return pow(X[0]-Y[0], 2)
        +  pow(X[1]-Y[1], 2)
        +  pow(X[2]-Y[2], 2);
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

void test_vs_brute_force(u32 n, collide3_backend be)
{
    if(n <= 1){
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
    collide3_info info = {};
    info.backend = be;
    fxx(collide3)(X, n, radius,
                  &info,
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

static int
timings(void)
{
    printf("--- Will perform some timing experiments\n");
    printf("\n");
    printf("| method |    N | t_construct [ms] | t_scan [ms] | t_total [ms] |  mem [kb] |\n");
    printf("|  ----  | ---: |     ---:         | ---:        | ---:         | ---:     |\n");
    for(i32 n = 16; n < 16*10000000; n*=2)
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
        printf("| collide3 | %'u | %'.3f | %'.3f | %'.3f | %'.0f |\n",
               n,
               t_construct / ntries,
               t_scan / ntries,
               t_total / ntries,
               mem_b / 1000.0);
    }
    printf("\n");
    printf("--done\n\n");
    return EXIT_SUCCESS;
}

static int
timings2(void)
{
    printf("N, t_ms\n");
    for(i32 n = 16; n < 16*10000000; n*=1.1)
    {
        fxx radius = 2.0/cbrt(n);
        fxx ntries = 0;
        double t_total = 0;
        while(t_total < 100.0)
        {
            fxx * X = random_points(n);
            collide3_info info = {};
            fxx(collide3)(X, n, radius,
                          &info,
                          NULL, NULL);
            free(X);
            t_total += info.t_create_ms + info.t_scan_ms;
            ntries++;
        }
        printf("%u, %e\n", n, t_total / ntries);
        fflush(stdout);
    }
    return EXIT_SUCCESS;
}


static int
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
    return EXIT_SUCCESS;
}


static int
validate(u32 ntest, u32 max_size, collide3_backend be)
{
    printf("--- Will perform %u tests comparing to brute force, n < %u\n", ntest, max_size);
    for(u32 kk = 0; kk < ntest; kk++)
    {
        u32 n = rand() % max_size;
        printf("\r %6u/%6u test_vs_brute_force, N = %u ",
               kk+1, ntest, n);
        fflush(stdout);

        test_vs_brute_force(n, be);
    }
    printf("\n");
    printf("--- done\n\n");
    return EXIT_SUCCESS;
}

// Find where brute force isn't the fastest method
static int
compare_backends(collide3_backend be1,
                 collide3_backend be2,
                 u32 from, u32 to, u32 step)
{
    printf("       N,  %s,  %s\n", backend_s(be1), backend_s(be2));
    for(u32 n = from; n < to; n += step)
    {
        double t_be1 = 0;
        double t_be2 = 0;
        double n_iter = 0;
        double dt = 0;
        struct timespec t0, t1;
        clock_gettime(CLOCK_REALTIME, &t0);

        while(dt < 0.1){
            fxx * X = random_points(n);
            fxx radius = 2.0/cbrt(n);

            collide3_info info1 = {};
            info1.backend = be1;
            if(fxx(collide3)(X, n, radius, &info1, NULL, NULL))
            {
                return EXIT_FAILURE;
            }

            collide3_info info2 = {};
            info2.backend = be2;
            if(fxx(collide3)(X, n, radius, &info2, NULL, NULL))
            {
                return EXIT_FAILURE;
            }

            t_be1 += info1.t_create_ms+info1.t_scan_ms;
            t_be2 += info2.t_create_ms+info2.t_scan_ms;
            n_iter++;
            clock_gettime(CLOCK_REALTIME, &t1);
            dt = timespec_diff(&t1, &t0);
            free(X);
        }
        t_be1 /= n_iter;
        t_be2 /= n_iter;

        fxx xf = 0;
        xf = t_be1 < t_be2 ? t_be2/t_be1 : t_be1/t_be2;

        printf("%6u %.2e %s %.2e %.1fX\n", n,
               t_be1,
               geqs(t_be1, t_be2),
               t_be2,
               xf);

    }
    return EXIT_SUCCESS;
}

void mark_red_cb(u32 u, u32 v, __attribute__((unused)) double d2, void * data)
{
    u32 * red = (u32*) data;
    red[u] = 1;
    red[v] = 1;
    return;
}

static int
gen_chimerax_file(u32 n)
{
    fxx * X = random_points(n);
    if(X == NULL){
        return EXIT_FAILURE;
    }
    fxx radius = 0.2/cbrt(n);
    u32 * red = calloc(n, sizeof(u32));
    if(red == NULL) { goto fail1; }

    collide3_info info = {};
    if(fxx(collide3)(X, n, 2.0*radius, &info, mark_red_cb, (void*) red)) {
        goto fail2;
    }

    double t_total = info.t_create_ms + info.t_scan_ms;
    printf("Found %lu collisions in %.3f ms (%.0f per second)\n",
           info.n_collisions,
           t_total,
           1000.0/t_total);

    FILE * fid = fopen("test.cmm", "w");
    if(fid == NULL) {
        goto fail2;
    }
    fprintf(fid, "<marker_set name=\"dump\">\n");
    for(u32 kk = 0; kk < n; kk++) {
        double r = 0.5;
        double g = 0.5;
        double b = 0.5;
        if(red[kk]){
            r = 1; g = 0; b = 0;;
        }
        fprintf(fid, "<marker id=\"%u\" x=\"%.3f\" y=\"%.3f\" z=\"%.3f\" r=\"%f\" g=\"%f\" b=\"%f\" radius=\"%f\" />\n",
                kk, // marker_id
                X[3*kk+0],X[3*kk+1],X[3*kk+2],
                r, g, b,
                radius);
    }
    fprintf(fid, "</marker_set>\n");
    fclose(fid);
    return EXIT_SUCCESS;

 fail2:
    free(X);
 fail1:
    free(red);
    return EXIT_FAILURE;

}

int main(int argc, char ** argv)
{
    setlocale(LC_NUMERIC, "");

    if(argc > 1) {
        if( strcmp(argv[1], "--timings") == 0) {
            return timings();
        }
        if( strcmp(argv[1], "--timings2") == 0) {
            return timings2();
        }
        if( strcmp(argv[1], "--validate") == 0) {
            validate(1000, 100, be_brute_force);
            validate(1000, 1234, be_spatial_lowmem);
            validate(1000000, 1234, be_auto);
            return EXIT_SUCCESS;
        }
        if(strcmp(argv[1], "--example1") == 0) {
            return example1();
        }
        if(strcmp(argv[1], "--backends") == 0) {
            compare_backends(be_spatial, be_spatial_lowmem, 100, 10000000, 100);
            compare_backends(be_brute_force, be_spatial, 2, 100, 1);
        }
        if(strcmp(argv[1], "--chimerax") == 0) {
            return gen_chimerax_file(10000);
        }
        printf("Options:\n");
        printf("--timings\n\t"
               "Time the code for some problem sizes.\n\t"
               "Will be displayed as a markdown table\n");
        printf("--validate\n\t"
               "Run some validation tests.\n");
        printf("--example1\n\t"
               "Run example1\n");
        printf("--chimerax\n\t"
               "generate a .cmm file that can be opened with chimerax");
        return EXIT_FAILURE;
    }

    printf("No option provided, running a quick validation\n");
    validate(1000, 1000, be_auto);
    printf("\n");

    return EXIT_SUCCESS;
}
