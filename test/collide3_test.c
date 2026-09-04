#include <assert.h>
#include <getopt.h>
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

typedef struct {
    collide3_backend be1;
    collide3_backend be2;
    u32 n; // number of points, interpretation depends on method
    int markdown_timings;
    int csv_timings;
    int chimera;
    int validate;
    int verbose;
    int timeout;
} config;

// Return a beads radius so that
// V_beads / V_domain = Vq
// V_domain = 2^3 = 8
// V_beads = n * (4/3)*pi*r^3
double get_bead_radius(double Vq, double n)
{
    double r = cbrt(Vq*6.0/(M_PI*(double) n));
    assert(fabs(Vq - ((double) n*4.0/3.0*M_PI*pow(r, 3))/8.0) < 1e-6);
    return r;
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

collide3_backend parse_backend(const char * str)
{
    if(strcmp(str, "lowmem") == 0)
    {
        return be_spatial_lowmem;
    }
    if(strcmp(str, "spatial") == 0)
    {
        return be_spatial;
    }
    if(strcmp(str, "bf") == 0)
    {
        return be_brute_force;
    }

    printf("Could not parse '%s' as a backend\n", str);
    return be_auto;
}


// Return a number in [-1, 1]
static fxx inbox(void)
{
    return 2*(0.5 -((fxx) rand() / (fxx) RAND_MAX));
}

// ||X-Y||^2
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
markdown_timings(config * conf)
{
    u32 min_tries = 3;
    u32 min_time = 100;
    printf("--- Will perform some timing experiments\n");
    printf("--- backend: %s\n", backend_s(conf->be1));
    printf("--- for each n, at least %u iterations are performed\n", min_tries);
    printf("--- and at least %u ms are used\n", min_time);
    printf("\n");
    printf("| method |    N | t_construct [ms] | t_scan [ms] | t_total [ms] |  mem [kb] |\n");
    printf("|  ----  | ---: |     ---:         | ---:        | ---:         | ---:     |\n");
    for(i32 n = 16; n < 16*10000000; n*=2)
    {
        fxx Vq = 0.1;
        fxx bead_radius = get_bead_radius(Vq, n);
        fxx detection_distance = 2.0*bead_radius;
        double ntries = 0;
        double t_total = 0;
        double t_scan = 0;
        double t_construct = 0;
        double mem_b = 0;
        while((ntries < min_tries) || (t_total < min_time))
        {
            fxx * X = random_points(n);
            collide3_info info = {};
            info.backend = conf->be1;
            fxx(collide3)(X, n, detection_distance,
                          &info,
                          NULL, NULL);
            free(X);
            t_total += info.t_total_ms;
            t_construct += info.t_create_ms;
            t_scan += info.t_scan_ms;
            mem_b = info.mem_alloc; // same each time
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
timings_csv(config * conf)
{
    double accumulation_time_ms = 100;
    printf("N, t_ms\n");
    for(i32 n = 16; n < 16*10000000; n*=1.1)
    {
        fxx radius = get_bead_radius(0.1, n);
        fxx detection_distance = 2.0*radius;
        fxx ntries = 0;
        double t_total = 0;
        while(t_total < accumulation_time_ms)
        {
            fxx * X = random_points(n);
            collide3_info info = {};
            info.backend = conf->be1;
            fxx(collide3)(X, n, detection_distance,
                          &info,
                          NULL, NULL);
            free(X);
            t_total += info.t_create_ms + info.t_scan_ms;
            ntries++;
        }
        printf("%u, %e\n", n, t_total / ntries);
        fflush(stdout);
        if(t_total > conf->timeout) {
            break;
        }
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
validate(config * conf)
{
    u32 ntest = 10000;
    u32 max_size = 1234;

    printf("--- Will perform %u tests comparing to brute force, n < %u\n", ntest, max_size);
    printf("--- Backend to be tested: %s\n", backend_s(conf->be1));
    for(u32 kk = 0; kk < ntest; kk++)
    {
        u32 n = rand() % max_size;
        printf("\r %6u/%6u test_vs_brute_force, N = %u ",
               kk+1, ntest, n);
        fflush(stdout);

        test_vs_brute_force(n, conf->be1);
    }
    printf("\n");
    printf("--- done\n\n");
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
gen_chimerax_file(config * conf)
{
    u32 n = conf->n;
    printf("Generating %u points\n", n);
    fxx * X = random_points(n);
    if(X == NULL){
        return EXIT_FAILURE;
    }
    fxx radius = 0.2/cbrt(n);
    printf("Bead radius = %.2e\n", radius);
    u32 * red = calloc(n, sizeof(u32));
    if(red == NULL) { goto fail1; }

    collide3_info info = {};
    if(fxx(collide3)(X, n, 2.0*radius, &info, mark_red_cb, (void*) red)) {
        goto fail2;
    }

    double t_total = info.t_create_ms + info.t_scan_ms;
    printf("Found %lu collisions in %.3f ms",
           info.n_collisions,
           t_total);
    if(info.n_collisions > 0){
        printf(" (%.0f per second)\n", 1000.0 / t_total);
    } else {
        printf("\n");
    }


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


static void usage(void){
    printf("Options:\n");
    printf("--backend str\n\t"
           "select backend. Valid options:\n\t"
           "- auto    -- which is the default\n\t"
           "- spatial -- spatial hashing\n\t"
           "- lowmem  -- low mem version of spatial hashing\n\t"
           "- bf      -- brute force\n");
    printf("--npoint n\n\t"
           "set number of points to use\n");
    printf("--validate\n\t"
           "Run some validation tests.\n");
    printf("--chimerax\n\t"
           "generate a .cmm file that can be opened with chimerax\n");
    printf("--timings\n\t"
           "perform some timings, output as csv to stdout\n");
    printf("--timeout ms\n\t"
           "abort the timings when more than this time has passed for\n\t"
           "a single measurement\n");
    printf("--verbose v\n\t"
           "set the verbosity level\n");
    printf("--help\n\t"
           "show this message\n");
    return;
}

config * parse_command_line(int argc, char ** argv)
{
    struct option longopts[] = {
        { "backend",  required_argument, NULL, '1' },
        { "markdown",  no_argument,       NULL, 'm' },
        { "validate", no_argument,        NULL, 'v' },
        { "chimera",  no_argument,        NULL, 'c' },
        { "timings",  no_argument,       NULL, 't' },
        { "npoint",   required_argument, NULL, 'n' },
        { "verbose",  required_argument, NULL, 'V' },
        { "help",     no_argument,       NULL, 'h' },
        {"timeout",   required_argument, NULL, 'T' },
        { NULL,           0,                 NULL,   0   }
    };

    config * conf = calloc(1, sizeof(config));
    conf->verbose = 1;
    conf->timeout = 1000;
    conf->n = 10000;

    char ch;
    while((ch = getopt_long(argc, argv,
                            "1:2:mvc",
                            longopts, NULL)) != -1) {
        switch(ch) {
        case '1':
            conf->be1 = parse_backend(optarg);
            break;
        case '2':
            conf->be1 = parse_backend(optarg);
            break;
        case 'm':
            conf->markdown_timings = 1;
            break;
        case 'v':
            conf->validate = 1;
            break;
        case 'c':
            conf->chimera = 1;
            break;
        case 'n':
            conf->n = atol(optarg);
            break;
        case 't':
            conf->csv_timings = 1;
            break;
        case 'h':
            usage();
            break;
        case 'V':
            conf->verbose = atoi(optarg);
            break;
        case 'T':
            conf->timeout = atol(optarg);
            break;
        }
    }
    return conf;
}

int main(int argc, char ** argv)
{
    setlocale(LC_NUMERIC, "");

    if(argc == 1){
        return example1();
    }

    config * conf = parse_command_line(argc, argv);

    if(conf->verbose > 1){
        printf("backend: %s\n", backend_s(conf->be1));
    }

    if(conf->markdown_timings){
        markdown_timings(conf);
    }

    if(conf->validate){
        validate(conf);
    }

    if(conf->csv_timings)
    {
        timings_csv(conf);
    }

    if(conf->chimera){
        return gen_chimerax_file(conf);
    }


    free(conf);
    return EXIT_SUCCESS;
}
