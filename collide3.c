// A collision detector for points in [-1, 1]^3
// Erik Wernersson 2026
//
// Version history:
//
// 1.0.4
// - Minor changes.
//
// 1.0.3
// - New: Added a low mem option which does not copy the input points
//
// 1.0.2
//
//- Improvement: Brute force scanning is faster for small
//  problems. The method switches brute force scanning when fewer than
//  64 points are provided.
//
// 1.0.1
//
// - Fix: The callback returns the distance between the points (not
//   the detection distance).

#include <assert.h>
#include <math.h> // for cbrt only
#include <stdio.h> // not needed
#include <time.h>
#include "collide3.h"

typedef uint8_t u8;
typedef int32_t i32;
typedef uint32_t u32;
typedef uint64_t u64;

#ifdef COLLIDE3_f64
typedef double fxx;
#define fxx(X) X##_f64
#else
typedef float fxx;
#define fxx(X) X##_f32
#endif

#define C3_SQUARE(x) ((x)*(x))

// hash table entry
typedef struct {
    fxx X[3]; // coordinate
    u32 idx; // original index
} entry;

// Euclidean distance between two 3D-vectors to the power of two
// i.e., ||A-B||^2
static fxx
eudist2(const fxx * A, const fxx * B)
{
    return C3_SQUARE(A[0]-B[0]) + C3_SQUARE(A[1]-B[1]) + C3_SQUARE(A[2]-B[2]);
}

static double
timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

// Hash function for a single coordinate
static u32
hash_coord(const i32 nDiv, const fxx x)
{
    const i32 bin = (i32) ((x+1.0)/2.0 * nDiv);
    if(bin <= 0)     { return 0; }
    if(bin >= nDiv-1){ return nDiv-1; }
    return bin;
}

// Hash function for a 3D point
static u32
hash(const i32 nDiv, const fxx * X)
{
    return hash_coord(nDiv, X[0]) +
        nDiv*hash_coord(nDiv, X[1]) +
        nDiv*nDiv*hash_coord(nDiv, X[2]);
}

// The default dummy callback used when no callback function is
// passed.
static void
dummy_cb(__attribute__((unused)) u32 u,
         __attribute__((unused)) u32 v,
         __attribute__((unused)) double d2,
         __attribute__((unused)) void * data)
{
    return;
}

// Used for small problem sizes
static int
collide3_brute_force(const fxx * D, const u32 N, const fxx d,
                     collide3_info * info,
                     collide3_cb cb, void * cb_data)
{
    if(cb == NULL) {
        cb = dummy_cb;
        cb_data = NULL;
    }

    double d2 = d*d;
    struct timespec t0, t1;
    if(info) {
        clock_gettime(CLOCK_REALTIME, &t0);
    }

    u64 counter = 0;
    for(u32 u = 0; u < N; u++) {
        for(u32 v = u+1; v < N; v++) {
            fxx pd2 = eudist2(D+3*u, D+3*v);
            if(pd2 < d2) {
                cb(u, v, pd2, cb_data);
                counter++;
            }
        }
    }

    if(info) {
        clock_gettime(CLOCK_REALTIME, &t1);
        info->t_scan_ms = 1000.0 * timespec_diff(&t1, &t0);
        info->n_collisions = counter;
        info->t_total_ms = info->t_scan_ms;
    }

    return EXIT_SUCCESS;
}


static int
collide3_spatial(const fxx * restrict D, const u32 N, const fxx d,
                 collide3_info * info,
                 collide3_cb cb, void * cb_data)
{
    struct timespec t0, t1, t2, t3;
    if(info) {
        clock_gettime(CLOCK_REALTIME, &t0);
    }

    if(cb == NULL) {
        cb = dummy_cb;
        cb_data = NULL;
    }

    // Make a decision about the number of buckets
    // per dimensions. Since we only care about the
    // a cubic domain, it will be the same for all dimensions

    int _nDiv = cbrt(N/7);

    _nDiv < 2 ? _nDiv = 2 : 0;
    const u32 nDiv = _nDiv;
    const u32 nDiv2 = nDiv*nDiv;
    const u32 nDiv3 = nDiv*nDiv*nDiv;

    //
    // Count sort to create the hash table, HT.
    //

    u32 * bucket_list = calloc(nDiv3+3, sizeof(u32));
    if(bucket_list == NULL) { return EXIT_FAILURE; }
    if(info) {
        info->mem_alloc = (nDiv3+3)*sizeof(u32);
    }

    // Count how many elements that will fall into each bucket.
    for(u32 kk = 0; kk<N; kk++) {
        bucket_list[hash(nDiv, D+3*kk) + 2]++;
    }

    // Integrate the list -- find the start position of each bucket
    for(u32 kk = 1; kk <= nDiv3; kk++) {
        bucket_list[kk+1] = bucket_list[kk+1] + bucket_list[kk];
    }

    // Create the actual hash table and insert copies
    // of all points tagged with their original indices
    // i.e. sort the points according to their bucket
    //
    // Points are copied to avoid cache misses in the later scanning
    // phase.
    //
    // Accumulates the the write positions
    // so that in the end writepos[kk] is startpos[kk+1]
    entry * HT = malloc(N*sizeof(entry));
    if(HT == NULL) {
        free(bucket_list);
        return EXIT_FAILURE;
    }
    if(info) {
        info->mem_alloc += N*sizeof(entry);
    }

    // Copy points and their indexes to the table
    for(u32 kk = 0; kk<N; kk++) {
        const u32 h = hash(nDiv, D+3*kk);
        const u32 ht_pos = bucket_list[h+1]; // hash table position
        bucket_list[h+1]++;
        HT[ht_pos].idx = kk;
        HT[ht_pos].X[0] = D[3*kk];
        HT[ht_pos].X[1] = D[3*kk+1];
        HT[ht_pos].X[2] = D[3*kk+2];
    }


    // To sort the points within each bucket reduces the scan time by
    // about 1% but increases the construction time by about 50%

    // The bucket list does now contain the start positions of the
    // bins and is sorted. This could be compressed quite much using an
    // Elias-Fano encoding or similar.

    if(info) {
        clock_gettime(CLOCK_REALTIME, &t1);
        info->t_create_ms = 1000.0 * timespec_diff(&t1, &t0);
        clock_gettime(CLOCK_REALTIME, &t2);
    }

    // Loop over the elements of the hash table, and use
    // the callback for each detected collision.
    //
    // The bins that we visit are those that could contain
    // a hit under the Manhattan distance, i.e., bins that
    // can't contain any hits might be included.
    //
    // Although it makes sense to discard bins using an axis-aligned
    // (aa) box vs sphere test, that costs more than scanning a few
    // extra bins unless the search radius is much larger than the
    // side length of each bin, i.e., 2/nDiv.
    u64 counter = 0;

    const fxx d2 = C3_SQUARE(d);

    for(u32 u = 0; u<N; u++) {
        // Figure out which bins might contain a hit
        // note that this requires that the hash function handles out of bounds
        const u32 ha_min = hash_coord(nDiv, HT[u].X[0] - d);
        const u32 ha_max = hash_coord(nDiv, HT[u].X[0] + d);
        const u32 hb_min = hash_coord(nDiv, HT[u].X[1] - d);
        const u32 hb_max = hash_coord(nDiv, HT[u].X[1] + d);
        const u32 hc = hash_coord(nDiv, HT[u].X[2]); // will never scan backward in the last dimension
        const u32 hc_max = hash_coord(nDiv, HT[u].X[2] + d);

        // Loop over the "worst" i.e. most strided dimension first
        for(u32 cc = hc; cc <= hc_max; cc++) {
            for(u32 bb = hb_min; bb <= hb_max; bb++) {

                // index of first element to compare with
                u32 v = bucket_list[ha_min + bb*nDiv + cc*nDiv2];

                // we only compare to later elements in the table
                // ensures no duplicate collision detection
                v <= u ? v = u+1 : 0;

                // compare to all elements before v_stop
                const u32 v_stop = bucket_list[ha_max + bb*nDiv + cc*nDiv2 + 1];

                for( ; v < v_stop; v++) {
                    const fxx pd2 = eudist2(HT[u].X, HT[v].X);
                    if(pd2 < d2) {
                        cb(HT[v].idx, HT[u].idx, pd2, cb_data);
                        counter++;
                    }
                }
            }
        }
    }

    free(bucket_list);

    free(HT);
    if(info) {
        clock_gettime(CLOCK_REALTIME, &t3);
        info->t_scan_ms = 1000.0 * timespec_diff(&t3, &t2);
        info->t_total_ms = 1000.0 * timespec_diff(&t3, &t0);
        info->n_collisions = counter;
    }
    return EXIT_SUCCESS; // == success
}

static int
collide3_spatial_lowmem(const fxx * D, const u32 N, const fxx d,
                        collide3_info * info,
                        collide3_cb cb, void * cb_data)
{
    struct timespec t0, t1, t2, t3;
    if(info) {
        clock_gettime(CLOCK_REALTIME, &t0);
    }

    if(cb == NULL) {
        cb = dummy_cb;
        cb_data = NULL;
    }

    const fxx d2 = d*d;

    int nDiv = cbrt(N/5);
    nDiv < 2 ? nDiv = 2 : 0;
    u32 n_buckets = nDiv*nDiv*nDiv;

    //
    // Count sort to create the hash table, HT.
    //

    // We will reference this array under several aliases and offsets,
    // use paper and pen to figure out!
    u32 * bucket_list = calloc(n_buckets+3, sizeof(u32));
    if(bucket_list == NULL) { return EXIT_FAILURE; }
    if(info) {
        info->mem_alloc = (n_buckets+3)*sizeof(u32);
    }
    u32 * bucket_size = bucket_list + 2;

    // Count how many elements that will fall into each bucket.
    for(u32 kk = 0; kk<N; kk++) {
        bucket_size[hash(nDiv, D+3*kk)]++;
    }

    // Integrate the list -- find the start position of each bucket
    u32 * bucket_writepos = bucket_list + 1;
    for(u32 kk = 1; kk<=n_buckets; kk++) {
        bucket_writepos[kk] = bucket_size[kk-1]+bucket_writepos[kk-1];
    }

    // Create the actual hash table and insert copies
    // of all points tagged with their original indices
    // i.e. sort the points according to their bucket
    //
    // Points are copied to avoid cache misses in the later scanning
    // phase.
    //
    // Accumulates the the write positions
    // so that in the end writepos[kk] is startpos[kk+1]
    u32 * HT = malloc(N*sizeof(u32));
    if(info) {
        info->mem_alloc += N*sizeof(u32);
    }
    if(HT == NULL) {
        free(bucket_list);
        return EXIT_FAILURE;
    }

    for(u32 kk = 0; kk<N; kk++) {
        u32 h = hash(nDiv, D+3*kk);
        u32 ht_pos = bucket_writepos[h]; // hash table position
        bucket_writepos[h]++;
        HT[ht_pos]= kk;
    }

    u32 * bucket_start = bucket_list;

    if(info) {
        clock_gettime(CLOCK_REALTIME, &t1);
        info->t_create_ms = 1000.0 * timespec_diff(&t1, &t0);
        clock_gettime(CLOCK_REALTIME, &t2);
    }

    // The loop is over the elements of the hash table,
    // i.e. not over the order that the elements were given
    // under the idea that this might reduce cache misses
    u64 counter = 0;
    for(u32 kk = 0; kk<N; kk++) {
        // Figure out which bins might contain a hit
        const fxx deps = d;
        const u32 u = HT[kk];
        const fxx * uX = D + 3*u;
        const u32 ha_min = hash_coord(nDiv, uX[0] - deps);
        const u32 ha_max = hash_coord(nDiv, uX[0] + deps);
        const u32 hb_min = hash_coord(nDiv, uX[1] - deps);
        const u32 hb_max = hash_coord(nDiv, uX[1] + deps);
        const u32 hc_min = hash_coord(nDiv, uX[2]); // will never scan backward in the last dimension
        const u32 hc_max = hash_coord(nDiv, uX[2] + deps);

        // Loop over the "worst" i.e. most strided dimension first
        for(u32 cc = hc_min; cc <= hc_max; cc++) {
            for(u32 bb = hb_min; bb <= hb_max; bb++) {

                // last bucket to visit
                const u32 hash1 =
                    ha_max +
                    bb*nDiv +
                    cc*C3_SQUARE(nDiv);

                if(bucket_start[hash1+1]-1 < kk) continue;

                // first bucket to visit
                const u32 hash0 =
                    ha_min +
                    bb*nDiv +
                    cc*C3_SQUARE(nDiv);

                // index of first element to compare with
                u32 ht_start = bucket_start[hash0];
                // we only compare to later elements in the table
                ht_start <= kk ? ht_start = kk+1 : 0;
                // index of last element to compare with
                const u32 ht_end = bucket_start[hash1+1];
                for(u32 pp = ht_start; pp < ht_end; pp++) {
                    u32 v = HT[pp];
                    const fxx * vX = D+3*v;
                    const fxx pd2 = eudist2(uX, vX);
                    if(pd2 < d2) { // squared distances
                        cb(u,v, pd2, cb_data);
                        counter++;
                    }
                }
            }
        }
    }
    free(bucket_list);
    free(HT);
    if(info) {
        clock_gettime(CLOCK_REALTIME, &t3);
        info->t_scan_ms = 1000.0 * timespec_diff(&t3, &t2);
        info->n_collisions = counter;
        info->t_total_ms = 1000.0 * timespec_diff(&t3, &t0);
    }
    return EXIT_SUCCESS; // == success
}

// Public API entry point
int
fxx(collide3)(const fxx * D, const u32 N, const fxx d,
              collide3_info * info,
              collide3_cb cb, void * cb_data)
{
    if(info && (info->backend == be_brute_force)){
        return collide3_brute_force(D, N, d, info, cb, cb_data);
    }

    if((info == NULL) || info->backend == be_spatial_lowmem) {
        if(N < 64){
            return collide3_brute_force(D, N, d, info, cb, cb_data);
        }
        return collide3_spatial_lowmem(D, N, d, info, cb, cb_data);
    }

    if((info == NULL) || info->backend == be_auto){
        if(N < 60){
            return collide3_brute_force(D, N, d, info, cb, cb_data);
        }
    }

    return collide3_spatial(D, N, d, info, cb, cb_data);
}
