/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 * See COPYRIGHT in top-level directory.
 */

#define _XOPEN_SOURCE 700
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <math.h>

#include <omp.h>
#include <mpi.h>
#define _DEFAULT_SOURCE

#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3

#define SEND 0
#define RECV 1

#define OLD 0
#define NEW 1

#define _x_ 0
#define _y_ 1


typedef unsigned int uint;

typedef uint    vec2_t[2]; // array used for 2D tuples
typedef double * restrict buffers_t[4]; // array of pointers (one per direction)

typedef struct {
    double   * restrict data; // a flattened 2D arrays with 1 cell halo on each side
    vec2_t     size; // interior size [xsize, ysize]
} plane_t; // grid


// add heat at specific coordinates on a given plane
extern int inject_energy ( const int, // periodic or not
                        const int, // Nsources
			            const vec2_t *, // sources
			            const double, // energy
                        plane_t *, // plane
                        const vec2_t); // N

// apply the 5-point-stencil once
extern int update_plane ( const int      , // periodic or not
                          const vec2_t   , // N=[N[_x_], N[_y_]]
                          const plane_t *, // oldplane
                                plane_t *  //newplane
                        ); 

// sums cells to get total heat
extern int get_total_energy( plane_t *, // plane
                             double  * ); // energy

// extern int dump ( const double *data, const uint size[2], const char *filename);

// function used to initialize everything (planes buffers, comm buffers etc.)
int initialize ( MPI_Comm *,
                int       ,
		        int       ,
		        int       ,
		        char    **,
                vec2_t   *,
                vec2_t   *,                 
		        int      *,
                int      *,
		        int      *,
		        int      *,
		        int      *,
		        int      *,
                vec2_t  **,
                double   *,
                plane_t  *,
                buffers_t * );

// free memory associated with the plane
int memory_release ( buffers_t *,
                      plane_t   *);

// computes and prints energy at a step
int output_energy_stat ( int      ,
                         plane_t *,
                         double   ,
                         int      ,
                         MPI_Comm *);

extern int fill_buffers(buffers_t *, 
                        const plane_t *,
                        int,
                        const vec2_t N);

extern int post_MPI_reqs(MPI_Request *,
                  buffers_t *,
                  const plane_t *,
                  int *,
                  MPI_Comm );

extern int copy_halos(buffers_t *buffers, 
              plane_t *plane,
              int* neigh,
              int,
            const vec2_t N);

// function used to inject energy (from heat sources) onto a plane
inline int inject_energy(const int      periodic,
                         const int      Nsources,
                         const vec2_t * restrict Sources,   // restrict helps alias analysis
                         const double   energy,
                         plane_t *plane,
                         const vec2_t   N)                  // process grid shape
{
    const unsigned int nx = plane->size[_x_];   // interior width
    const unsigned int ny = plane->size[_y_];   // interior height
    const size_t sizex   = (size_t)nx + 2;      // interior + 2 halos
    double * restrict data = plane->data;

    // If periodic and we have only one rank along an axis, we must wrap locally.
    const int wrapX = periodic && (N[_x_] == 1);
    const int wrapY = periodic && (N[_y_] == 1);

    for (int s = 0; s < Nsources; ++s) {
        const unsigned int x = Sources[s][_x_]; // expected 1..nx
        const unsigned int y = Sources[s][_y_]; // expected 1..ny
        const size_t row = (size_t)y * sizex;

        // center point
        data[row + x] += energy;

        if (wrapX) {
            if (x == 1)              data[row + (nx + 1)] += energy;  // west interior → east halo
            else if (x == nx)        data[row + 0]         += energy;  // east interior → west halo
        }

        if (wrapY) {
            if (y == 1)              data[(size_t)(ny + 1) * sizex + x] += energy; // north interior → south halo
            else if (y == ny)        data[0 * sizex + x]                += energy; // south interior → north halo
        }

        // Corner halo when source is on both edges and both axes wrap locally.
        if (wrapX && wrapY && (x == 1 || x == nx) && (y == 1 || y == ny)) {
            const size_t halo_x = (x == 1) ? (size_t)(nx + 1) : 0;
            const size_t halo_y = (y == 1) ? (size_t)(ny + 1) : 0;
            data[halo_y * sizex + halo_x] += energy;
        }
    }
    return 0;
}


inline int update_plane(const int periodic,         // periodic boundaries?
                        const vec2_t N,             // MPI process grid shape
                        const plane_t *oldplane,    // read from
                        plane_t *newplane)          // write to
{
    const unsigned int nx   = oldplane->size[_x_];      // interior
    const unsigned int ny   = oldplane->size[_y_];
    const size_t       sx   = (size_t)nx + 2;           // stride incl. halos

    // If you can guarantee 64B alignment at allocation time,
    // uncomment the assume_aligned lines below.
    double * restrict oldp = /*(double*)__builtin_assume_aligned(*/ oldplane->data /*,64)*/;
    double * restrict newp = /*(double*)__builtin_assume_aligned(*/ newplane->data /*,64)*/;

    const double c0 = 0.5;    // center weight
    const double c1 = 0.125;  // neighbor weight (sum of 4)

    // Parallelize over rows; inner loop is contiguous -> good SIMD
    // If ny is tiny and nx is huge, keep as-is. If both are large-ish and balanced,
    // you can try: #pragma omp parallel for collapse(2) schedule(static)
    #pragma omp parallel for schedule(static)
    for (unsigned int j = 1; j <= ny; ++j) {
        // row pointers (use halos at j-1 and j+1; they exist)
        const double * restrict up    = &oldp[(size_t)(j-1) * sx];
        const double * restrict mid   = &oldp[(size_t) j    * sx];
        const double * restrict down  = &oldp[(size_t)(j+1) * sx];
              double * restrict out   = &newp[(size_t) j    * sx];

        // Hint vectorization; add 'aligned(mid,out,up,down:64)' if you can guarantee it
        #pragma omp simd /* aligned(mid,out,up,down:64) */ linear(i:1) safelen(16)
        for (unsigned int i = 1; i <= nx; ++i) {
            // new[j,i] = 0.5*old + 0.125*(left+right+up+down)
            out[i] = c0 * mid[i]
                   + c1 * (mid[i-1] + mid[i+1] + up[i] + down[i]);
        }

        // optional: light prefetch for next row (tune distance if helpful)
        // if (j+1 <= ny) {
        //     __builtin_prefetch(&oldp[(size_t)(j+2)*sx + 16], 0, 1);
        //     __builtin_prefetch(&oldp[(size_t)(j+1)*sx + 16], 0, 1);
        //     __builtin_prefetch(&newp[(size_t)(j+1)*sx + 16], 1, 1);
        // }
    }

    if (periodic) {
        // If there is only one MPI rank along X, wrap left/right locally
        if (N[_x_] == 1) {
            for (unsigned int j = 1; j <= ny; ++j) {
                const size_t r = (size_t)j * sx;
                newp[r + 0     ] = newp[r + nx];  // west halo  <- east edge
                newp[r + (nx+1)] = newp[r + 1 ];  // east halo  <- west edge
            }
        }
        // If there is only one MPI rank along Y, wrap top/bottom locally
        if (N[_y_] == 1) {
            // rows are contiguous → memcpy is faster than a loop
            memcpy(&newp[0 * sx      + 1], &newp[(size_t)ny    * sx + 1], nx * sizeof(double));  // north halo  <- south edge
            memcpy(&newp[(size_t)(ny+1)*sx + 1], &newp[1 * sx  + 1],     nx * sizeof(double));  // south halo  <- north edge
        }
        // Corners not needed for 5-point stencil (no diagonal reads).
    }

    return 0;
}


// function to compute the total energy on the plane 
inline int get_total_energy( plane_t *plane,
                             double  *energy )
/*
 * NOTE: this routine a good candiadate for openmp
 *       parallelization
 */
{

    const int register xsize = plane->size[_x_];
    const int register ysize = plane->size[_y_];
    const int register fsize = xsize+2;

    double * restrict data = plane->data;
    
    #define IDX( i, j ) ( (j)*fsize + (i) )

    #if defined(LONG_ACCURACY)    
        long double totenergy = 0;
    #else
        double totenergy = 0;    
    #endif

    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    #pragma omp parallel for reduction(+:totenergy) schedule(static)
    for ( int j = 1; j <= ysize; j++ ) {
        for ( int i = 1; i <= xsize; i++ ) {
            totenergy += data[ IDX(i, j) ];
        }
    }

    #undef IDX

    *energy = (double)totenergy;
    return 0;
}



