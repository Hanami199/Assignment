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
inline int inject_energy( const int      periodic, // if the boundary is periodic
                        const int      Nsources, // number of sources for the rank
			            const vec2_t  *Sources, // positions of those sources
			            const double   energy, // energy per source
                        plane_t *plane, // plane on which to inject the energu
                        const vec2_t   N) // process grid shape
{
    // plane->size[_x_]=nx (interior length), plane->size[_y_]=ny (interior height)
    const uint register sizex = plane->size[_x_]+2; // interior + halos
    double * restrict data = plane->data;
    
    // used to ready from flattened data array
    #define IDX( i, j ) ( (j)*sizex + (i) )

    // for every source of this rank
    const int nx = plane->size[_x_];
    const int ny = plane->size[_y_];
    const int wrap_x = periodic && (N[_x_] == 1);
    const int wrap_y = periodic && (N[_y_] == 1);

    // loop over sources
    for (int s = 0; s < Nsources; ++s) {
        const int x = Sources[s][_x_];
        const int y = Sources[s][_y_];

        // always add to the interior point
        data[ IDX(x, y) ] += energy;

        // only do halo propagation when needed
        if (wrap_x) {
            // note: keep as two independent if's in case nx==1 (both can be true)
            if (x == 1)  data[ IDX(nx + 1, y) ] += energy; // east interior mirrors to west halo
            if (x == nx) data[ IDX(0,y) ] += energy; // west interior mirrors to east halo
        }
        if (wrap_y) {
            if (y == 1) data[ IDX(x, ny + 1) ] += energy; // north interior -> south halo
            if (y == ny)  data[ IDX(x, 0)] += energy; // south interior -> north halo
        }
    }
    #undef IDX
    return 0;
}


inline int update_plane ( const int periodic, // toggle for periodic-nonperiodic boundaries 
                        const vec2_t N,         // the grid of MPI tasks
                        const plane_t *oldplane, // plane to take values for stencil computation from
                        plane_t *newplane // plane to write results of stencil computation
                        ) 

{
    uint register fxsize = oldplane->size[_x_]+2; // interior + halo
    uint register fysize = oldplane->size[_y_]+2;
    
    uint register xsize = oldplane->size[_x_]; // just interior
    uint register ysize = oldplane->size[_y_];
    
    #define IDX( i, j ) ( (j)*fxsize + (i) ) // preprocessor macro to access flattened array

    double * restrict old = oldplane->data; // we will read from this
    double * restrict new = newplane->data; // and write on this
    
    #pragma omp parallel for schedule(static)
    for (uint j = 1; j <= ysize; j++) {
        for ( uint i = 1; i <= xsize; i++)
        {
            new[ IDX(i,j) ] = old[ IDX(i,j) ] * 0.5 + ( old[IDX(i-1, j)] + old[IDX(i+1, j)] +
                                old[IDX(i, j-1)] + old[IDX(i, j+1)] ) * 0.125;
            
        }
    }

    if ( periodic ) {
        // we need to propagate only when we are alone row-wise or column-wise in the
        // grid of processors
            
        if ( N[_x_] == 1 ) {
            // propagate the boundaries as needed
            // check the serial version
            for (uint j = 1; j <= ysize; j++) {
                new[ IDX(0,j) ] = new[ IDX(xsize,j) ]; // we propagate east edge on west halo
                new[ IDX(xsize+1,j) ] = new[ IDX(1,j) ]; // we propagate west edge on east halo
            }
        }

        if ( N[_y_] == 1 ) {
                // propagate the boundaries as needed
                // check the serial version
            for (uint i = 1; i <= xsize; i++) {
                new[ IDX(i, 0) ] = new[ IDX(i, ysize) ]; // propagate south edge on the north halo
                new[ IDX(i, ysize+1) ] = new[ IDX(i,1) ]; // propagate north edge on the south halo
            }
        }

        
    }

    #undef IDX
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



