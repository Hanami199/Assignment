/*

/*
 *
 *  mysizex   :   local x-extendion of your patch
 *  mysizey   :   local y-extension of your patch
 *
 */

#define _XOPEN_SOURCE 700
#include "stencil_template_parallel.h"


// ------------------------------------------------------------------
// ------------------------------------------------------------------

// EVERY rank will execute main():

int main(int argc, char **argv) {

  /* -----------1. INITIALISATION ---------------- */
  MPI_Comm myCOMM_WORLD; // communicator
  int  Rank, Ntasks; // rank=this process ID, Ntasks=total # of MPI ranks
  uint neighbours[4]; // neighbour ranks (NORTH-SOUTH-EAST-WEST)

  int  Niterations; // maximum number of iterations
  int  periodic; // 0-1 flag to define whether the boundary is periodic
  vec2_t S, N; // S=global grid, N=process grid
  
  int      Nsources; // total number of heat sources (global)
  int      Nsources_local; // number of heat sources in THIS rank
  vec2_t  *Sources_local; // local sources coordinates
  double   energy_per_source; // energy injected per source

  
  plane_t   planes[2]; // array for plane[NEW] and plane[OLD]
  buffers_t buffers[2]; // communication buffers
  
  int output_energy_stat_perstep; // 0-1: print local energy at each step?

  /* ------------- initialize MPI envrionment ----------------*/
  int level_obtained;
  
  // NOTE: change MPI_FUNNELED if appropriate

  // with FUNNELED we ALLOW multithreading (through OpenMP) but
  // only the master thread makes the MPI calls, so
  // they do compute in parallel, then funnel all MPI to master thread 
  MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );

  // if FUNNELED is not possible, we return error
  if ( level_obtained < MPI_THREAD_FUNNELED ) {
    printf("MPI_thread level obtained is %d instead of %d\n",
    level_obtained, MPI_THREAD_FUNNELED );
    MPI_Finalize();
    exit(1); 
  }
  
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank); // get my ID (of the rank)
  MPI_Comm_size(MPI_COMM_WORLD, &Ntasks); // get total number of ranks
  MPI_Comm_dup(MPI_COMM_WORLD, &myCOMM_WORLD); // make a private communicator
  // we will have the same members (the other ranks, since we called MPI_COMM_WORLD) 
  // but a different context,
  // so that we don't interfere with others (same members, new context)
  
  
  /* ------------------ argument checking and setting -----------------------*/

  // we call the setup routing:
  int ret = initialize( &myCOMM_WORLD, Rank, Ntasks, argc, argv, &S, &N, &periodic, 
      &output_energy_stat_perstep,
			neighbours, &Niterations,
			&Nsources, &Nsources_local, &Sources_local, &energy_per_source,
			&planes[0], &buffers[0] );
      // sending &planes[0] or &planes is the same since they are arrays 
  printf("Rank %d signals initialization completed. Beginning computation part.\n", Rank);
  // if the init failed (ret != 0), then print a message, shut down MPI and exit
  if ( ret ) {
    printf("task %d is opting out with termination code %d\n",
      Rank, ret );
    MPI_Finalize();
    return 0;
  }
  


  /* ----------- 2. COMMUNICATION & COMPUTATION PART ---------------- */

  int current = OLD;  // set the buffer at "read" -> read from planes[OLD], write on NEW
  double t1 = MPI_Wtime();   /* take wall-clock time */
  double comm_sum = 0.0; // this will use to compute mean comm waiting time
  double copy_sum      = 0.0; // this for copy wait time
  int dump_every = 5;

  // to be reapeated for Niterations times:
  for (int iter = 0; iter < Niterations; ++iter) { 
    MPI_Request reqs[8]; // vector of handles returned by MPI calls
    // we store it so that we can complete them with MPI_Wait
    // it's 8 because we will have 4 receives and 4 sends
    // 0=N-se, 1=S-se, 2=E-se, 3=W-se, 4=N-re, 5=S-re, 6=E-re, 7=W-re

    int flag = 0; // 8 flags for SEND/RECV

    // we initialize all reqs to MPI_REQUEST_NULL
    for (int k = 0; k < 8; ++k) reqs[k] = MPI_REQUEST_NULL;
    
    /* we inject new energy from sources */
    double t0 = MPI_Wtime();
    inject_energy( periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N );


    /* ------------ COMMUNICATIONS & COMPUTATIONS -------------------------- */

    // [A] fill the buffers, and/or make the buffers' pointers pointing to the correct position

    fill_buffers(buffers, &planes[current], periodic);

    // We initialize Isend, Irecv (since before #pragma omp parallel we only have the master thread)
    post_MPI_reqs(reqs, buffers, &planes[current], neighbours, myCOMM_WORLD);
    
    #pragma omp parallel /*shared(flag)*/
    {
      if (omp_get_thread_num() == 0) {
        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);

        // time to copy the halos
        double t1 = MPI_Wtime();
        copy_halos(buffers, &planes[current], neighbours, periodic);
        double t2 = MPI_Wtime();

        // only master writes these 
        comm_sum += t1 - t0;
        copy_sum += t2 - t1;

      } 

      #pragma omp barrier
      update_plane(periodic, N, &planes[current], &planes[!current]);
    }

    /* output if needed */
    if ( output_energy_stat_perstep )
      output_energy_stat( iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD);

    char filename[100];
    sprintf( filename, "data_parallel/%d_plane_%05d.bin", Rank, iter );
    int dump_status = dump(planes[!current].data, planes[!current].size, filename);
    if (dump_status != 0)
    {
      fprintf(stderr, "Error in dump_status. Exit with %d\n", dump_status);
    }

    /* swap plane indexes for the new iteration */
    current = !current; 
    }
  
  // after all the comms-computes
  t1 = MPI_Wtime() - t1;
  printf("Total execution time for rank %d: %f\n", Rank, t1);

  // compute mean waiting comm time and copy time
  double comm_mean = comm_sum/(double)Niterations;
  double copy_mean = copy_sum/(double)Niterations;

  // make a global mean 
  double comm_sum_all = 0.0, copy_sum_all = 0.0;
  MPI_Reduce(&comm_sum, &comm_sum_all, 1, MPI_DOUBLE, MPI_SUM, 0, myCOMM_WORLD);
  MPI_Reduce(&copy_sum, &copy_sum_all, 1, MPI_DOUBLE, MPI_SUM, 0, myCOMM_WORLD);

  if (Rank == 0) {
    int P; MPI_Comm_size(myCOMM_WORLD, &P);
    double comm_mean_all = comm_sum_all/(Niterations * (double)P);
    double copy_mean_all = copy_sum_all/(Niterations * (double)P);

    printf("Average comm waiting time and copy time (averaged across %d ranks): wait=%.6fs, copy=%.6fs\n",
          P, comm_mean_all, copy_mean_all);
  }
  output_energy_stat ( -1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );


  memory_release( buffers, planes );
  
  MPI_Finalize();
  return 0;
}


/* ==========================================================================
   =                                                                        =
   =   routines called within the integration loop                          =
   ========================================================================== */

   
// Space for subroutines ... if only i had any....

// function to fill the buffers (the only one that need filling are the EAST and WEST ones)
int fill_buffers(buffers_t *buffers, 
                      const plane_t *plane,
                      int periodic) {

  const unsigned int nx = plane->size[_x_];
  const unsigned int ny = plane->size[_y_];
  const size_t fx = (size_t)nx + 2;

  #define IDX(i,j) ((j)*fx + (i))

  // rows are contiguos so we just point in the correct position:
  buffers[SEND][NORTH] = &plane->data[IDX(1,1)]; // north starts from position (1,1) (since we also have the halos in data)
  buffers[SEND][SOUTH] = &plane->data[IDX(1,ny)]; // north starts from position (1,ny)

  // we will receive in the plane directly
  buffers[RECV][NORTH] = &plane->data[IDX(1,0)]; // north starts from position (1,1) (since we also have the halos in data)
  buffers[RECV][SOUTH] = &plane->data[IDX(1,ny+1)]; 

  if (periodic) {
    buffers[RECV][NORTH] = &plane->data[IDX(1,ny+1)];  // north starts from position (1,1) (since we also have the halos in data)
    buffers[RECV][SOUTH] = &plane->data[IDX(1,0)]; 
  }

  // west column (all the values on the left side, so where x=1 and for y=1,...,ny)
  if (buffers[SEND][WEST]) {
    for (uint j = 0; j < ny; j++)
      buffers[SEND][WEST][j] = plane->data[IDX(1, j+1)];
  }

  // east column (all the values on the right side, so where x=nx and for y=1,...,ny)
  if (buffers[SEND][EAST]) {
    for (uint j = 0; j < ny; j++)
      buffers[SEND][EAST][j] = plane->data[IDX(nx, j+1)];
  }

  #undef IDX
  return 0;
}

// function to post all the non-blocking calls needed for computation
int post_MPI_reqs(MPI_Request *reqs, 
                  buffers_t *buffers, 
                  const plane_t *plane, 
                  int *my_neigh, 
                  MPI_Comm comm) {

  const unsigned int nx = plane->size[_x_];
  const unsigned int ny = plane->size[_y_];

  // first 4 are 0=N-se, 1=S-se, 2=E-se, 3=W-se, then 4=N-re, 5=S-re, 6=E-re, 7=W-re
  // also chose a different tag for direction (0=N, 1=S, 2=E, 3=W)
  for (uint d=0; d<2; d++) {
    if (my_neigh[d] != MPI_PROC_NULL) {
      MPI_Irecv(buffers[RECV][d], nx, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d+4]);
      MPI_Isend(buffers[SEND][d], nx, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d]); 
    }
  }

  // P.S East and West have different dimension (ny, not nx)
  for (uint d=2; d<4; d++) {
    if (my_neigh[d] != MPI_PROC_NULL) {
      MPI_Irecv(buffers[RECV][d], ny, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d+4]); 
      MPI_Isend(buffers[SEND][d], ny, MPI_DOUBLE, my_neigh[d], 0,  comm, &reqs[d]);  
    }
  }

  return 0;
}

// function to copy halos from buffers to planes when needed
int copy_halos(buffers_t *buffers, 
              plane_t *plane,
              int* neigh,
              int periodic) {

  const unsigned int nx = plane->size[_x_];
  const unsigned int ny = plane->size[_y_];
  const size_t fx = (size_t)nx + 2;

  #define IDX(i,j) ((j)*fx + (i))

  if (neigh[WEST] != MPI_PROC_NULL) {
  for (uint j=0; j<ny; j++){
    plane->data[IDX(0, j+1)] = buffers[RECV][WEST][j];
    }
  }

  if (neigh[EAST] != MPI_PROC_NULL) {
    for (uint j=0; j<ny; j++){
    plane->data[IDX(nx+1, j+1)] = buffers[RECV][EAST][j];
    }
  }

  // in the periodic case we are working in a sort of torus
  // if (periodic) {
  //   if (neigh[WEST] != MPI_PROC_NULL) {
  //     for (uint j=0; j<ny; j++){
  //     plane->data[IDX(nx+1, j+1)] = buffers[RECV][WEST][j];
  //     }
  //   }

  //   if (neigh[EAST] != MPI_PROC_NULL) {
  //     for (uint j=0; j<ny; j++){
  //     plane->data[IDX(0, j+1)] = buffers[RECV][EAST][j];
  //     }
  //   }
  // }

  #undef IDX
  return 0;
}


/* ==========================================================================
   =                                                                        =
   =   initialization                                                       =
   ========================================================================== */


uint simple_factorization( uint, int *, uint ** );

int initialize_sources( int       ,
			int       ,
			MPI_Comm  *,
			uint      [2],
			int       ,
			int      *,
			vec2_t  ** );


int memory_allocate ( const int       *,
		      const vec2_t     ,
		            buffers_t *,
		            plane_t   * );
		      

int initialize ( MPI_Comm *Comm,                // communicator
                  int      Me,                  // the rank of the calling process
                  int      Ntasks,              // the total number of MPI ranks
                  int      argc,                // the argc from command line
                  char   **argv,                // the argv from command line
                  vec2_t  *S,                   // the size of the plane
                  vec2_t  *N,                   // two-uint array defining the MPI tasks' grid
                  int     *periodic,            // periodic-boundary tag
                  int     *output_energy_stat,
                  int     *neighbours,          // four-int array that gives back the neighbours of the calling task
                  int     *Niterations,         // how many iterations
                  int     *Nsources,            // how many heat sources
                  int     *Nsources_local,
                  vec2_t **Sources_local,
                  double  *energy_per_source,   // how much heat per source
                  plane_t *planes,
                  buffers_t *buffers
                  ) {

  int halt = 0;
  int ret; // this will be returned at the end, to show whether the init was successful
  int verbose = 0;
  
// ··········1. SET DEFAULT VALUES·······················

  // Default grid will be 10000x10000
  (*S)[_x_]         = 10000;
  (*S)[_y_]         = 10000;

  // We set (by default):
  *periodic         = 0; // non-periodic boundary
  *Nsources         = 4; // number of global heat sources
  *Nsources_local   = 0; // number of local (rank) heat sources (will be set with initialize_sources)
  *Sources_local    = NULL; // positions of local heat sources
  *Niterations      = 1000; // number of iterations
  *energy_per_source = 1.0; // amount of energy per source (injected at every step)


  if ( planes == NULL ) {
    // manage the situation
    fprintf(stderr, "initialisation error: planes pointer is NULL\n");
    return 1;
  }

  // interior size "not set yet" so set to 0
  planes[OLD].size[0] = planes[OLD].size[1] = 0; // both 0 and 1??
  planes[NEW].size[0] = planes[NEW].size[1] = 0;
  
  // we set no neighbour by default until we add one
  for ( int i = 0; i < 4; i++ ) {
    neighbours[i] = MPI_PROC_NULL;
  }

  // all comm buffers are NULL, because they are not allocated yet
  for ( int b = 0; b < 2; b++ ) { // SEND and RECV
    for ( int d = 0; d < 4; d++ ) { // N/S/E/W
      buffers[b][d] = NULL; 
    }
  }
  
  // ··································································
  // process the commadn line
  // 
  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":hx:y:e:E:n:o:p:v:")) != -1) {
      switch( opt )
        {
        case 'x': (*S)[_x_] = (uint)atoi(optarg);
          break;

        case 'y': (*S)[_y_] = (uint)atoi(optarg);
          break;

        case 'e': *Nsources = atoi(optarg);
          break;

        case 'E': *energy_per_source = atof(optarg);
          break;

        case 'n': *Niterations = atoi(optarg);
          break;

        case 'o': *output_energy_stat = (atoi(optarg) > 0); // 0 or 1
          break;

        case 'p': *periodic = (atoi(optarg) > 0); // 0 or 1
          break;

        case 'v': verbose = atoi(optarg);
          break;

        case 'h': {
          if ( Me == 0 )
            printf( "\nvalid options are ( values btw [] are the default values ):\n"
              "-x    x size of the plate [10000]\n"
              "-y    y size of the plate [10000]\n"
              "-e    how many energy sources on the plate [4]\n"
              "-E    how many energy sources on the plate [1.0]\n"
              "-n    how many iterations [1000]\n"
              "-p    whether periodic boundaries applies  [0 = false]\n\n");
          halt = 1; }
          break;
          
        case ':': printf( "option -%c requires an argument\n", optopt);
          break;
          
        case '?': printf(" -------- help unavailable ----------\n");
          break;
        }
      }

    if ( opt == -1 )
      break;
  }

  if ( halt )
    return 1;

  
  // ··················· parameter checks ···············································
  /* here we should check for all the parms being meaningful */

  if ( (*S)[_x_] <= 0 || (*S)[_y_] <=0) {
    printf( "invalid grid input, please select a positive non-zero integer");
    return 1; 
  }
  
  if (*Niterations <= 0) {
    printf( "invalid number of iterations, please select a positive non-zero integer");
    return 1;
  }

  if (*Nsources < 0) {
    printf( "invalid number of heat sources, please select a positive integer");
    return 1;
  }

  if (verbose < 0) {
    printf( "invalid verbose option, please select a positive integer");
    return 1;
  }

  // output energy and periodic are already coerced to boolean

  
// ·························· DOMAIN DECOMPOSITION ········································
/*
  * find a suitable domain decomposition
  * very simple algorithm, you may want to
  * substitute it with a better one
  *
  * the plane Sx x Sy will be solved with a grid
  * of Nx x Ny MPI tasks
  */

vec2_t Grid; // it's the obtained process grid
// (*S)[_x_], (*S)[_y_] is the global grid size
// Ntasks is the total MPI ranks P (what is passed with mpirun -np P)

// 1. PRIME FACTOR HEURISTIC

// formfactor is used to decide if we can enough ranks to tile in 2D
double formfactor = ((*S)[_x_] >= (*S)[_y_] 
                    ? (double)(*S)[_x_]/(*S)[_y_] 
                    : (double)(*S)[_y_]/(*S)[_x_] );

// if Ntasks is less than the lower floor of formfactor+1, we use a 1D split                    
int dimensions = 2 - (Ntasks <= ((int)formfactor+1) );

// they will be put row or column wise
if ( dimensions == 1 ) {
    if ( (*S)[_x_] >= (*S)[_y_] ) {
      Grid[_x_] = Ntasks, Grid[_y_] = 1;} // row
    else {
      Grid[_x_] = 1, Grid[_y_] = Ntasks;} // column
}

// else we use a simple prime factorization (defined below)
else {
  int   Nf;
  uint *factors;
  uint  first = 1;
  ret = simple_factorization( Ntasks, &Nf, &factors );
  
  for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ ) 
    first *= factors[i];

    if ( (*S)[_x_] > (*S)[_y_] ) 
      Grid[_x_] = Ntasks/first, Grid[_y_] = first;
    else 
      Grid[_x_] = first, Grid[_y_] = Ntasks/first;
  }

  (*N)[_x_] = Grid[_x_];
  (*N)[_y_] = Grid[_y_];
  

  // ··································································
  // Coordinates of a process in the process grid
  //
  int X = Me % Grid[_x_];
  int Y = Me / Grid[_x_];

  // ··································································
  // Find the neighbours of a process (some could be NULL if we are in the non-periodic case)

  if ( Grid[_x_] > 1 )
    {  
      if ( *periodic ) {       
        neighbours[EAST]  = Y*Grid[_x_] + (Me + 1) % Grid[_x_];
        neighbours[WEST]  = (X%Grid[_x_] > 0 ? Me-1 : (Y+1)*Grid[_x_]-1); }
      
      else {
        neighbours[EAST]  = ( X < Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
        neighbours[WEST]  = ( X > 0 ? (Me-1)%Ntasks : MPI_PROC_NULL ); }  
    }

  if ( Grid[_y_] > 1 )
    {
      if ( *periodic ) {      
        neighbours[NORTH] = (Ntasks + Me - Grid[_x_]) % Ntasks;
        neighbours[SOUTH] = (Ntasks + Me + Grid[_x_]) % Ntasks; }

      else {    
        neighbours[NORTH] = ( Y > 0 ? Me - Grid[_x_]: MPI_PROC_NULL );
        neighbours[SOUTH] = ( Y < Grid[_y_]-1 ? Me + Grid[_x_] : MPI_PROC_NULL ); }
    }



  // ··············  PATCH SIZE ··································
  // the size of my patch

  /*
   * every MPI task determines the size sx x sy of its own domain
   * REMIND: the computational domain will be embedded into a frame
   *         that is (sx+2) x (sy+2)
   *         the outern frame will be used for halo communication or
   */
  
  vec2_t mysize; // interior size

  // S is the global problem size
  uint s = (*S)[_x_] / Grid[_x_]; // base chuck each column gets
  uint r = (*S)[_x_] % Grid[_x_]; // leftover columns
  mysize[_x_] = s + (X < r); // the first r columns get one extra

  // ex. Gx=10, Px=3 -> s=3, r=1: each columns will have 3 chunks apart from the first, which will have 4

  // same for columns
  s = (*S)[_y_] / Grid[_y_];
  r = (*S)[_y_] % Grid[_y_];
  mysize[_y_] = s + (Y < r);

  // mysize={sx, sy} is the rank's interior size with no halos

  // both read plane (OLD) and write plane (NEW) have the same interior size
  planes[OLD].size[0] = mysize[0];
  planes[OLD].size[1] = mysize[1];
  planes[NEW].size[0] = mysize[0];
  planes[NEW].size[1] = mysize[1];
  // the allocated arrays will be (sx+2) x (sy+2) (to include the halo cells) (to be done in memory_allocate)

  if ( verbose > 0 ) {// if option was added from commandline
    if ( Me == 0 ) {
      printf("Tasks are decomposed in a grid %d x %d\n\n", Grid[_x_], Grid[_y_] ); // we just print the chosen greed
      fflush(stdout);
    }

    MPI_Barrier(*Comm); // so that one rank prints a time
    
    // each rank shows:
    // - its grid coordinates (X,Y) inside px x py
    // - the neighbor rank IDs in the N/E/S/W directions
    for ( int t = 0; t < Ntasks; t++ ) {
      if ( t == Me ) {
        printf("Task %4d: "
        "\tgrid coordinates : %3d, %3d\n"
        "\tneighbours: N %4d    E %4d    S %4d    W %4d\n",
        Me, X, Y,
        neighbours[NORTH], neighbours[EAST],
        neighbours[SOUTH], neighbours[WEST] );
        printf("My (rank %d) patch size is %d x %d\n", Me, mysize[0], mysize[1]);
        fflush(stdout);
      }

    

    MPI_Barrier(*Comm); 
    }
  }
  
  // ··································································
  // allocate the needed memory
  //
  ret = memory_allocate(neighbours, // who the N/S/E/W neighbours are
                        *N, // process grid shpae
                        buffers, // communication buffers
                        planes); //planes (OLD and NEW)
  

  // ··································································
  // allocate the heat sources
  //
  ret = initialize_sources( Me, Ntasks, Comm, mysize, 
                            *Nsources, Nsources_local, Sources_local );
  
  return 0;  
}



// we want to build a process grid px x py whose product is A and whose
// aspect ratio isn't skinnier than the domain -> to be used in initialize()
uint simple_factorization( uint A, int *Nfactors, uint **factors ) {
/*
 * rought factorization for A;
 * assumes that A is small, of the order of <~ 10^5 max,
 * since it represents the number of tasks
 */

  int N = 0; // counter for prime factors
  int f = 2; // starting from the smallest prime
  uint _A_ = A; //working copy

  // we first count how many prime factors there are
  while ( f < A ) { // we try every f from 2 to A-1
      while( _A_ % f == 0 ) { // while f divides what is left
	      N++; // count this factor
	      _A_ /= f; 
      } // divide it and keep goind
      
      // update the next possible divisor:  
      f++; 
    }

  *Nfactors = N;
  uint *_factors_ = (uint*)malloc( N * sizeof(uint) );

  // reset for second pass
  N   = 0;
  f   = 2;
  _A_ = A;

  // we fill the factors array in ascending order
  while ( f < A ) {
      while( _A_ % f == 0 ) {
	      _factors_[N++] = f;
	      _A_ /= f; }
      f++; }

  *factors = _factors_;
  return 0;
}



// function to iniatialize heat resources depending on rank
int initialize_sources( int Me,
			int       Ntasks,
			MPI_Comm *Comm,
			vec2_t    mysize,
			int       Nsources,
			int      *Nsources_local,
			vec2_t  **Sources ) {

  srand48(time(NULL) ^ Me); // get a different seed per rank

  // allocate the buffer
  int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
  
  // rank 0 randomly assigns each source to a rank
  if ( Me == 0 ) {
      for ( int i = 0; i < Nsources; i++ ) {
	      tasks_with_sources[i] = (int)lrand48() % Ntasks; }
  }
  // tasks_with_sources[s] tells which rank owns source s
  
  //we broadcast the assignment to all ranks
  MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm );

  // we count how many sources this rank owns
  int nlocal = 0;
  for ( int i = 0; i < Nsources; i++ ) {
    nlocal += (tasks_with_sources[i] == Me); 
  }
  *Nsources_local = nlocal;
  
  // if this rank owns any sources, allocate and place them
  if ( nlocal > 0 ) {
      vec2_t * restrict helper = (vec2_t*)malloc( nlocal * sizeof(vec2_t) );      
      for ( int s = 0; s < nlocal; s++ ) {
        // positions are uniform random inside the rank's interior
        helper[s][_x_] = 1 + lrand48() % mysize[_x_];
        helper[s][_y_] = 1 + lrand48() % mysize[_y_]; }

      *Sources = helper;
    }
  
  free( tasks_with_sources );

  return 0;
}



int memory_allocate(const int *neighbours,
                    const vec2_t N,
                    buffers_t *buffers_ptr ,
                    plane_t *planes_ptr) {
    /*
      here you allocate the memory buffers that you need to
      (i)  hold the results of your computation
      (ii) communicate with your neighbours

      The memory layout that I propose to you is as follows:

      (i) --- calculations
      you need 2 memory regions: the "OLD" one that contains the
      results for the step (i-1)th, and the "NEW" one that will contain
      the updated results from the step ith.

      Then, the "NEW" will be treated as "OLD" and viceversa.

      These two memory regions are indexed by *plate_ptr:

      planew_ptr[0] ==> the "OLD" region
      plamew_ptr[1] ==> the "NEW" region

      (ii) --- communications

      you may need two buffers (one for sending and one for receiving)
      for each one of your neighnours, that are at most 4:
      north, south, east amd west.      

      To them you need to communicate at most mysizex or mysizey
      daouble data.

      These buffers are indexed by the buffer_ptr pointer so
      that

      (*buffers_ptr)[SEND][ {NORTH,...,WEST} ] = .. some memory regions
      (*buffers_ptr)[RECV][ {NORTH,...,WEST} ] = .. some memory regions
      
      --->> Of course you can change this layout as you prefer
      
     */

  if (planes_ptr == NULL ) {
    // an invalid pointer has been passed
    // manage the situation
    fprintf(stderr, "invalid planes");
    return 1;}

  if (buffers_ptr == NULL ) {
    // an invalid pointer has been passed
    // manage the situation
    fprintf(stderr, "invalid buffers");
    return 1;}
    
  // ··················································
  // allocate memory for data
  // we allocate the space needed for the plane plus a contour frame
  // that will contains data form neighbouring MPI tasks
  unsigned int frame_size = (planes_ptr[OLD].size[_x_]+2) * (planes_ptr[OLD].size[_y_]+2); //(nx+2) x (ny+2)

  planes_ptr[OLD].data = (double*)malloc( frame_size * sizeof(double) );

  if ( planes_ptr[OLD].data == NULL ) {
    // manage the malloc fail
    fprintf(stderr, "memory allocation of the planes (OLD) failes");
    return 1;
  }

  memset ( planes_ptr[OLD].data, 0, frame_size * sizeof(double) );


  planes_ptr[NEW].data = (double*)malloc( frame_size * sizeof(double) );

  if ( planes_ptr[NEW].data == NULL ) {
    // manage the malloc fail
    fprintf(stderr, "memory allocation of the planes (NEW) failes");

    // since we failed the planes[NEW] allocation, we also need to free 
    // the memory we allocated for the planes[OLD]
    free(planes_ptr[OLD].data);

    return 1;
  }

  memset ( planes_ptr[NEW].data, 0, frame_size * sizeof(double) );


  // ··················································
  // buffers for north and south communication 
  // are not really needed
  //
  // in fact, they are already contiguous, just the
  // first and last line of every rank's plane
  //
  // you may just make some pointers pointing to the
  // correct positions
  //

  // or, if you preer, just go on and allocate buffers
  // also for north and south communications

  // ··················································
  // allocate buffers (only for EAST and WEST)
  // ex. buffers_pt[SEND][WEST] = buffer used to send leftmost column to left neighbour
  // we need one for both SEND and RECEIVE and for each neighbour 
  // N.B neighbours could be MPI_PROC_NULL

  const uint sizex = planes_ptr[OLD].size[_x_];
  const uint sizey = planes_ptr[OLD].size[_y_];

  #define IDX( i, j ) ( (j)*(sizex+2)+ (i) )
  

  // N.B. pointers are already initialized to NULL

  if (neighbours[WEST] != MPI_PROC_NULL) {

    buffers_ptr[SEND][WEST] = (double*)malloc( sizey * sizeof(double) );
    buffers_ptr[RECV][WEST] = (double*)malloc( sizey * sizeof(double) );

    if (!buffers_ptr[SEND][WEST] || !buffers_ptr[RECV][WEST]) {
      fprintf(stderr, "Memory allocation for west buffers failed");

      if ( !buffers_ptr[SEND][WEST] && buffers_ptr[SEND][EAST]) { free(buffers_ptr[SEND][EAST]); }
      else if (!buffers_ptr[RECV][EAST] && buffers_ptr[RECV][WEST]) { free(buffers_ptr[RECV][WEST]); }

      return 1;
    } 
  }

  if (neighbours[EAST] != MPI_PROC_NULL) {
    
    buffers_ptr[SEND][EAST] = (double*)malloc( sizey * sizeof(double) );
    buffers_ptr[RECV][EAST] = (double*)malloc( sizey * sizeof(double) );

    if (!buffers_ptr[SEND][EAST] || !buffers_ptr[RECV][EAST]) {
      fprintf(stderr, "Memory allocation for east buffers failed");

      if (!buffers_ptr[SEND][EAST] && buffers_ptr[SEND][EAST]) { free(buffers_ptr[SEND][EAST]); }
      else if (!buffers_ptr[RECV][EAST] && buffers_ptr[RECV][EAST]) { free(buffers_ptr[RECV][EAST]); }

      return 1;
    } 
  }

  #undef IDX
  return 0;
}


// need to release all the allocated memory (so planes and comm buffers)
int memory_release ( buffers_t *buffer_ptr,
                      plane_t   *planes)
{
  if ( planes != NULL ) {
      if ( planes[OLD].data != NULL ){
	      free(planes[OLD].data);}
      
      if ( planes[NEW].data != NULL ){
	      free(planes[NEW].data);}
    }

  // we free only EAST and WEST buffers, since for NORTH and SOUTH we already points to data
  // in stencil_template_parallel.h -> #define EAST  2 and #define WEST  3

  for (int d=2; d<4; d++){
    if ( buffer_ptr[SEND][d] != NULL ){
      free(buffer_ptr[SEND][d]);}
    
    if ( buffer_ptr[RECV][d] != NULL ){
      free(buffer_ptr[RECV][d]);}
  }
  
  return 0;
}


// function to print the energy (only for rank 0)
int output_energy_stat ( int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm )
{
  double system_energy = 0;
  double tot_system_energy = 0;
  get_total_energy ( plane, &system_energy );
  
  MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
  /* with this: 
  - each rank in the comm gives its local energy
  - the rank 0 stores the global sum in tot_system_energy
  */

  if ( Me == 0 ) {
      if ( step >= 0 ){
	      printf(" [ step %4d ] ", step ); fflush(stdout);}

      printf( "total injected energy is %g, "
	            "system energy is %g "
	            "( in avg %g per grid point)\n",
	            budget,
	            tot_system_energy,
	            tot_system_energy / (plane->size[_x_]*plane->size[_y_]) );}
  
  return 0;
}

int dump ( const double *data, const uint size[2], const char *filename)
{
  if ( (filename != NULL) && (filename[0] != '\0') ){
    FILE *outfile = fopen( filename, "w" );
    if ( outfile == NULL )
      return 2;
    
    float *array = (float*)malloc( size[0] * sizeof(float) );

    for ( int j = 1; j <= size[1]; j++ ) {      
      const double * restrict line = data + j*(size[0] + 2);
      for ( int i = 1; i <= size[0]; i++ ) {
        //int cut = line[i] < 100;
        array[i-1] = (float) line[i];
        
      }
      //printf("\n");
      fwrite( array, sizeof(float), size[0], outfile );
    }

    free( array );

    fclose( outfile );
    return 0;
  }

  return 1;
  
}

