/* lb.h: Assembly of various functions for simulating a lattice
 *       Boltzmann dynamics with BGK collision term on a two-
 *       dimenional D2Q9 lattice. The code is generic: every
 *       lattice site can have a different dynamics.
 *       Other collision terms than BGK can easily be added.
 */


#ifndef LB_H
#define LB_H

#include "stddef.h"

/*------------------------add by Yuankun------------------------*/
#include <stdint.h>
#include <assert.h>
#include <errno.h>
#include <string.h>

#define clean_errno() (errno == 0 ? "None" : strerror(errno))
#define log_error(M, ...) fprintf(stderr, "[ERROR] (%s:%d: errno: %s) " M "\n", __FILE__, __LINE__, clean_errno(), ##__VA_ARGS__)
#define assertf(A, M, ...) if(!(A)) {log_error(M, ##__VA_ARGS__); assert(A); }
#define ERROR_RETURN(retval) { fprintf(stderr, "Error %d %s:line %d: \n", retval,__FILE__,__LINE__);  exit(retval); }

// #define ZGB

//#define SAVE

#ifdef _OPENMP
  #include <omp.h>
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

extern int tile;
extern double *myrho1, *myrho2;
extern int my_domain_H; // each thread's data domain's height
extern int NUM_THREADS;
/*------------------------End add by Yuankun------------------------*/

/* struct Dynamics                                               */
/*****************************************************************/
/*   emulation of a class that defines the dynamics of a lattice
 *   site. The function dynamicsFun contains the algorithm of the
 *   collision, and selfData points to the function arguments
 *   (for example the local viscosity)
 */
typedef struct {
    void   (*dynamicsFun) (double* fPop, void* selfData);
    void* selfData;
} Dynamics;


/* struct Node                                                   */
/*****************************************************************/
/*  a lattice node, containing the data (distribution functions)
 *  for a D2Q9 lattice, and a pointer to the local dynamics, i.e.
 *  the collision term. Two "methods" are added to construct and
 *  initialize a node.
 */
typedef struct {
    double    fPop[9];
    Dynamics* dynamics;
} Node;

void constructNode(Node* node);
void iniEquilibrium(Node* node, double rho, double ux, double uy);


/* struct Simulation                                             */
/*****************************************************************/
/*  a full D2Q9 LB simulation, containing the allocated memory
 *  for the lattice. Some "methods" are added to initiate the
 *  dynamics.
 */
typedef struct {
    int lx, ly;               // lx*ly lattice
    Node*  memoryChunk;       // contiguous raw memory
    Node*  tmpMemoryChunk;    // contiguous raw tmp memory
    Node** lattice;           // lattice, points to raw memory
    Node** tmpLattice;        // tmp lasttice, points to raw memory
} Simulation;

/*---------------- Global variables -----------------*/
extern int lx, ly;
extern int obst_x, obst_y, obst_r;
extern double uMax, Re, nu, omega;
extern int numIter, tSave, warmUpIter;

extern Dynamics       bulkDynamics;
extern Dynamics       bounceBackDynamics;

extern Dynamics upperBoundary;
extern Dynamics lowerBoundary;
extern Dynamics* leftBoundary;   // Those two objects are initialized in the
extern Dynamics* rightBoundary;  //   function iniData()

  // The main object containing the simulation
extern Simulation sim;
/*---------------- End original Global variables -----------------*/

void constructSim(Simulation* sim, int lx, int ly);
void destructSim(Simulation* sim);
void setDynamics(Simulation* sim, int iX, int iY, Dynamics* dyn);
void collide(Simulation* sim);
void propagate(Simulation* sim);
void makePeriodic(Simulation* sim);
void saveVel(Simulation* sim, char fName[]);
void saveF(Simulation* sim, int iPop, char fName[]);

/*------------------------add by Yuankun------------------------*/
extern inline void collideNode(Node* node);
extern inline void collide_stream_buf1_to_buf2(Simulation* sim, int iX, int iY);
extern inline void collide_stream_buf2_to_buf1(Simulation* sim, int iX, int iY);

void collideOMP(Simulation* sim);
void propagateOMP(Simulation* sim);

void collideStream(Simulation* sim);
void collideStreamOMP(Simulation* sim);
void swapLattice(Simulation* sim);

void collideStreamTile(Simulation* sim);
void collideStreamTileOMP(Simulation* sim);

void updateZeroGradientBoundary();
void step2CollideStream(Simulation* sim);
void step2CollideStreamOMP(Simulation* sim);

void step2CollideStreamTile(Simulation* sim);
void step2CollideStreamTileOMP(Simulation* sim);

void step3CollideStream(Simulation* sim);
void step3CollideStreamOMP(Simulation* sim);

void step3CollideStreamTile(Simulation* sim);
void step3CollideStreamTileOMP(Simulation* sim);
/*------------------------End add by Yuankun------------------------*/

/* some free helper functions                                    */
/*****************************************************************/

  // compute density and velocity from the f's 
void computeMacros(double* f, double* rho, double* ux, double* uy);

  // compute local equilibrium from rho and u
double computeEquilibrium(int iPop, double rho,
                          double ux, double uy, double uSqr);
  // bgk collision term
void bgk(double* fPop, void* selfData);


// init
void setConstants(int argc, char *argv[]);
void iniData();
void freeData();
void iniGeometry();
void Hello();
#endif
