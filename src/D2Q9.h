#ifndef D2Q9_H
#define D2Q9_H
/* D2Q9 lattice constants                                        */
/*****************************************************************/

  // opposite directions, for bounce back implementation
static const int oppositeOf[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

  // lattice weights
static const double t[9] = { 4./9., 1./9., 1./9., 1./9., 1./9.,
                             1./36., 1./36., 1./36., 1./36. };
  // lattice velocities
static const int c[9][2] = {
    {0,0},
    {1,0}, {0,1}, {-1,0}, {0,-1},
    {1,1}, {-1,1}, {-1,-1}, {1,-1}
};
#endif