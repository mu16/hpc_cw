#ifndef BURGERS_H
#define BURGERS_H
#include <mpi.h>
#include "Model.h"
#include <math.h>  

class burgers
{
public:
    burgers(Model* m);
    //   ~burgers();
    
    // Declaring functions
    void Initial_velocity();
    void Integrate_velocity();
    void energy_field();
    void exchange();
    void print_velocity();
    
private:

    Model* m;
    
    // Numerics
    double x;
    double y;
    double x0;
    double y0;
    double Lx;
    double Ly;
    double T;
    int  Nx;
    int  Ny;
    int Nt;
    double dx;
    double dy;
    double dt;
    
    // Physics
    double ax;
    double ay;
    double b;
    double c;
    
    // Velocity field arrays
    double*u;
    double*v;
    double*un;
    double*vn;
    
    // MPI
    int myrank;
    int down;
    int up;
    int left;
    int right;
    int Px;
    int Py;
    int NxSub;
    int NySub;
    int NyRem;
    int NxRem;
    double* uLeftMine;
    double* vLeftMine;
    double* uUpMine;
    double* vUpMine;
    double* uRightMine;
    double* vRightMine;
    double* uDownMine;
    double* vDownMine;
    
    double* uLeftUrs;
    double* vLeftUrs;
    double* uUpUrs;
    double* vUpUrs;
    double* uRightUrs;
    double* vRightUrs;
    double* uDownUrs;
    double* vDownUrs;
    
    int positionX;
    int positionY;
    int row;
    int col;
    
};

#endif




