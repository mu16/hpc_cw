#ifndef CLASS_MODEL
#define CLASS_MODEL
#include <iostream>
#include <cctype>
// #include <stdexcept>
//#include <mpi.h>


using namespace std;

class Model {
public:
    
    Model(int argc, char* argv[]);
    //void PrintParameters();
    //bool IsValid(int argc, char* argv[]);
    //~Model(){};

    // Getters
    bool   IsVerbose() const { return verbose; }
    bool   IsHelp()    const { return help; }
    double GetX0()     const { return x0; }
    double GetY0(   )     const { return y0; }
    double GetLx()     const { return Lx; }
    double GetLy()     const { return Ly; }
    double GetT()      const { return T; }
    int    GetNx()     const { return Nx; }
    int    GetNy()     const { return Ny; }
    int    GetNt()     const { return Nt; }
    double GetDx()     const { return dx; }
    double GetDy()     const { return dy; }
    double GetDt()     const { return dt; }
    double GetAx()     const { return ax; }
    double GetAy()     const { return ay; }
    double GetB()      const { return b; }
    double GetC()      const { return c; }
    int GetPx()      const { return Px; }
    int GetPy()      const { return Py; }
    
private:
    
    void Introduction();
    void GetHelp();
    void PrintParameters();
    void ParseParameters(int argc, char* argv[]);
    bool NumberOfParameters(int argc, char* argv[]);
    
    bool verbose;
    bool help;
    
    // Numerics
    double x0;
    double y0;
    double Lx=10;
    double Ly=10;
    double T=1;
    int    Nx=2001;
    int    Ny=2001;
    int    Nt=4000;
    double dx=Lx/(Nx-1);
    double dy=Ly/(Ny-1);
    double dt=T/Nt;
    
    // Physics
    double ax;
    double ay;
    double b;
    double c;
    
    // MPI
    int Px;
    int Py;
    
};

#endif