#include <chrono>
#include <stdio.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
using namespace std;
#include "Model.h"
#include "burgers.h"

int main(int argc, char* argv[]) {
    // if "-h" is present, run help function
    
    Model m(argc, argv);
    
    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();
    burgers b(&m);
    MPI_Init(&argc, &argv);

   // Call code to initialise the problem here
    b.Initial_velocity();
    // Call code to perform time integration here
    b.Integrate_velocity();
    // Print velocity matrix 
    //b.print_velocity();
    // Calculate final energy
    b.energy_field();

    hrc::time_point end = hrc::now();
    
    auto time_spent = end - start;
    cout << "Time past: " << chrono::duration <double, milli> (time_spent).count() << "ms" << endl; 
    
    MPI_Finalize();
    return 0;
}