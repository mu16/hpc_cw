#include "burgers.h"
#include "Model.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

burgers::burgers(Model* m_)
{
    m = m_;
    Nx = m->GetNx();
    Ny = m->GetNy();
    ax = m->GetAx();
    ay = m->GetAy();
    b = m->GetB();
    c = m->GetC();
    Nt = m->GetNt();
    dx = m->GetDx();
    dy = m->GetDy();
    dt = m->GetDt();
    T = m->GetT();
    x0 = m->GetX0();
    y0 = m->GetY0();
    Lx = m->GetLx(); 
    Ly = m->GetLy(); 
    Px = m->GetPx(); 
    Py = m->GetPy(); 
}

void burgers::Initial_velocity()
{

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    NxSub = Nx / Px;
    NySub = Ny / Py;

    NxRem = Nx % Px;
    NyRem = Ny % Py;

    row = myrank / Px + 1;
    col = myrank % Px + 1;

    if(NyRem == 0 && NxRem == 0) {
        positionX = (col - 1) * NxSub;
        positionY = (row - 1) * NySub;
    } else {
        if(row <= NyRem) {
            positionY = (row - 1) * NySub + (row - 1);
        } else if (row>NyRem) {
            positionY = (row - 1) * NySub + NyRem;
        }
        
        if(col <= NxRem) {
            positionX = (col - 1) * NxSub + (col - 1);
        } else if (col > NxRem){
            positionX = (col - 1) * NxSub + NxRem;
        }
    }

    if(col <= NxRem) NxSub++;

    if(row <= NyRem) NySub++;

    // Initilizing the arrays for the space and velocity
    double x[NxSub];
    double y[NySub];
    u = new double[NxSub * NySub];
    v = new double[NxSub * NySub];

    // initilizing space matrices
    for(int i = 0; i < NxSub; i++) x[i] = x0  -Lx / 2 + dx * (i + positionX);

    for(int i = 0; i < NySub; i++) y[i] = y0  -Ly / 2 + dy * (i+positionY);

    double r;
    // Initial Conditions for velocity at t=0
    for(int i = 0; i < NxSub; i++) {
        for(int j = 0; j < NySub; j++) {
            r = sqrt((x[i] * x[i] + y[j] * y[j]));

            if(r <= 1) {
                u[i * NySub + j] = 2.0 * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (4.0 * r + 1.0);
                v[i * NySub + j] = 2.0 * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (4.0 * r + 1.0);
            } else if(r>1) {
                u[i * NySub + j] = 0;
                v[i * NySub + j] = 0;
            }
        }
    }
 // Left Top Corner
} // End initial_velocity

void burgers::Integrate_velocity()
{
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    // Initilising values
    int lowerBoundX=0;
    int upperBoundX=0;
    int lowerBoundY=0;
    int upperBoundY=0;
    
    int flagLeftTop=0;
    int flagRightTop=0;
    int flagLeftBot=0;
    int flagRightBot=0;
    int flagTopMid=0;
    int flagBotMid=0;
    int flagRightMid=0;
    int flagLeftMid=0;
    int flagOrdinary=0;
    
    int jPlusFlag=0;
    int jMinusFlag=0;
    int iPlusFlag=0;
    int iMinusFlag=0;
    
    auto* unew = new double[NxSub*NySub];
    auto* vnew = new double[NxSub*NySub];
    double unext;
    double vnext;
    
    uLeftMine = new double[NySub];
    vLeftMine = new double[NySub];
    
    uUpMine = new double[NxSub];
    vUpMine = new double[NxSub];
    
    uDownMine = new double[NxSub];
    vDownMine = new double[NxSub];
    
    uRightMine = new double[NySub];
    vRightMine = new double[NySub];
    
    uLeftUrs = new double[NySub];
    vLeftUrs = new double[NySub];

    uRightUrs = new double[NySub];
    vRightUrs = new double[NySub];

    uUpUrs = new double[NxSub]; 
    vUpUrs = new double[NxSub]; 

    uDownUrs = new double[NxSub]; 
    vDownUrs = new double[NxSub];
    
    double iMinusTermU;
    double iPlusTermU;
    double jPlusTermU;
    double jMinusTermU;
    
    double iMinusTermV;
    double iPlusTermV;
    double jPlusTermV;
    double jMinusTermV;
      
    // Coefficinets for forward time integration //
    // coefiMinus: i-1 term
    // coefiPlus: i+1 term
    // coef : i,j term
    // coefjPlus: j+1 term
    // coefjMinus: j-1 term
    double coefiMinus =  ax / dx + c / dx / dx;
    double coefiPlus = c / (dx * dx);
    double coef  = 1/dt - 2.0*c*(1/(dx * dx) + 1/(dy * dy)) - ax/dx - ay/dy;
    double coefjPlus =  ay / dy + c / (dy * dy);
    double coefjMinus = c / (dy * dy);
    
    // Looking at the whole array, the code below figures ou where which a curent rank is //
    // flag is equal to 1 if it is true, sort of like boolian
    if(col == 1 && row == Py) flagLeftBot=1;
    if( col != 1 && col != Px && row == 1) flagTopMid=1;
    if(col == Px && row != 1 && row != Py) flagRightMid=1;
    if (col!=1 && col!=Px && row!=1 && col!=Px) flagOrdinary=1;
    if(col == 1 && row!= 1 && row != Py) flagLeftMid=1;   
    if(col != 1 && col!= Px && row == Py) flagBotMid=1;
    if(col == Px && row == Py) flagRightBot=1;
    if(col == Px && row == 1) flagRightTop=1;
    if(col == 1 && row == 1) flagLeftTop=1;
        
        
    // After determining the location, the bounds are then computed for each case
    // The bounds for x and y form a smaller rectangle inside the whole domain
    // The size of the smaller rectangle is smaller by 1 in each direction
    if (flagLeftBot==1){
        lowerBoundX = 1;
        upperBoundX= NxSub;
        
        lowerBoundY = 0;
        upperBoundY = NySub - 1;
    }
    
    if (flagTopMid==1){
        lowerBoundX = 0;
        upperBoundX = NxSub;
        
        lowerBoundY= 1;
        upperBoundY = NySub;
    } 
    
    if (flagRightMid==1){
        lowerBoundX = 0;
        upperBoundX = NxSub-1;
        
        lowerBoundY = 0;
        upperBoundY = NySub;
    }
    
      if (flagOrdinary==1){
        lowerBoundX = 0;
        upperBoundX = NxSub;
        
        lowerBoundY = 0;
        upperBoundY = NySub;
    }
        
        
    if (flagLeftTop==1){
        lowerBoundX = 1;
        upperBoundX = NxSub;
        
        lowerBoundY = 1;
        upperBoundY = NySub;
    }

    if (flagRightTop==1){
        lowerBoundX = 0;
        upperBoundX  = NxSub - 1;

        lowerBoundY  = 1;
        upperBoundY = NySub;
    }

    if (flagRightBot==1){
        lowerBoundX = 0;
        upperBoundX= NxSub - 1;
        
        lowerBoundY= 0;
        upperBoundY = NySub - 1;
    }

    if (flagBotMid==1){
        lowerBoundX = 0;
        upperBoundX = NxSub;
        
        lowerBoundY = 0;
        upperBoundY = NySub - 1;
    }  
    
    if (flagLeftMid==1){
        lowerBoundX = 1;
        upperBoundX = NxSub;
        
        lowerBoundY  = 0;
        upperBoundY = NySub;
    }
    
    
    // Time integration
    for (int timeStep=1; timeStep<Nt; timeStep++)
        {
            exchange(); // exchanging values
                for (int j=lowerBoundY; j < upperBoundY; j++)                    
                    {
                    for (int i = lowerBoundX; i < upperBoundX; i++)
                        {
                            
                            if (i==0 || i==NxSub-1 || j==0 || j==NySub-1){
                                
                            if (i==0 && j==0){
            
                            iMinusTermU = uLeftMine[j]; 
                            iMinusTermV = vLeftMine[j];
                            
                            jMinusTermU=uUpMine[i];
                            jMinusTermV=vUpMine[i];
                            
                            iPlusTermU=u[(i+1)*NySub+j];
                            iPlusTermV=v[(i+1)*NySub+j];
                            
                            jPlusTermU=u[i*NySub+(j+1)]; 
                            jPlusTermV=v[i*NySub+(j+1)];
                            
                            
                        } else if (i==NxSub-1 && j==0){

                            iPlusTermU=uRightMine[j];
                            iPlusTermV=vRightMine[j];
                            
                            jMinusTermU=uUpMine[i];
                            jMinusTermV=vUpMine[i];
                            
                            iMinusTermU=u[(i-1)*NySub+j];
                            iMinusTermV=v[(i-1)*NySub+j];
                            
                            jPlusTermU=u[i*NySub+(j+1)]; 
                            jPlusTermV=v[i*NySub+(j+1)];
                            
                            
                        } else if (i==NxSub-1 && j==NySub-1){

                            iPlusTermU=uRightMine[j];
                            iPlusTermV=vRightMine[j];
                            
                            jPlusTermU=uDownMine[i];
                            jPlusTermV=vDownMine[i];
                            
                            iMinusTermU=u[(i-1)*NySub+j];
                            iMinusTermV=v[(i-1)*NySub+j];
                            
                            jMinusTermU=u[i*NySub+(j-1)];
                            jMinusTermV=u[i*NySub+(j-1)];
                            
                        }  else if (i==0 && j==NySub-1 ){

                            jPlusTermU=uDownMine[i];
                            jPlusTermV=vDownMine[i];
                            
                            iMinusTermU=uLeftMine[j]; 
                            iMinusTermV=vLeftMine[j];
                            
                            iPlusTermU=u[(i+1)*NySub+j];
                            iPlusTermV=v[(i+1)*NySub+j];
                            
                            jMinusTermU=u[i*NySub+(j-1)];
                            jMinusTermV=u[i*NySub+(j-1)];
                        

                        } else if(j == 0 && i!=0 && i!=NxSub-1){
                            
                            jMinusTermU=uUpMine[i];
                            jMinusTermV=vUpMine[i];
                            
                            jPlusTermV=v[i*NySub+(j+1)];
                            jPlusTermU=u[i*NySub+(j+1)];
                            
                            iMinusTermU=u[(i-1)*NySub+j];
                            iMinusTermV=v[(i-1)*NySub+j];
                            
                            iPlusTermU=u[(i+1)*NySub+j];
                            iPlusTermV=v[(i+1)*NySub+j];
                            
                        } else if (i==0 && j!=0 && j!=NySub-1){
                            
                            iMinusTermU=uLeftMine[j]; 
                            iMinusTermV=vLeftMine[j];
                            
                            jPlusTermV=v[i*NySub+(j+1)];
                            jPlusTermU=u[i*NySub+(j+1)];
                            
                            iPlusTermU=u[(i+1)*NySub+j];
                            iPlusTermV=v[(i+1)*NySub+j];
                        
                            jMinusTermU=u[i*NySub+(j-1)];
                            jMinusTermV=v[i*NySub+(j-1)];
                            
                        } else if(i!=0 && j == NySub-1 && i!=NySub-1){

                            jPlusTermU=uDownMine[i];
                            jPlusTermV=vDownMine[i];
                            
                            iMinusTermU=u[(i-1)*NySub+j];
                            iMinusTermV=v[(i-1)*NySub+j];
                            
                            iPlusTermU=u[(i+1)*NySub+j];
                            iPlusTermV=v[(i+1)*NySub+j];
                        
                            jMinusTermU=u[i*NySub+(j-1)];
                            jMinusTermV=v[i*NySub+(j-1)];
                            
                        } else if (  j!=0 && i==NxSub-1 &&  j!=NySub){

                            iPlusTermU=uRightMine[j];
                            iPlusTermV=vRightMine[j];
                            
                            jPlusTermV=v[i*NySub+(j+1)];
                            jPlusTermU=u[i*NySub+(j+1)];
                            
                            iMinusTermU=u[(i-1)*NySub+j];
                            iMinusTermV=v[(i-1)*NySub+j];
                            
                            jMinusTermU=u[i*NySub+(j-1)];
                            jMinusTermV=v[i*NySub+(j-1)];
                        }
                            
                        } else {
                            
                            jPlusTermV=v[i*NySub+(j+1)];
                            jPlusTermU=u[i*NySub+(j+1)];
                            
                            iMinusTermU=u[(i-1)*NySub+j];
                            iMinusTermV=v[(i-1)*NySub+j];
                            
                            iPlusTermU=u[(i+1)*NySub+j];
                            iPlusTermV=v[(i+1)*NySub+j];
                        
                            jMinusTermU=u[i*NySub+(j-1)];
                            jMinusTermV=v[i*NySub+(j-1)];
                            
                        }                        
                        
                            unew[i*NySub+j]=dt* (coefiMinus  * iMinusTermU + coefiPlus  * iPlusTermU+ coef * u[i*NySub+j] + coefjMinus  * jMinusTermU +  coefjPlus * jPlusTermU + b/dx * u[i*NySub+j] * (-u[i*NySub+j]+iMinusTermU) + b/dy * v[i*NySub+j] * (jPlusTermU-u[i*NySub+j] )); 
                            vnew[i*NySub+j]=dt*(coefiMinus  * iMinusTermV + coefiPlus  * iPlusTermV + coef * v[i*NySub+j] + coefjMinus  * jMinusTermV + coefjPlus * jPlusTermV + b/dx * v[i*NySub+j] * (-v[i*NySub+j]+jPlusTermV) + b/dy * v[i*NySub+j] * (iMinusTermV-v[i*NySub+j] )); 
                            
                            }
                        }
                    
        // Reassigning values to the velocity arrays
        for(int jnew=lowerBoundY; jnew < upperBoundY; jnew++)
            {
            for (int inew=lowerBoundX; inew<upperBoundX; inew++)
                {
                u[inew*NySub+jnew]=unew[inew*NySub+jnew];
                v[inew*NySub+jnew]=vnew[inew*NySub+jnew];
                }
            }

            // Outputting time step
            if(myrank == 0) {
                if(timeStep % 25 == 0) {
                    cout << "Nt is currently at: " << timeStep << endl;
                }
            }
    }
        } // End integration


void burgers::print_velocity()
{
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    string filenameu=std::to_string(myrank) + "u" + ".txt";
    ofstream uprint(filenameu, fstream::trunc);
    for (int i=0; i < NxSub; i++){
        for (int j=0; j<NySub;j++){
            uprint.precision(5);
            uprint<< fixed << u[i*NySub+j] << " ";
        }
    }
    uprint.close();
    
    string filenamev=std::to_string(myrank) + "v" + ".txt";
    ofstream vprint(filenamev, fstream::trunc);
    for (int i=0; i < NxSub; i++){
        for (int j=0; j<NySub;j++){
            vprint.precision(5);
            vprint<< fixed << v[i*NySub+j] << " ";
        }
    }
    vprint.close();

}

void burgers::energy_field()
{
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        
        // Initialising values
        double energy;
        double energyRecv;
        
        // Performing energy integration
        for (int i=0; i<NxSub*NySub ;i++)
        {
            energy+=(u[i]*u[i]+v[i]*v[i]);
        }
        
        // Every rank except rank 0 sends data to rank 0
        if (myrank!=0){
        MPI_Send(&energy,1,MPI_DOUBLE, 0, 123, MPI_COMM_WORLD); // passes energy
        cout << "Sent energy to rank: " << myrank << endl;
        } 
        
        // Rank 0 performs the calculation of the energy
        if (myrank==0){
            for (int ranking = 1; ranking < Px*Py; ranking++){
                MPI_Recv(&energyRecv,1, MPI_DOUBLE, ranking, 123,  MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                energy += energyRecv;
            }
        }
        
        // Final calculation by rank 0
        if (myrank==0){
            energy *= 0.5 * dx *dy;
            cout<< "Final Energy: " << energy <<endl;
        }
        
} // End energy_field

void burgers::exchange(){
    
    // This function sends and recieves u and v data with neighbouring arrays
    
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    // Right side of the matrix to be sent
    for(int j = 0; j < NySub; j++) {
        uRightUrs[j]=u[(NxSub-1)*NySub+j];
        vRightUrs[j]=v[(NxSub-1)*NySub+j];
    }
    
    // Left side of the matrix to be sent
    for(int j = 0; j < NySub; j++) {
        uLeftUrs[j]=u[0*NySub+j];
        vLeftUrs[j]=v[0*NySub+j];
    }

    // Down side of the matrix to be sent
    for(int i = 0; i < NxSub; i++) {
        uDownUrs[i]=u[i*NySub+NySub-1]; 
        vDownUrs[i]=v[i*NySub+NySub-1];
    }
    
    // Upper side of the matrix to be sent
        for(int i = 0; i < NxSub; i++) {
        uUpUrs[i]=u[i*NySub+0];
        vUpUrs[i]=v[i*NySub+0];
    }
    
    // Finding the rank by a position in the whole domain
    // left - left
    // right - right
    // down - one below
    // up - one above
    
    left = myrank - 1;
    right = myrank + 1;
    down = myrank + Px;
    up = myrank - Px;
    
    // Check if there is anyhting to the left
    if(col - 1 > 0) {
        MPI_Recv(uLeftMine, NySub, MPI_DOUBLE, left, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(vLeftMine, NySub, MPI_DOUBLE, left, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(uLeftUrs, NySub, MPI_DOUBLE, left, 300, MPI_COMM_WORLD);
        MPI_Send(vLeftUrs, NySub, MPI_DOUBLE, left, 400, MPI_COMM_WORLD);
    }

    // Check if there is anyhting to the right
    if(col+ 1 <= Px) {
        MPI_Send(uRightUrs, NySub, MPI_DOUBLE, right, 100, MPI_COMM_WORLD);
        MPI_Send(vRightUrs, NySub, MPI_DOUBLE, right, 200, MPI_COMM_WORLD);

        MPI_Recv(uRightMine, NySub, MPI_DOUBLE, right, 300, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(vRightMine, NySub, MPI_DOUBLE, right, 400, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Check if there is anyhting below
    if(row + 1 <= Py) {
        MPI_Send(uDownUrs, NxSub, MPI_DOUBLE, down, 1300, MPI_COMM_WORLD);
        MPI_Send(uDownUrs, NxSub, MPI_DOUBLE, down, 1400, MPI_COMM_WORLD);

        MPI_Recv(uDownMine, NxSub, MPI_DOUBLE, down, 11100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(vDownMine, NxSub, MPI_DOUBLE, down, 12100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    // Check if there is anyhting above
    if(row - 1 > 0) {
        MPI_Recv(uUpMine, NxSub, MPI_DOUBLE, up, 1300, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(vUpMine, NxSub, MPI_DOUBLE, up, 1400, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Send(uUpUrs, NxSub, MPI_DOUBLE, up, 11100, MPI_COMM_WORLD);
        MPI_Send(vUpUrs, NxSub, MPI_DOUBLE, up, 12100, MPI_COMM_WORLD);
    }
    
}
    

 