#include "Model.h"

using namespace std;

Model::Model(int argc, char* argv[])
{
    void GetHelp();

    Introduction();
    ParseParameters(argc, argv);
    bool boolian;
    boolian = NumberOfParameters(argc, argv);
    if(boolian == true) {
        PrintParameters();
    }
}

void Model::PrintParameters()
{
    cout << "----------------Input is correct---------------" << endl;
    cout << "----------------Printing Parameters---------------" << endl;
    cout << "Numerics" << endl;
    cout << "Lx= " << Lx << endl;
    cout << "Ly= " << Ly << endl;
    cout << "T= " << T << endl;
    cout << "Nx= " << Nx << endl;
    cout << "Ny= " << Ny << endl;
    cout << "Nt= " << Nt << endl;
    cout << "dx= " << dx << endl;
    cout << "dy= " << dy << endl;
    cout << "dt= " << dt << endl;

    cout << "Physics" << endl;
    cout << "ax= " << ax << endl;
    cout << "ay= " << ay << endl;
    cout << "b= " << b << endl;
    cout << "c= " << c << endl;

    cout << "MPI" << endl;
    cout << "Px= " << Px << endl;
    cout << "Py= " << Py << endl;
}

bool Model::NumberOfParameters(int argc, char* argv[])
{

    bool NoParameters = true;
    string arg = string(argv[1]);
    if(arg == "help")
        NoParameters = false;

    if(argc > 7 & arg != "help") {
        NoParameters = false;
        cout << "Error: Too many arguments!" << endl;
        GetHelp();
    } else if(argc == 7 & arg != "help") {
        NoParameters = true;
    } else if(argc < 7 & arg != "help" ){
        NoParameters = false;
    }

    if(NoParameters == true & arg != "help") {
        for(int i = 1; i < argc; i++) {
            if(isdigit(*argv[i])) {
                NoParameters = true;
            } else {
                NoParameters = false;
                cout << "Error: Please enter digits!" << endl;
                GetHelp();
            }
        }
    }

    return NoParameters;
}

void Model::ParseParameters(int argc, char* argv[])
{
    string arg = string(argv[1]);
    if(arg == "help") {
        GetHelp();
    } else if(argc < 7 & arg != "help") {
        cout << "Error: Too little arguments!" << endl;
        GetHelp();
    } else if(argc >= 7 & arg != "help") {
        ax = stod(argv[1]);      // Convert: sting -> double
        ay = stod(argv[2]);      // Convert: sting -> double
        b = stod(argv[3]);       // Convert: sting -> double
        c = stod(argv[4]);       // Convert: sting -> double
        Px = int(stod(argv[5])); // Convert: sting -> double -> integer
        Py = int(stod(argv[6])); // Convert: sting -> double -> integer
    }
}

void Model::Introduction()
{
    cout << "----------------Introduction----------------" << endl;
    cout << "Hello!" << endl;
    cout << "Welcome to a general Burgers equation solver." << endl;
}

void Model::GetHelp()
{
    cout << "Help is below" << endl;
    cout << "Input order: ax, ay, b, c, Px, Py" << endl;
    cout << "Data type: double, double, double, double, int, int" << endl;
}