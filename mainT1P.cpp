#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>

#include "P5lib.h"

using namespace std;
double* LocalEnergyT1P(int N_P, int D, double** &pos, double psi_curr, double h_R, double alpha, double beta, double omega,
                    double (*trialFunction)(int, int, double** &, double, double, double) );
// output file
ofstream outfile;

int main(int argc, char* argv[])
  {
    string filename;
    int N_P, D, MCcycles;
    double var_init, var_final, h_var, h_R, omega, fix_value;
    string fix_specify;

    if ( argc > 1 )
      {
        filename    =      argv[1];
        N_P         = atoi(argv[2]);            // Number of particles
        D           = atoi(argv[3]);            // Dimension
        var_init    = atof(argv[4]);
        var_final   = atof(argv[5]);
        h_var       = atof(argv[6]);            // Variational step-size
        fix_specify = argv[7];                  // Specify alpha or beta for fixing
        fix_value   = atof(argv[8]);            // Specify value of fixed parameter
        omega       = atof(argv[9]);
        MCcycles    = pow(10, atoi(argv[10]) );
      }
    // Customize file name with particle number and dimension
    filename = filename.append("NP");
    filename = filename.append(to_string(N_P));
    filename = filename.append("D");
    filename = filename.append(to_string(D));
    filename = filename.append("MC");
    filename = filename.append(argv[10]);

    outfile.open(filename);

    int runNumber = 1;
    int N_StoreValues = 7;    // Number of values to be stored
    double* StoreValues = new double[N_StoreValues];

    for(double var_param = var_init; var_param <= var_final; var_param += h_var)
      {
        // Keeping track of run number
        cout << "Run " << runNumber << endl;

        // Reset expectationvalues
        for(int i = 0; i < N_StoreValues; i++) StoreValues[i] = 0;

        if( fix_specify == "alpha" )
            // Alpha fixed and beta varies
          {
            // Run short Monte Carlo computation to find optimal step-size
            h_R = LocateOptimal_h_R(N_P, D, fix_value, var_param, omega, TrialWaveFunction1);
            cout << "Optimal h_R" << endl;
            cout << h_R << "\n--------------" << endl;

            // Start main Monte Carlo computation
            MetropolisSampling(N_P, D, MCcycles, StoreValues, fix_value, var_param, omega, h_R,
              TrialWaveFunction1, LocalEnergyT1P, E_Local_Analytic_T1P);
            // End main Monte Carlo computation

            // Write expectationvalues to file
            output(N_P, D, MCcycles, StoreValues, fix_value, var_param, omega);

          }
        if( fix_specify == "beta")
            // Beta fixed and alpha varies
          {
            // Run short Monte Carlo computation to find optimal step-size
            h_R = LocateOptimal_h_R(N_P, D, var_param, fix_value, omega, TrialWaveFunction1);
            cout << "Optimal h_R" << endl;
            cout << h_R << "\n--------------" << endl;

            // Start main Monte Carlo computation
            MetropolisSampling(N_P, D, MCcycles, StoreValues, var_param, fix_value, omega, h_R,
              TrialWaveFunction1, LocalEnergyT1P, E_Local_Analytic_T1P);
            // End main Monte Carlo computation

            // Write expectationvalues to file
            output(N_P, D, MCcycles, StoreValues, var_param, fix_value, omega);
          }


        runNumber += 1;
        // End current variational iteration
      }// End variational loop

    // Deallocate memory
    delete[] StoreValues;
    // close output file
    outfile.close();

    return 0;
  }

// Function to calculate the local energy with num derivative
double* LocalEnergyT1P(int N_P, int D, double** &pos, double psi_curr, double h_R, double alpha, double beta, double omega,
                    double (*trialFunction)(int, int, double** &, double, double, double) )
  {
    double* E_Local = new double[2];

    E_Local[0] = Kinetic_Local(N_P, D, pos, psi_curr, h_R, alpha, beta, omega, trialFunction);
    E_Local[1] = U_Harmonic_osc(N_P, D, omega, pos) + 1/sqrt(r_12_squared(N_P, D, pos));

    return E_Local;
  }

void output(int N_P, int D, int MCcycles, double* &StoreValues, double alpha, double beta, double omega)
  {
    double norm = 1.0/((double) (MCcycles)); // divided by number of cycles
    double T_ExpectationValue  = StoreValues[0]*norm;
    double U_ExpectationValue  = StoreValues[1]*norm;
    double E_ExpectationValue  = StoreValues[2]*norm;
    double E2_ExpectationValue = StoreValues[3]*norm;
    double E_analytical        = StoreValues[4]*norm;
    double AcceptedRatio = StoreValues[5]*norm;
    double Dist_r1r2 = StoreValues[6]*norm;
    // all expectation values are per particle
    double Evariance = (E2_ExpectationValue - E_ExpectationValue*E_ExpectationValue);
    outfile << setiosflags(ios::showpoint | ios::uppercase);  //Indices below for readability
    outfile << setw(15) << setprecision(8) << alpha;  //0
    outfile << setw(15) << setprecision(8) << beta;   //1
    outfile << setw(15) << setprecision(8) << omega;  //2
    outfile << setw(15) << setprecision(8) << T_ExpectationValue; //3
    outfile << setw(15) << setprecision(8) << U_ExpectationValue; //4
    outfile << setw(15) << setprecision(8) << E_ExpectationValue; //5
    outfile << setw(15) << setprecision(8) << E_analytical;       //6
    outfile << setw(15) << setprecision(8) << Evariance;          //7
    outfile << setw(15) << setprecision(8) << Dist_r1r2;          //8
    outfile << setw(15) << setprecision(8) << AcceptedRatio;      //9
    outfile << setw(15) << setprecision(8) << MCcycles << endl;   //10
  } // end output function
