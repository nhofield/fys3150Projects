/*Function library*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
using namespace std;

#include "P5lib.h"
/*--------------------------------
|    Main-Algorithm functions    |
/*------------------------------*/
void MetropolisSampling(int N_P, int D, int MCcycles, double* &StoreValues, double alpha, double beta, double omega, double h_R,
                        double (*trialFunction)(int, int, double** &, double, double, double),
                        double* (*localEnergy) (int, int, double** &, double, double, double, double, double,
                        double (*trialFunction)(int, int, double** &, double, double, double)),
                        double (*localEnergyAnalytic)(int, int, double** &, double, double, double) )
  {
    //Prepare mersenne twister
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> positionDistribution(-1.0, 1.0);  // When assigning initial positions
    std::uniform_real_distribution<double> uniformDistribution(0, 1);        // When performing metropolis test
    std::uniform_real_distribution<double> moveDistribution(-0.5, 0.5);      // When performing translation of particles

    // Allocate memory for old and new positions matrices
    double** r_curr = MatrixAlloc(N_P, D);
    double** r_new  = MatrixAlloc(N_P, D);

    double psi_new, psi_curr;     // To keep tabs of wavefunction evaluation

    for(int i = 0; i < N_P; i++)
      {
        for(int j = 0; j < D; j++) r_new[i][j] = r_curr[i][j] = 0;
      }

    // initialize quantities
    double* E_Local_arr = new double[2];                  // Array for unloading kinetic and potential energy

    double T_Energy     = 0; double U_Energy = 0;         // For expectation value of kinetic and potential energy
    double Energy       = 0; double Energy2  = 0;         // For expectation value of total energy and total energy squared
    double E_analytical = 0; double Delta_E  = 0;         // For analytical local energy and temporary storage of energy
    double Delta_E_analytical = 0; double Dist_r1r2 = 0;  // For temporary storage of analytical total energy, and distance between particles

    double Probability_ratio = 0;

    int Accepted = 0;                                     //To store number of accepted transitions

    // Set initial positions of particles randomly within interval [-a,a], a real number.
    for (int i = 0; i < N_P; i++)
      {
        for (int j = 0; j < D; j++) r_curr[i][j] = positionDistribution(gen);
      }

    psi_curr = trialFunction(N_P, D, r_curr, alpha, beta, omega);

    // Begin Monte Carlo cycle
    for (int cycle = 1; cycle <= MCcycles; cycle++)
      {
        // Loop over particles and dimensions
        for(int particle = 0; particle < N_P; particle++)
          {
            for(int dimension = 0; dimension < D; dimension++)
              {
                // Translate particles randomly
                r_new[particle][dimension] = r_curr[particle][dimension] + moveDistribution(gen) * h_R;
              }
          }
        // Evaluate the trial wavefunction at the new positions
        psi_new = trialFunction(N_P, D, r_new, alpha, beta, omega);


        Probability_ratio = (psi_new*psi_new)/(psi_curr*psi_curr);

        // Metropolis test of probability ratio
        if( Probability_ratio >= uniformDistribution(gen) )
        // Accept
          {
            // update quantities
            for(int i = 0; i < N_P; i++)
              {
                for(int j = 0; j < D; j++) r_curr[i][j] = r_new[i][j];
              }
            psi_curr = psi_new;
            // Compute local energy and store kinetic energy as first element and potential as second
            E_Local_arr = localEnergy(N_P, D, r_curr, psi_curr, h_R, alpha, beta, omega, trialFunction);
            // Store local energy in Delta_E
            Delta_E = E_Local_arr[0] + E_Local_arr[1];
            // Compute analytical local energy
            Delta_E_analytical = localEnergyAnalytic(N_P, D, r_curr, alpha, beta, omega);

            T_Energy     += E_Local_arr[0];    U_Energy += E_Local_arr[1];
            Energy       += Delta_E;            Energy2 += Delta_E*Delta_E;
            E_analytical += Delta_E_analytical;

            Accepted  += 1;
            Dist_r1r2 += sqrt( r_12_squared(N_P, D, r_curr) );
          }

       else
       // Keep position, but add contributions to update
        {
          // Compute local energy and store kinetic energy as first element and potential as second
          E_Local_arr = localEnergy(N_P, D, r_curr, psi_curr, h_R, alpha, beta, omega, trialFunction);
          // Store local energy in Delta_E
          Delta_E = E_Local_arr[0] + E_Local_arr[1];
          // Compute analytical local energy
          Delta_E_analytical = localEnergyAnalytic(N_P, D, r_curr, alpha, beta, omega);

          T_Energy     += E_Local_arr[0];    U_Energy += E_Local_arr[1];
          Energy       += Delta_E;            Energy2 += Delta_E*Delta_E;
          E_analytical += Delta_E_analytical;

          Dist_r1r2 += sqrt( r_12_squared(N_P, D, r_curr) );        }

      }
    // update expectation values for local node
    StoreValues[0] = T_Energy;      StoreValues[1] = U_Energy;
    StoreValues[2] = Energy;        StoreValues[3] = Energy2;
    StoreValues[4] = E_analytical;  StoreValues[5] = Accepted;
    StoreValues[6] = Dist_r1r2;
    //Deallocate memory
    MatrixDeAlloc(r_curr, N_P);
    MatrixDeAlloc(r_new, N_P);
  } // end of Metropolis sampling function

double LocateOptimal_h_R(int N_P, int D, double alpha, double beta, double omega,
  double (*trialFunction)(int, int, double** &, double, double, double) )
  {
    //Prepare mersenne twister
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> positionDistribution(-1.0, 1.0);
    std::uniform_real_distribution<double> uniformDistribution(0, 1);
    std::uniform_real_distribution<double> moveDistribution(-0.5, 0.5);

    double h_R = 15;
    int MC = 1000;
    double epsilon = 0.05;    //Five percent
    int Accepted = 0; //To store number of accepted transitions
    double AcceptRatio = Accepted/( (double) MC );
    double Probability_ratio = 0;
    double psi_new, psi_curr;
    bool condition = (0.5 - epsilon) < AcceptRatio && AcceptRatio < (0.5 + epsilon);
    int runNumber = 0;

    while( condition != 1 && h_R > 0 )
    {
      // Allocate memory for old and new positions matrices
      double** r_curr = MatrixAlloc(N_P, D);
      double** r_new  = MatrixAlloc(N_P, D);

      for(int i = 0; i < N_P; i++)
        {
          for(int j = 0; j < D; j++) r_new[i][j] = r_curr[i][j] = 0;
        }

      // Set initial positions of particles randomly within [-2,2]
      for (int i = 0; i < N_P; i++)
        {
          for (int j = 0; j < D; j++) r_curr[i][j] = positionDistribution(gen);
        }

      psi_curr = trialFunction(N_P, D, r_curr, alpha, beta, omega);

    // Begin Monte Carlo cycle
    for (int cycle = 1; cycle <= MC; cycle++)
      {
        // Loop over particles and dimensions
        for(int particle = 0; particle < N_P; particle++)
          {
            for(int dimension = 0; dimension < D; dimension++)
              {
                // Translating particles randomly
                r_new[particle][dimension] = r_curr[particle][dimension] + moveDistribution(gen) * h_R;
              }
          }
        // Evaluate the trial wavefunction at the new positions
        psi_new = trialFunction(N_P, D, r_new, alpha, beta, omega);

        // Metropolis test of probability ratio
        Probability_ratio = (psi_new*psi_new)/(psi_curr*psi_curr);

        if(Probability_ratio >= uniformDistribution(gen))
        // Accept
          {
            // update energies and position
            for(int i = 0; i < N_P; i++)
              {
                for(int j = 0; j < D; j++) r_curr[i][j] = r_new[i][j];
              }
            Accepted += 1;
            psi_curr = psi_new;
          }

        else
        // keep position
          {

          }

      }

    MatrixDeAlloc(r_curr, N_P);
    MatrixDeAlloc(r_new, N_P);
    AcceptRatio = Accepted/( (double) MC );
    h_R -= 0.001;
    Accepted = 0;
    condition = AcceptRatio > 0.5 + epsilon;
    }
  return h_R;

  }



/*---------------------------------------------
|   Energy-specific computational functions   |
/*-------------------------------------------*/
double U_Harmonic_osc(int N_P, int D, double omega, double** &pos)
    // harmonic oscillator potential
  {
    return 0.5*omega*omega*sum_r1r2_squared(N_P, D, pos);
  }

double U_ee_Repulsion(int N_P, int D, double** &pos)
    // harmonic oscillator potential
  {
    double r_12 = 0; double U_ee = 0;
    for (int i = 0; i < N_P - 1; i++)
      {
        for (int j = i + 1; j < N_P; j++)
          {
            r_12 = 0;
            for (int k = 0; k < D; k++)
              {
                r_12 += ( pos[i][k] - pos[j][k] )*( pos[i][k]-pos[j][k] );
              }
            U_ee += 1/sqrt(r_12);
          }
    }
  }


double Kinetic_Local(int N_P, int D, double** &pos, double psi_curr, double h_R, double alpha, double beta, double omega,
            double (*trialFunction)(int, int, double** &, double, double, double) )

  {
    int i, j;
    double psi_minus, psi_plus, T_Local;
    double hh_R = h_R*h_R;

    // allocate matrices which contain the position of the particles
    double** r_plus  = MatrixAlloc(N_P, D);
    double** r_minus = MatrixAlloc(N_P, D);

    for (i = 0; i < N_P; i++)
      {
        for ( j = 0; j < D; j++)
          {
            r_plus[i][j] = r_minus[i][j] = pos[i][j];
          }
      }
    // compute the kinetic energy
    T_Local = 0;
    for (i = 0; i < N_P; i++)
      {
        for (j = 0; j < D; j++)
          {
            r_plus[i][j]  = pos[i][j] + h_der;
            r_minus[i][j] = pos[i][j] - h_der;
            psi_minus = trialFunction(N_P, D, r_minus, alpha, beta, omega);
            psi_plus  = trialFunction(N_P, D, r_plus,  alpha, beta, omega);
            T_Local = T_Local - (psi_minus + psi_plus - 2*psi_curr);
            r_plus[i][j]  = pos[i][j];
            r_minus[i][j] = pos[i][j];
          }
      }
    T_Local = 0.5*h2_der*T_Local/(psi_curr);
    MatrixDeAlloc(r_plus,  N_P); // free memory
    MatrixDeAlloc(r_minus, N_P);
    return T_Local;
  }

double TrialWaveFunction1(int N_P, int D, double** &pos, double alpha, double beta, double omega)
  {
    double argument, wavefunction, single_particle;
    argument = wavefunction = 0;

    for(int i = 0; i < N_P; i++)
      {
        single_particle = 0;
        for(int j = 0; j < D; j++)
          {
            single_particle += pos[i][j]*pos[i][j];
          }
        argument += single_particle;
      }
    wavefunction = exp(-0.5*alpha*omega*argument);
    return wavefunction;
  }

double TrialWaveFunction2(int N_P, int D, double** &pos, double alpha, double beta, double omega)
  {
    double r_12 = sqrt( r_12_squared(N_P, D, pos) );
    double factor = 1.0 + beta*r_12;
    double Jastrow = exp( r_12/(2*factor) );
    double wf1 = TrialWaveFunction1(N_P, D, pos, alpha, beta, omega);
    return wf1*Jastrow;
  }


double r_12_squared(int N_P, int D, double** &pos)
// Distance between electrons squared
  {
    double value = 0;
    for(int i = 0; i < D; i++)
      {
        value += ( pos[0][i] - pos[1][i] ) * ( pos[0][i] - pos[1][i] );
      }
    return value;
  }

double sum_r1r2_squared(int N_P, int D, double** &pos)
// sum of the norms of positions^2
  {
    double r_sq = 0;
    for(int i = 0; i < N_P; i++)
      {
        r_sq += (pos[i][0]*pos[i][0] + pos[i][1]*pos[i][1] + pos[i][2]*pos[i][2]);
      }
    return r_sq;
  }
double E_Local_Analytic_T1U(int N_P, int D, double** &pos, double alpha, double beta, double omega)
  {
    return 0.5*omega*omega*sum_r1r2_squared(N_P, D, pos)*(1 - alpha*alpha) + 3*alpha*omega;
  }
double E_Local_Analytic_T1P(int N_P, int D, double** &pos, double alpha, double beta, double omega)
  {
    double r_12 = sqrt(r_12_squared(N_P, D, pos));
    return E_Local_Analytic_T1U(N_P, D, pos, alpha, beta, omega) + 1.0/r_12;
  }

double E_Local_Analytic_T2P(int N_P, int D, double** &pos, double alpha, double beta, double omega)
  {
    double r_12 = sqrt( r_12_squared(N_P, D, pos) );
    double arg = 1 + beta*r_12;
    double arg2 = 1/(2*arg*arg);
    double term = arg2 * ( alpha*omega*r_12 - arg2 -2/r_12 + (2*beta)/arg );
    double ET1P = E_Local_Analytic_T1P(N_P, D, pos, alpha, beta, omega);
//    cout << ET1P + term;
    return ET1P + term;
  }

double** MatrixAlloc(int N_P, int D)
/*Function which returns dynamically allocates memory for a square nxn matrix*/
  {
    double** Matrix = new double*[N_P];
    for(int i = 0; i < N_P; i++)
      {
        Matrix[i] = new double[D];
      }
    return Matrix;
  }
void MatrixDeAlloc(double** &Matrix, int N_P)
/*Void for deallocating dynamic matrix*/
  {
    for(int i = 0; i < N_P; i++) delete[] Matrix[i];
    delete[] Matrix;
  }

void PrintMatrix(double** &Matr, int N_P, int D)
  {
      for(int i = 0; i < N_P ; i++)
        {
          for(int j = 0; j < D; j++) cout << setw(7) << Matr[i][j] << " ";
          cout << endl;
        }
    }
