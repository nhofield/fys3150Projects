#ifndef P5LIB_ADD_H
#define P5LIB_ADD_H
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
using namespace std;

#define h_der 0.0001     //Step-size for the differentiation
#define h2_der 100000000  //Inverse squared Step-size for the differentiation

/*--------------------------
|   File-output function   |
/*------------------------*/
void output(int N_P, int D, int MCcycles, double* &StoreValues, double alpha, double beta, double omega);

/*--------------------------------
|    Main-Algorithm functions    |
/*------------------------------*/
void MetropolisSampling(int N_P, int D, int MCcycles, double* &StoreValues, double alpha, double beta, double omega, double h_R,
                        double (*trialFunction)(int, int, double** &, double, double, double),
                        double* (*localEnergy) (int, int, double** &, double, double, double, double, double,
                        double (*trialFunction)(int, int, double** &, double, double, double)),
                        double (*localEnergyAnalytic)(int, int, double** &, double, double, double) );

double LocateOptimal_h_R(int N_P, int D, double alpha, double beta, double omega,
                        double (*trialFunction)(int, int, double** &, double, double, double) );

/*-------------------------
|   Trial wavefunctions   |
/*-----------------------*/
double TrialWaveFunction1(int N_P, int D, double** &pos, double alpha, double beta, double omega);
double TrialWaveFunction2(int N_P, int D, double** &pos, double alpha, double beta, double omega);

/*---------------------------------------------
|   Energy-specific computational functions   |
/*-------------------------------------------*/
double Kinetic_Local(int N_P, int D, double** &pos, double psi_curr, double h_R, double alpha, double beta, double omega,
          double (*trialFunction)(int, int, double** &, double, double, double) );
double U_Harmonic_osc(int N_P, int D, double omega, double** &pos);
double U_ee_Repulsion(int N_P, int D, double** &pos);

/*----------------------------------------
|   Analytical local energy functions    |
/*--------------------------------------*/
double E_Local_Analytic_T1U(int N_P, int D, double** &pos, double alpha, double beta, double omega);
double E_Local_Analytic_T1P(int N_P, int D, double** &pos, double alpha, double beta, double omega);
double E_Local_Analytic_T2P(int N_P, int D, double** &pos, double alpha, double beta, double omega);

/*-------------------------------------------------
|   Positional-specific computational functions   |
/*-----------------------------------------------*/
double r_12_squared(int N_P, int D, double** &pos);      // returns distance between two particles |vec{r_1}-vec{r_2}|
double sum_r1r2_squared(int N_P, int D, double** &pos);  // returns cartesian sum of squared positions r_1^2 + r_2^2

/*-------------------------------
|   Matrix specific functions   |
/*-----------------------------*/
double** MatrixAlloc(int, int);
void MatrixDeAlloc(double** &, int);
void PrintMatrix(double** &, int, int);

#endif
