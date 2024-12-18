#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include <thread>
#include"Models.hpp"
#include<iostream>
#include<cmath>
#include<algorithm>
#include<vector>
#include<ctime>
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include<fstream>  
#include<string>
#include <numeric>
using namespace std;

// Normal generator of a single normal random number
double GenNormal(boost::variate_generator<boost::mt19937, boost::normal_distribution<> >& generator);

// Function to do Cholesky decomposition
vector< vector<double> > Cholesky(vector< vector<double> >& data);

double CaclVariance(vector<double>& data);

// Monte Carlo pricing method
double MonteCarlo(int NumPath, int NumIncre, Option& Opt, Heston& Hest, double* ThreadSum = nullptr, int thread = -1);

// Monte Carlo pricing method with multithreading
double MonteCarloMultiThread(int NumThreads, int NumPath, int NumIncre, Option& Opt, Heston& Hest, bool Antithetic = false);

// Monte Carlo pricing method with Antithetic Variates
double MonteCarloAntitheticVariates(int NumPath, int NumIncre, Option& Opt, Heston& Hest, double* ThreadSum = nullptr, int thread = -1);

#endif 