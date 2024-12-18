#ifndef MODELS_HPP
#define MODELS_HPP

#include<algorithm>
#include<iostream>
#include<vector>
#include"Options.hpp"
#include<cmath>
using namespace std;

class Heston
{

public:
	double Kappa;
	double Theta;
	double Epslon;
	double rho;
	const Option& Opt;

	// Constructor and destructor
	Heston(const double& kap, const double& ta, const double& eps, const double& ro, const Option& op);
	virtual ~Heston();

	// Calculate the volatility path
	void GetVolPath(const vector<double>& DW2, vector<double>& v_t);
	// Calculate the asset price path
	void GetUnderlyingPath(const vector<double>& DW1, const vector<double>& v_t,
		vector<double>& S_t) const;

};

#endif 