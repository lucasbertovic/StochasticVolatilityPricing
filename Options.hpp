#ifndef OPTIONS_HPP
#define OPTIONS_HPP

#include <algorithm>
#include <string>
#include <vector>
#include <numeric>
using namespace std;

enum class OptionType {
	EuroCall,
	EuroPut,
	AsianCall,
	AsianPut,
	LookbackCall,
	LookbackPut
};

class Option
{
public:
	// Input parameters
	double K;
	double r;
	double T;
	double S_0;
	double v_0;
	OptionType type;
	// Constructor and destructor
	Option(double s_K, double s_r, double s_T, double s_S0, double s_v0, OptionType optionType);
	virtual ~Option();

	double payoff(const vector<double>& stockPath) const;

};

#endif 