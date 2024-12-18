#include"Simulation.hpp"


// Normal generator of a single normal random number
double GenNormal(boost::variate_generator<boost::mt19937, boost::normal_distribution<> >& generator)
{
	return generator();
}

// Function to do Cholesky decomposition
vector< vector<double> > Cholesky(vector< vector<double> >& data)
{
	int n = data.size();
	vector< vector<double> > mat(n, vector<double>(n));
	double sum1 = 0.0;
	double sum2 = 0.0;
	double sum3 = 0.0;
	// Initialize the first element
	mat[0][0] = sqrt(data[0][0]);

	// First elements of each row
	for (int j = 1; j <= n - 1; j++)
	{
		mat[j][0] = data[j][0] / mat[0][0];
	}
	for (int i = 1; i <= (n - 2); i++)
	{
		for (int k = 0; k <= (i - 1); k++)
		{
			sum1 += pow(mat[i][k], 2);
		}
		mat[i][i] = sqrt(data[i][i] - sum1);
		for (int j = (i + 1); j <= (n - 1); j++)
		{
			for (int k = 0; k <= (i - 1); k++)
			{
				sum2 += mat[j][k] * mat[i][k];
			}
			mat[j][i] = (data[j][i] - sum2) / mat[i][i];
		}
	}
	for (int k = 0; k <= (n - 2); k++)
	{
		sum3 += pow(mat[n - 1][k], 2);
	}
	mat[n - 1][n - 1] = sqrt(data[n - 1][n - 1] - sum3);

	return mat;
}

double CaclVariance(vector<double>& data)
{
	double sum = std::accumulate(data.begin(), data.end(), 0.0);
	double mean = sum / data.size();
	double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
	double stdev = std::sqrt(sq_sum / data.size() - mean * mean);
	return stdev;
}

// Monte Carlo pricing method
double MonteCarlo(int NumPath, int NumIncre, Option& Opt, Heston& Hest, double* ThreadSum, int thread)
{

	// Create a thread-specific seed
	unsigned int seed = static_cast<unsigned int>(time(0));
	if (thread >= 0) {
		seed += thread * 10000; // Add a large multiple of thread ID for unique seeds
	}

	// Boost normal generator
	boost::mt19937 rng(seed);
	boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
		generator(rng, boost::normal_distribution<>());


	// Get corr matrix
	vector< vector<double> > Corr(2, vector<double>(2));
	Corr[0][0] = 1.0, Corr[0][1] = Hest.rho, Corr[1][0] = Hest.rho, Corr[1][1] = 1.0;
	// Cholesky decomposition of corr matrix
	vector< vector<double> > CorrMatrix = Cholesky(Corr);


	// Generate two uncorrelated vectors
	vector< vector<double> > rnd(NumIncre, vector<double>(2));
	// Generate two correlated vectors
	vector<double> DW1(NumIncre), DW2(NumIncre);


	// Monte Carlo simulation
	double SumPayOff = 0.0;
	// The final asset price
	double ST;
	// Create the spot and vol paths
	vector<double> S_t(NumIncre, Opt.S_0);   // Vector of initial spot prices
	vector<double> v_t(NumIncre, Opt.v_0);   // Vector of initial vol prices

	// Create each path every time
	for (int s = 0; s < NumPath; ++s)
	{
		// Create random numbers
		for (int t = 0; t < NumIncre; ++t)
		{
			rnd[t][0] = GenNormal(generator);
			rnd[t][1] = GenNormal(generator);
			DW1[t] = rnd[t][0] * CorrMatrix[0][0] + rnd[t][1] * CorrMatrix[1][0];
			DW2[t] = rnd[t][0] * CorrMatrix[0][1] + rnd[t][1] * CorrMatrix[1][1];
		}
		Hest.GetVolPath(DW2, v_t);
		Hest.GetUnderlyingPath(DW1, v_t, S_t);

		ST = Opt.payoff(S_t);

		SumPayOff += ST;

	}
	// Discounted payoff
	double OptionPrice = SumPayOff / (NumPath * exp(Opt.r * Opt.T));

	if (ThreadSum != nullptr)
	{
		*ThreadSum = SumPayOff;
	}

	return OptionPrice;
}

// Monte Carlo pricing method with multithreading
double MonteCarloMultiThread(int NumThreads, int NumPath, int NumIncre, Option& Opt, Heston& Hest, bool Antithetic)
{
	vector<thread> Threads(NumThreads);
	vector<double> SumValues(NumThreads);

	for (int i = 0; i < NumThreads; ++i)
	{
		int PerPaths = NumPath / NumThreads;
		if (i == 0) {
			PerPaths += NumPath % NumThreads; // Thread 0 gets remaining paths
		}

		// Clone Option and Heston for each thread
		Option LocalOpt = Opt;
		Heston LocalHest = Hest;

		Threads[i] = std::thread([&, i, PerPaths, LocalOpt, LocalHest]() mutable {
			if (!Antithetic) {
				MonteCarlo(PerPaths, NumIncre, LocalOpt, LocalHest, &SumValues[i], i);
			}
			else {
				MonteCarloAntitheticVariates(PerPaths, NumIncre, LocalOpt, LocalHest, &SumValues[i], i);
			}
			});
	}

	// Join all threads
	for (auto& thread : Threads) {
		if (thread.joinable()) {
			thread.join();
		}
	}

	// Calculate option price (sum all thread results)
	double TotalSum = std::accumulate(SumValues.begin(), SumValues.end(), 0.0);
	double OptionPrice = TotalSum / NumPath;
	return OptionPrice;
}

// Monte Carlo pricing method with Antithetic Variates
double MonteCarloAntitheticVariates(int NumPath, int NumIncre, Option& Opt, Heston& Hest, double* ThreadSum, int thread)
{

	// Create a thread-specific seed
	unsigned int seed = static_cast<unsigned int>(time(0));
	if (thread >= 0) {
		seed += thread * 10000; // Add a large multiple of thread ID for unique seeds
	}

	// Boost normal generator
	boost::mt19937 rng(seed);
	boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
		generator(rng, boost::normal_distribution<>());

	// Get corr matrix
	vector< vector<double> > Corr(2, vector<double>(2));
	Corr[0][0] = 1.0, Corr[0][1] = Hest.rho, Corr[1][0] = Hest.rho, Corr[1][1] = 1.0;
	// Cholesky decomposition of corr matrix
	vector< vector<double> > CorrMatrix = Cholesky(Corr);


	// Generate two uncorrelated vectors
	vector< vector<double> > rnd(NumIncre, vector<double>(2));
	// Generate two correlated vectors
	vector<double> DW1(NumIncre), DW2(NumIncre);
	vector<double> DW3(NumIncre), DW4(NumIncre);

	// Monte Carlo simulation
	double SumPayOff = 0.0;
	// The final asset price
	double ST1 = 0.0, ST2 = 0.0;
	// Create the spot and vol paths
	vector<double> S_t1(NumIncre, Opt.S_0);   // Vector of initial spot prices
	vector<double> v_t1(NumIncre, Opt.v_0);   // Vector of initial vol prices
	vector<double> S_t2(NumIncre, Opt.S_0);   // Vector of initial spot prices
	vector<double> v_t2(NumIncre, Opt.v_0);   // Vector of initial vol prices

	// Create each path every time
	for (int s = 0; s < NumPath / 2; ++s)
	{
		// Create random numbers
		for (int t = 0; t < NumIncre; ++t)
		{
			rnd[t][0] = GenNormal(generator);
			rnd[t][1] = GenNormal(generator);
			DW1[t] = rnd[t][0] * CorrMatrix[0][0] + rnd[t][1] * CorrMatrix[1][0];
			DW2[t] = rnd[t][0] * CorrMatrix[0][1] + rnd[t][1] * CorrMatrix[1][1];
			DW3[t] = -rnd[t][0] * CorrMatrix[0][0] - rnd[t][1] * CorrMatrix[1][0];
			DW4[t] = -rnd[t][0] * CorrMatrix[0][1] - rnd[t][1] * CorrMatrix[1][1];
		}
		Hest.GetVolPath(DW2, v_t1);
		Hest.GetUnderlyingPath(DW1, v_t1, S_t1);
		Hest.GetVolPath(DW4, v_t2);
		Hest.GetUnderlyingPath(DW3, v_t2, S_t2);
		// Calculate the final payoff of the call option
		ST1 = Opt.payoff(S_t1);
		ST2 = Opt.payoff(S_t2);
		SumPayOff += ST1 + ST2;
	}
	// Discounted payoff
	double OptionPrice = SumPayOff / (NumPath * exp(Opt.r * Opt.T));

	if (ThreadSum != nullptr)
	{
		*ThreadSum = SumPayOff;
	}

	return OptionPrice;
}