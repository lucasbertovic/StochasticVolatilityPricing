#include"Models.hpp"
#include"Simulation.hpp"
#include <chrono>
#include <utility> 

int main()
{
	// Timer
	clock_t start, end; // typedef long clock_t
	start = clock();
	// Input the number of paths and number of increments
	int NumPath = 1000, NumIncre = 10000;

	// Correlation of underlying and volatility
	double rho = -0.6;
	// Option type
	OptionType OptType = OptionType::EuroCall;

	// Inputs for the option pricing simulation
	double S_0 = 50.0;    // Initial spot price
	double K = 50.0;      // Strike price
	double r = 0.03;     // Risk-free rate
	double v_0 = 0.30; // Initial volatility 
	double T = 1.00;       // One year until expiry

	double Kappa = 1.15;   // Mean-reversion rate
	double Theta = 0.20;  // Long run average volatility
	double Epslon = 0.39;      // Vol of vol

	//Number of threads for multithreading
	int NumThreads = 5;

	Option Opt(K, r, T, S_0, v_0, OptType);
	Heston Hest(Kappa, Theta, Epslon, rho, Opt);

	// Test for effciency of ordinary vs antithetic
	int TrailNum = 100;
	vector<pair<double, double>> OptionPricesMC;    // Store (price, time) for MC
	vector<pair<double, double>> OptionPricesMCMT; // Store (price, time) for MC Multi Thread
	vector<pair<double, double>> OptionPricesAnti; // Store (price, time) for Antithetic MC
	vector<pair<double, double>> OptionPricesAntiMT; // Store (price, time) for Antithetic MC Multi Thread

	for (int i = 0; i < TrailNum; ++i)
	{
		auto startMC = chrono::high_resolution_clock::now(); // Start time for MC
		double TempPriceMC = MonteCarlo(NumPath, NumIncre, Opt, Hest);
		auto endMC = chrono::high_resolution_clock::now();   // End time for MC
		double durationMC = chrono::duration<double>(endMC - startMC).count();

		auto startMCMT = chrono::high_resolution_clock::now(); // Start time for MC Multi-Thread
		double TempPriceMCMT = MonteCarloMultiThread(NumThreads, NumPath, NumIncre, Opt, Hest);
		auto endMCMT = chrono::high_resolution_clock::now();   // End time for MC Multi-Thread
		double durationMCMT = chrono::duration<double>(endMCMT - startMCMT).count();

		auto startAnti = chrono::high_resolution_clock::now(); // Start time for Antithetic MC
		double TempPriceAnti = MonteCarloAntitheticVariates(NumPath, NumIncre, Opt, Hest);
		auto endAnti = chrono::high_resolution_clock::now();   // End time for Antithetic MC
		double durationAnti = chrono::duration<double>(endAnti - startAnti).count();

		auto startAntiMT = chrono::high_resolution_clock::now(); // Start time for Antithetic MC
		double TempPriceAntiMT = MonteCarloMultiThread(NumThreads, NumPath, NumIncre, Opt, Hest, true);
		auto endAntiMT = chrono::high_resolution_clock::now();   // End time for Antithetic MC
		double durationAntiMT = chrono::duration<double>(endAntiMT - startAntiMT).count();

		// Store price and time in respective vectors
		OptionPricesMC.push_back(make_pair(TempPriceMC, durationMC));
		OptionPricesMCMT.push_back(make_pair(TempPriceMCMT, durationMCMT));
		OptionPricesAnti.push_back(make_pair(TempPriceAnti, durationAnti));
		OptionPricesAntiMT.push_back(make_pair(TempPriceAntiMT, durationAntiMT));
	}

	// Separate prices for variance and mean calculations
	vector<double> PricesMC, PricesMCMT, PricesAnti, PricesAntiMT;
	for (const auto& p : OptionPricesMC) PricesMC.push_back(p.first);
	for (const auto& p : OptionPricesMCMT) PricesMCMT.push_back(p.first);
	for (const auto& p : OptionPricesAnti) PricesAnti.push_back(p.first);
	for (const auto& p : OptionPricesAntiMT) PricesAntiMT.push_back(p.first);

	// Calculate variance and mean
	double VarMC = CaclVariance(PricesMC);
	double VarMCMT = CaclVariance(PricesMCMT);
	double VarAnti = CaclVariance(PricesAnti);
	double VarAntiMT = CaclVariance(PricesAntiMT);

	double MeanMC = accumulate(PricesMC.begin(), PricesMC.end(), 0.0) / PricesMC.size();
	double MeanMCMT = accumulate(PricesMCMT.begin(), PricesMCMT.end(), 0.0) / PricesMCMT.size();
	double MeanAnti = accumulate(PricesAnti.begin(), PricesAnti.end(), 0.0) / PricesAnti.size();
	double MeanAntiMT = accumulate(PricesAntiMT.begin(), PricesAntiMT.end(), 0.0) / PricesAntiMT.size();

	// Calculate total time
	double TotalTimeMC = 0, TotalTimeMCMT = 0, TotalTimeAnti = 0, TotalTimeAntiMT = 0;
	for (const auto& p : OptionPricesMC) TotalTimeMC += p.second;
	for (const auto& p : OptionPricesMCMT) TotalTimeMCMT += p.second;
	for (const auto& p : OptionPricesAnti) TotalTimeAnti += p.second;
	for (const auto& p : OptionPricesAntiMT) TotalTimeAntiMT += p.second;

	// Output results
	cout << "Ordinary MC simulation has an average price of: " << MeanMC
		<< " and standard deviation of: " << VarMC << endl;

	cout << "Ordinary MC Multi-Thread simulation has an average price of: " << MeanMCMT
		<< " and standard deviation of: " << VarMCMT << endl;

	cout << "Antithetic MC simulation has an average price of: " << MeanAnti
		<< " and standard deviation of: " << VarAnti << endl;

	cout << "Antithetic MC Multi Thread simulation has an average price of: " << MeanAntiMT
		<< " and standard deviation of: " << VarAntiMT << endl;

	cout << "Total time for Ordinary MC simulations: " << TotalTimeMC << " s" << endl;
	cout << "Total time for Ordinary MC Multi-Thread simulations: " << TotalTimeMCMT << " s" << endl;
	cout << "Total time for Antithetic MC simulations: " << TotalTimeAnti << " s" << endl;
	cout << "Total time for Antithetic MC Multi Thread simulations: " << TotalTimeAntiMT << " s" << endl;

	// Pause
	getchar();
	return 0;
}