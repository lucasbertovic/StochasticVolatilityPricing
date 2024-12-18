#include"Options.hpp"
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <stdexcept>
using namespace std;


Option::Option(double s_K, double s_r, double s_T, double s_S0, double s_v0, OptionType s_type)
{
    K = s_K;
    r = s_r;
    T = s_T;
    S_0 = s_S0;
    v_0 = s_v0;
    type = s_type;
}

double Option::payoff(const vector<double>& stockPath) const
{
    if (stockPath.empty()) {
        throw std::invalid_argument("Stock path cannot be empty.");
    }

    double finalPrice = stockPath.back();
    double averagePrice = std::accumulate(stockPath.begin(), stockPath.end(), 0.0) / stockPath.size();

    double minPrice = *std::min_element(stockPath.begin(), stockPath.end());
    double maxPrice = *std::max_element(stockPath.begin(), stockPath.end());

    auto max_payoff = [](double x) { return std::max(x, 0.0); };

    switch (type) {
    case OptionType::EuroCall:
        return max_payoff(finalPrice - K);
    case OptionType::EuroPut:
        return max_payoff(K - finalPrice);
    case OptionType::AsianCall:
        return max_payoff(averagePrice - K);
    case OptionType::AsianPut:
        return max_payoff(K - averagePrice);
    case OptionType::LookbackCall:
        return max_payoff(finalPrice - minPrice);  // Lookback Call payoff formula
    case OptionType::LookbackPut:
        return max_payoff(maxPrice - finalPrice);  // Lookback Put payoff formula
    default:
        throw std::logic_error("Invalid option type.");
    }
}

Option::~Option() {}