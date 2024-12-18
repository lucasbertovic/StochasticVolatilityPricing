#include"Models.hpp"

// Consturctor and destructor
Heston::Heston(const double& kap, const double& ta, const double& eps, const double& ro, const Option& op)
    : Kappa(kap), Theta(ta), Epslon(eps), rho(ro), Opt(op) {
    // No copying occurs; Opt refers to the same Option object passed in.
}

Heston::~Heston() {}

// Function to get the volatility path 
void Heston::GetVolPath(const vector<double>& DW2, vector<double>& v_t)
{
    int numSteps = DW2.size();
    double dt = Opt.T / static_cast<double>(numSteps);

    for (int i = 1; i < numSteps; ++i)
    {
        // Ensure non-negative variance before computing the square root
        double clampedVol = std::max(v_t[i - 1], 0.0);
        double sqrtTerm = sqrt(std::max(clampedVol, 1e-8));

        // Euler-Maruyama update
        double nextVol = v_t[i - 1] + Kappa * dt * (Theta - clampedVol)
            + Epslon * sqrtTerm * sqrt(dt) * DW2[i - 1];

        // Reflecting boundary: ensure non-negative variance
        v_t[i] = std::max(nextVol, 0.0);
    }
}


// Calculate the asset price path
void Heston::GetUnderlyingPath(const vector<double>& DW1, const vector<double>& v_t,
    vector<double>& S_t) const
{
    int numSteps = DW1.size();
    double dt = Opt.T / static_cast<double>(numSteps); // Time increment
    double sqrt_dt = sqrt(dt);                         // Precompute sqrt(dt)

    for (int i = 1; i < numSteps; ++i)
    {
        // Ensure non-negative volatility for numerical stability
        double clampedVol = std::max(v_t[i - 1], 1e-8);
        double sqrtVol = sqrt(clampedVol) * sqrt_dt;

        // Update the underlying asset price
        S_t[i] = S_t[i - 1] * exp((Opt.r - 0.5 * clampedVol) * dt + sqrtVol * DW1[i - 1]);
    }
}
