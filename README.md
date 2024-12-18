# Stochastic Volatility Option Pricing Using Monte Carlo Simulations

This project uses Monte Carlo simulations applied to the Heston Stochastic Volatility model to price both vanilla European options, as well as Asian and Lookback exotic options. A performance comparison is made between the following methods:

1. Ordinary Monte Carlo Simulation
2. Antithetic Variate Monte Carlo Simulation
3. Multi-Threaded Ordinary Monte Carlo Simulation
4. Multi-Threaded Antithetic Variate Monte Carlo Simulation

## Heston Stochastic Volatility Model

The Heston model is a well-known stochastic volatility model that allows for the modeling of option prices by assuming that the volatility of the underlying asset follows a mean-reverting square-root process. This model is useful because it can capture the observed phenomena of volatility clustering and leverage effects, unlike the simpler constant volatility models (e.g., Black-Scholes).

The Heston model is described by the following system of stochastic differential equations (SDEs):

$$
dS_t = \mu S_t dt + \sqrt{v_t} S_t dW_t^{(1)}
$$

$$
dv_t = \kappa (\theta - v_t) dt + \sigma \sqrt{v_t} dW_t^{(2)}
$$

where:

- $S_t$ is the underlying asset price at time $t$,
- $v_t$ is the variance (volatility squared) at time $t$,
- $\mu$ is the rate of return of the asset,
- $\kappa$ is the rate at which the volatility reverts to the long-term mean $\theta$,
- $\sigma$ is the volatility of volatility (or how volatile the variance is),
- $W_t^{(1)}$ and $W_t^{(2)}$ are two correlated Brownian motions, with correlation $\rho$.

These SDEs govern the evolution of the asset price and its volatility over time. A key characteristic of this model is the correlation between the two Wiener processes, allowing for a realistic capture of volatility dynamics.

## Discretised Version of the Heston Model and Monte Carlo Simulation

To use the Heston model for option pricing, we need to discretise the SDEs in a way that makes them amenable to numerical simulation. The discretisation can be done using the Euler-Maruyama method, which approximates the continuous SDEs using time steps of size $\Delta t$. The discretised versions of the equations are as follows:

For the asset price $S_t$:
$$
S_{t+1} = S_t \cdot \exp\left( \left( \mu - \frac{1}{2} v_t \right) \Delta t + \sqrt{v_t} \cdot \epsilon_1 \sqrt{\Delta t} \right)
$$

For the variance process $v_t$:
$$
v_{t+1} = \max\left( v_t + \kappa (\theta - v_t) \Delta t + \sigma \sqrt{v_t} \cdot \epsilon_2 \sqrt{\Delta t}, 0 \right)
$$

where $\epsilon_1$ and $\epsilon_2$ are standard normal random variables that are correlated with correlation $\rho$, and $\max(\cdot, 0)$ ensures that the variance $v_t$ does not become negative.

## Monte Carlo Simulation

Monte Carlo simulation is used to estimate the price of options under the Heston model by simulating multiple paths of the underlying asset and its volatility. The process begins by generating random numbers from a standard normal distribution, which are then transformed into correlated Brownian motions using Cholesky decomposition.

Cholesky decomposition factors a covariance matrix into a lower triangular matrix, enabling the generation of correlated random variables from independent standard normal variables. These correlated Brownian motions are used to simulate the asset price and variance over time.

The final option price is obtained by averaging the discounted payoffs from all simulated paths, providing an estimate of the option’s price under the Heston model.

## Performance Comparison

Below are the performance results, calculated over 100 trials with 1,000 simulations each. The multi-threaded antithetic variate method outperforms the others, exhibiting both the smallest standard deviation and the shortest runtime.

```
Ordinary MC simulation has an average price of: 13.2713 and standard deviation of: 0.760764
Ordinary MC Multi-Thread simulation has an average price of: 13.6665 and standard deviation of: 0.666989
Antithetic MC simulation has an average price of: 13.4394 and standard deviation of: 0.662268
Antithetic MC Multi Thread simulation has an average price of: 13.8337 and standard deviation of: 0.646228
Total time for Ordinary MC simulations: 448.408 s
Total time for Ordinary MC Multi-Thread simulations: 193.18 s
Total time for Antithetic MC simulations: 294.811 s
Total time for Antithetic MC Multi Thread simulations: 126.749 s
```

## Dependencies

- C++17 or later
- [Boost](https://www.boost.org/) (for random number generation and statistical functions)



