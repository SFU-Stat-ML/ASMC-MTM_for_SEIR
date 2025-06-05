# Annealed Sequential Monte Carlo with Adaptive Multiple-Try Metropolis Kernel and Applications to Disease Transmission Models

**Bayesian inference for infectious disease transmission models using advanced Monte Carlo algorithms**

This repository explores efficient Bayesian inference techniques for compartmental epidemic models (SIR/SEIR). We propose enhanced sampling strategies to overcome the inefficiency of traditional MCMC and SMC methods, especially for complex models and high-dimensional data.

## Background

Mathematical and statistical epidemiology uses compartmental models, often expressed as systems of ordinary differential equations (ODEs), to simulate the spread of diseases like COVID-19. Bayesian inference is a popular approach for estimating unknown model parameters using observed data.

While **Markov Chain Monte Carlo (MCMC)** and **Sequential Monte Carlo (SMC)** methods are standard, they can be computationally intensive. To address this, we:

- Integrated a **Multiple-Try Metropolis (MTM)** kernel into the **Annealed SMC (ASMC)** framework (implemented in R and Julia).
- Are applying the **Compound Auxiliary Metropolis (CAM)** method by Doig, R. (2025) in **Julia**, a more advanced and efficient sampling algorithm.


## Methods Implemented

| Model     | Method              | Language | Status      |
|-----------|---------------------|----------|-------------|
| SIR/SEIR  | Adaptive ASMC       | R/ Julia | ✅ Complete |
| SIR/SEIR  | ASMC + MTM          | R/ Julia | ✅ Complete |
| SEIR      | CAM                 | Julia    | ✅ Complete |
