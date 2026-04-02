# Monte Carlo Methods for Statistical Physics and Dynamical Systems

This repository contains implementations of Monte Carlo methods applied to the Ising model, with a focus on both **equilibrium thermodynamics** and **dynamical behaviour of stochastic systems**.

---

## 🔹 Project Overview

I study the 2D Ising model using a combination of:

- **Wang–Landau sampling** → to estimate the density of states \( g(E) \)  
- **Glauber dynamics** → to generate stochastic trajectories  
- **Koopman operator approximation (EDMD)** → to analyse system dynamics  

This allows me to connect:
- equilibrium properties (e.g. specific heat)  
- with dynamical properties (e.g. relaxation time, spectral gap)  

---

## 🔹 Key Results

- Estimated density of states \( g(E) \) using Wang–Landau sampling  
- Recovered thermodynamic quantities such as:
  - partition function \( Z(T) \)  
  - internal energy  
  - specific heat \( C_v(T) \)  

- Approximated the **Koopman operator** from simulation data using EDMD  
- Observed a **reduction in spectral gap near the critical temperature**, indicating **critical slowing down**  
- Analysed how dynamical behaviour changes with temperature (parametric dependence)  

---

## 🔹 Methods

### Wang–Landau Algorithm
- Estimates the density of states \( g(E) \)  
- Uses adaptive updates to achieve a flat histogram
- Non-Markovian

### Glauber Dynamics
- Generates a Markov chain at fixed temperature  
- Used to generate 'snapshots' of data for EDMD of Koopman Operator at each temperature

### Koopman Operator (EDMD)
- Data-driven approximation of the evolution operator  
- Extracts spectral information from trajectories  
- Used to estimate relaxation times via spectral gap  

---

## 🔹 Code Structure

- `ising_wl_1d.py`  
  → Wang–Landau simulation (1D)

- `ising_wl_2d.py`  
  → Wang–Landau simulation (2D)

- `ising_koop.py`  
  → Glauber dynamics + Koopman operator (EDMD)

---

## 🔹 Outputs

- Log density of states \( \log g(E) \)  
- Specific heat curves showing phase transition  
- Koopman spectral gap vs temperature  
- Relaxation time vs temperature
  
