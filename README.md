# Koopman Operator Analysis of Relaxation Dynamics in the 2D Ising Model

This repository contains the code and results from an undergraduate research project investigating **relaxation dynamics in the two-dimensional Ising model** using **Koopman operator theory**.

The central objective of this project is to examine whether the spectral properties of a finite-dimensional approximation to the Koopman operator capture signatures of **critical slowing down** near the ferromagnetic phase transition.

---

## Background

The two-dimensional Ising model exhibits a second-order phase transition at the critical temperature

\[
T_c = \frac{2}{\ln(1+\sqrt{2})} \approx 2.269.
\]

Near this temperature, the relaxation time of the system diverges, leading to the phenomenon of **critical slowing down**.

Instead of analysing relaxation solely through conventional autocorrelation functions, this project adopts an **operator-theoretic viewpoint**, approximating the Koopman operator using **Extended Dynamic Mode Decomposition (EDMD)**.

---

## Methodology

### Simulation

Snapshots were generated from the dynamics of a

- 16 × 16 square Ising lattice
- Periodic boundary conditions
- Glauber dynamics
- Monte Carlo evolution

For each temperature, trajectories were collected after equilibration.

---

### Observables

The EDMD basis consists of two physically meaningful observables:

- **Nearest-neighbour spin correlation**
- **Next-nearest-neighbour spin correlation**

These observables form the feature space from which the finite-dimensional Koopman operator is approximated.

---

### Koopman Approximation

Using the collected trajectory data,

- EDMD constructs a finite-dimensional approximation of the Koopman operator,
- the eigenvalue spectrum is computed,
- and the **spectral gap**

\[
\Delta = |\lambda_1| - |\lambda_2|
\]

is evaluated as a function of temperature.

The spectral gap serves as an indicator of the rate at which observables relax toward equilibrium.

---

## Results

The primary result of this work is the temperature dependence of the Koopman spectral gap.

As the system approaches the Onsager critical temperature,

- the spectral gap decreases significantly,
- indicating increasingly slow relaxation,
- consistent with the phenomenon of **critical slowing down**.

The repository includes the generated spectral-gap plot.

---

## Repository Structure
