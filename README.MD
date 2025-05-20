# Applied Macroeconomic Modelling – A Crash Course on Climate Economics

This repository accompanies the lecture **"A Crash Course on Climate Economics"**, covering both long-run integrated assessment models (IAMs) like DICE and short-run stochastic dynamic general equilibrium (DSGE) climate models.

It includes code, handouts, and slides for replicating and extending the results presented in the course.

---

## ⚠️ Requirements

> **The code in this repository runs exclusively with:**
>
> - **MATLAB**
> - **Dynare version 6.xx**
>
> It is not compatible with Octave or earlier versions of Dynare.

---

## Repository Structure

### `/A_functions/`

Utility functions to solve nonlinear perfect foresight problems using **extended path methods**.

---

### `/B_RBC/` – Baby Stochastic IAM (Real Business Cycle version)

This folder contains a stylized dynamic general equilibrium model for climate-economy interaction, inspired by DICE but adapted to the RBC framework.

- `get_Z.m`: Computes long-run emissions and total factor productivity (TFP) trends.
- `guess_path.m`: Generates initial guess paths for the solver (not mandatory but improves convergence speed).
- `model_file.mod`: The core Dynare file for the RBC-based IAM model.
- `run1_simple_bau.mod`: Simulates the **Business-as-Usual** deterministic path, similar in spirit to DICE.
- `run2_simple_bau.mod`: Compares deterministic paths under **Paris Agreement** vs **BAU**.
- `run3_compare_nz0_bau_stochastic.mod`: Simulates the **stochastic IAM** under productivity shocks.

---

### `/C_NKC/` – New Keynesian Climate Model

This directory contains the New Keynesian Climate (NKC) model developed by **Sahuc, Smets, and Vermandel (2025)**.

- `model_file.mod`: Dynare code for the NKC model.
- `run_estimation_filtered_data.mat`: MATLAB file with the filtered data used in model estimation.
- `run_simulations.mod`: Runs **transition scenario simulations** (e.g., climate policies, shocks).

---

## Teaching Materials

- `Slides.pdf`: Presentation slides used in the lecture.
  - Covers both the **DICE** model (long-run, deterministic IAM) and **E-DSGE** models (stochastic, policy-oriented models).
  - Includes methodological insights on solving dynamic models far from steady state.

- `Handout_BabyIAM.pdf`: A step-by-step handout to guide you through the simulation and understanding of the Baby IAM model (RBC version).

---

## Citation

If you use or adapt this material, please cite:

> Sahuc, J.-G., Smets, F., & Vermandel, G. (2025). *The New Keynesian Climate Model*.

---

For questions or contributions, contact [gauthier@vermandel.fr](mailto:gauthier@vermandel.fr)
