# MOPSODTC-Algorithm

This repository contains the MATLAB source code for the algorithm **A Multi-Objective Particle Swarm Optimization Algorithm with Dynamic Two-Population Co-evolution (MOPSODTC)**.

This code was developed to accompany the paper submission and to ensure the reproducibility of our research results.

## Algorithm Overview

MOPSODTC is a multi-objective particle swarm optimization (MOPSO) algorithm designed to improve the balance between convergence and diversity. Its main features include:
- A **dynamic two-population co-evolutionary framework** that divides the swarm into an Exploration Swarm (ES) and a Development Swarm (DS).
- **Tailored update strategies** for each sub-population to focus on diversity (ES) and convergence (DS) respectively.
- The use of **evolutionary operators** (Simulated Binary Crossover and Polynomial Mutation) to enhance particle quality.
- An **adaptive external archive maintenance strategy** using a modified Shift-based Mean (SM) density estimation to preserve high-quality, well-distributed solutions.

## Requirements

- **MATLAB** (tested on R2020b and newer).
- The true Pareto Front data files (`ZDT1_PF.mat`, `ZDT3_PF.mat`, etc.) are required for plotting comparisons. These can be obtained from academic sources like the PlatEMO platform.

## How to Run the Code

1.  **Clone or Download:** Download all the `.m` files from this repository into a single folder on your computer.
2.  **Open MATLAB:** Open MATLAB and navigate to the folder where you saved the files.
3.  **Run the Main Script:** Execute the `run_MOPSODTC.m` script from the MATLAB command window or editor.

    ```matlab
    run_MOPSODTC
    ```

4.  **Configuration:** You can easily change the target problem and algorithm parameters directly within the `run_MOPSODTC.m` file. For example, to switch from ZDT3 to ZDT1, change this line:
    ```matlab
    problem_name = 'ZDT3'; % Change to 'ZDT1', 'ZDT2', etc.
    ```

## File Structure

- `run_MOPSODTC.m`: The main script to execute an experiment.
- `MOPSODTC.m`: The core implementation of the MOPSODTC algorithm.
- `TestProblems.m`: Definitions for the benchmark test functions (ZDT series).
- `UpdateArchive.m`: Implements the external archive management, including the pruning strategy.
- `Dominates.m`: A helper function to determine Pareto dominance.
- `SBX.m`: Implementation of the Simulated Binary Crossover operator.
- `PM.m`: Implementation of the Polynomial Mutation operator.

## Citation

If you use this code in your research, please cite our paper:

```
[Your Paper's Citation Information will go here once it is published]
```
