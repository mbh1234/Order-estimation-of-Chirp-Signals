# README: Estimation of Number of Components for Chirp and Sinusoidal Signals

## Project Overview

This project aims to estimate the number of components in chirp and sinusoidal signals using various statistical estimators under different noise conditions. The implemented estimators include:

- **MAP (Maximum A Posteriori)**
- **AIC (Akaike Information Criterion)**
- **BIC (Bayesian Information Criterion)**
- **Corrected AIC & BIC**
- **PAL (Penalized Likelihood)**
- **GIC (Generalized Information Criterion)**

Noise distributions considered in the simulations include:

- Normal (Gaussian)
- Sum of Gaussians (mixture)
- t-distribution

### Files in the Repository

1. **`order_est_code.m`**: The main code for estimating the number of components in the signal.
2. **`data_nor.m`**: Generates data with noise following a normal distribution.
3. **`data_n_mix.m`**: Generates data with noise following a mixture of normal distributions.
4. **`data_t.m`**: Generates data with noise following a t-distribution.
5. **`obj_L2_func.m`**: Objective function used in component estimation.
6. **`per_maxm.m`**: Code for performing maximum likelihood estimation.

## Getting Started

### Prerequisites

- MATLAB (Tested on R2021b, but should work on most recent versions)
- Basic understanding of MATLAB scripting and statistical estimators

### Directory Structure

Place all the provided files in the same directory to ensure proper execution of the scripts.

## Usage Instructions

### Step 1: Data Generation

Generate datasets with different noise distributions. Use one of the following scripts based on the noise type you want:

- **For Normal noise**:

  ```matlab
  data_nor;
  ```

  This script will generate a dataset with normal noise and save it to the workspace.

- **For Mixture of Normals noise**:

  ```matlab
  data_n_mix;
  ```

  This script will generate a dataset with a mixture of normal distributions.

- **For t-distribution noise**:

  ```matlab
  data_t;
  ```

  This script will generate a dataset with t-distributed noise.

### Step 2: Component Estimation

Run the main script `order_est_code.m` to estimate the number of components in the signal. Ensure the dataset generated in Step 1 is loaded in the workspace.

```matlab
order_est_code;
```

This script will:

- Load the dataset from the workspace
- Apply various statistical estimators
- Output the estimated number of components for each method

### Step 3: Understanding Results

The output of `order_est_code.m` provides the estimated number of components using the following estimators:

- AIC
- BIC
- MAP
- Corrected AIC
- Corrected BIC
- PAL
- GIC

### Additional Functions

- **`obj_L2_func.m`**: Computes the objective function for each estimator. This is called internally by `order_est_code.m`.
- **`per_maxm.m`**: Performs parameter maximization and likelihood calculations for the estimators.

## Example Execution

Here’s an example sequence to estimate components for a signal with normal noise:

1. Generate the dataset:

   ```matlab
   data_nor;
   ```

2. Estimate the components:

   ```matlab
   order_est_code;
   ```

3. View the results in the MATLAB console to interpret the number of components estimated by each statistical method.

## Notes

- Ensure proper random seed initialization for reproducibility.
- Modify the parameters within `data_nor.m`, `data_n_mix.m`, or `data_t.m` to test different noise levels and scenarios.

## Troubleshooting

- **Error: Undefined variable or function**: Ensure all files are in the same directory and properly named.
- **Unexpected results**: Verify the noise generation script parameters and dataset consistency.

## Contact

For further assistance, please contact the project author or refer to the MATLAB documentation for any function-related queries.
