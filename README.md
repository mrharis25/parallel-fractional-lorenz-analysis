# Parallel Fractional-Order Lorenz Analysis

MATLAB implementation for parallel Lyapunov exponent and bifurcation analysis of a fractional-order Lorenz system.

## Description

This repository provides a parallel MATLAB implementation for computing the Lyapunov exponent spectrum and bifurcation diagram of a fractional-order Lorenz chaotic system.

The parameter scan is parallelised using MATLAB's `parfor` structure. The Lyapunov exponents are calculated using the Gram–Schmidt reorthonormalisation procedure, while the bifurcation diagram is constructed from the local maxima of the state variable \(x\).

## Requirements

- MATLAB
- Parallel Computing Toolbox

## Files

- `Main_Parfor_LEandBif.m`: Main simulation and parameter-scan code
- `memo.m`: Fractional-order memory-term calculation
- `GSR.m`: Gram–Schmidt reorthonormalisation procedure

## Usage

Download all files and place them in the same MATLAB working directory.

Run:


Main_Parfor_LEandBif


## Citation

If you use this code, its algorithms, or a modified version in an academic publication, please cite the following articles:

1. Li, H., Shen, Y., Han, Y., Dong, J., & Li, J. (2023). Determining Lyapunov exponents of fractional-order systems: A general method based on memory principle. Chaos, Solitons & Fractals, 168, 113167.

2. H. Calgan and Hang Li,  
   “Parallel Bifurcation and Lyapunov Analysis of Commensurate and Incommensurate Fractional-Order Lorenz Systems,”  
   Manuscript under review.

## Licence

This project is licensed under the MIT Licence. See the [LICENSE](LICENSE) file for details.
