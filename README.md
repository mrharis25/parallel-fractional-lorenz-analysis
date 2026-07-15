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
