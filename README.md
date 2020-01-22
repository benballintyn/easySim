# easySim
Simulator for network simulations of exponential-LIF neurons that utilizes MATLAB's gpuCoder

# Description
easySim provides an easy to use Object-oriented format for creating networks of neurons. It then uses MATLAB's GPU Coder
to compile optimized CUDA code to actually run the network. While more benchmarking is necessary, currently the simulator runs at ~.33-.5x real time meaning simulating 10s of a network takes ~20-30s (this was done for a network with 4000 cells). Currently "EVLIF" neurons (exponential leaky integrate and fire with variable threshold) and "AEVLIF" (adaptive exponential leaky integrate and fire with variable threshold) are supported although others (LIF, ELIF) will gain support in the near future.

# Installation
easySim requires a number of packages in order to utilize the GPU Coder functionality. I recommend updating MATLAB to the 2019b version and installing the 10.1 version of CUDA toolkit (this is compatible with the 2019b MATLAB). For additional installation instructions for GPU Coder please see [GPU Coder docs](https://www.mathworks.com/help/pdf_doc/gpucoder/gpucoder_gs.pdf).
