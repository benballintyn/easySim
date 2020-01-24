# easySim
Simulator for network simulations of LIF neuron variants that utilizes MATLAB's gpuCoder

# Description
easySim provides an easy to use Object-oriented format for creating networks of neurons. It then uses MATLAB's GPU Coder
to compile optimized CUDA code to actually run the network. While more benchmarking is necessary, currently the simulator runs at ~.33-.5x real time meaning simulating 10s of a network takes ~20-30s (this was done for a network with 4000 cells). Currently "EVLIF" neurons (exponential leaky integrate and fire with variable threshold) and "AEVLIF" (adaptive exponential leaky integrate and fire with variable threshold) are supported although others (LIF, ELIF) will gain support in the near future.

# Installation
easySim requires a number of packages in order to utilize the GPU Coder functionality. I recommend updating MATLAB to the 2019b version and installing the 10.1 version of CUDA toolkit (this is compatible with the 2019b MATLAB). For additional installation instructions for GPU Coder please see [GPU Coder docs](https://www.mathworks.com/help/pdf_doc/gpucoder/gpucoder_gs.pdf).

For CPU only use, neither CUDA nor GPU Coder are required. Only MATLAB's embedded coder is required.

# Getting Started
To get started, check out the examples in the 'test' folder. There are examples utilizing either the CPU or GPU as well as the currently supported EVLIF or AEVLIF networks.

# Documentation
There is currently no comprehensive documentation page although I intend to make one in the near future. However, there is extensive documentation within the code itself. From the MATLAB command line simply type help 'class/function name here' for in depth information about each class and function (e.g. help EVLIFnetwork).

# Planned additions
At the moment, easySim does not support short or long-term changes in synaptic strength. The additions of synaptic facilitation and depression as well as STDP based learning rules will be added in the coming weeks/months.

Networks supporting simple LIF or Exponential-LIf neurons will be added soon.

If you have any suggestions for additions, please create an issue topic and I will do my best to add it.
