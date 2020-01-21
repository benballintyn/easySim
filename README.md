# easySim
Simulator for network simulations of exponential-LIF neurons that utilizes MATLAB's gpuCoder

# Description
easySim provide an easy to use Object-oriented format for creating networks of neurons. It then uses MATLAB's GPU Coder
to compile optimized CUDA code to actually run the network. While more benchmarking is necessary, currently the simulator runs at ~.33-.5x real time meaning simulating 10s of a network takes ~20-30s (this was done for a network with 4000 cells). 
