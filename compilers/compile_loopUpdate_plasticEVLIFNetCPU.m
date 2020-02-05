function [] = compile_loopUpdate_plasticEVLIFNetCPU()
% This function compiles the loopUpdate_plasticEVLIFNetCPU function into C
% code and produces a mex file loopUpdate_plasticEVLIFNetCPU_mex to run a
% simulation on the CPU.
% compile_loopUpdateEVLIFNetCPU()
maxN = 100000;
maxSpikeGen = 100000;
%% Define argument types for entry-point 'loopUpdateNetGPU'.
ARGS = cell(1,1);
ARGS{1} = cell(42,1);
ARGS{1}{1} = coder.typeof(0,[maxN   1],[true false]); % V
ARGS{1}{2} = coder.typeof(0,[maxN   1],[true false]); % Vreset
ARGS{1}{3} = coder.typeof(0,[maxN   1],[true false]); % tau_ref
ARGS{1}{4} = coder.typeof(0,[maxN   1],[true false]); % Vth
ARGS{1}{5} = coder.typeof(0,[maxN   1],[true false]); % Vth0
ARGS{1}{6} = coder.typeof(0,[maxN   1],[true false]); % Vth_max
ARGS{1}{7} = coder.typeof(0,[maxN   1],[true false]); % VsynE
ARGS{1}{8} = coder.typeof(0,[maxN   1],[true false]); % VsynI
ARGS{1}{9} = coder.typeof(0,[maxN   1],[true false]); % GsynE
ARGS{1}{10} = coder.typeof(0,[maxN   1],[true false]); % GsynI
ARGS{1}{11} = coder.typeof(0,[maxN maxN+maxSpikeGen],[true true]); % GsynMax
ARGS{1}{12} = coder.typeof(0,[maxN+maxSpikeGen 1],[true false]); %p0
ARGS{1}{13} = coder.typeof(0,[maxN   1],[true false]); % tau_synE
ARGS{1}{14} = coder.typeof(0,[maxN   1],[true false]); % tau_synI
ARGS{1}{15} = coder.typeof(0,[maxN   1],[true false]); % Cm
ARGS{1}{16} = coder.typeof(0,[maxN   1],[true false]); % Gl
ARGS{1}{17} = coder.typeof(0,[maxN   1],[true false]); % El
ARGS{1}{18} = coder.typeof(0,[maxN   1],[true false]); % dth
ARGS{1}{19} = coder.typeof(0,[maxN   1],[true false]); % Iapp
ARGS{1}{20} = coder.typeof(0,[maxN  1],[true false]); % std_noise
ARGS{1}{21} = coder.typeof(0); % dt
ARGS{1}{22} = coder.typeof(false,[maxN+maxSpikeGen   1],[true false]); % ecells
ARGS{1}{23} = coder.typeof(false,[maxN+maxSpikeGen   1],[true false]); % icells
ARGS{1}{24} = coder.typeof(0,[maxSpikeGen 1],[true false]); % spikeGenProbs
ARGS{1}{25} = coder.typeof(0,[maxN 1],[true false]); % cells2record
ARGS{1}{26} = coder.typeof(false,[maxN maxN+maxSpikeGen],[true true]); % is_plastic logical matrix
ARGS{1}{27} = coder.typeof('string',[1 20],[false true]); % plasticity_type string
ARGS{1}{28} = coder.typeof(false,[maxN maxN+maxSpikeGen],[true true]); % Connectivity matrix
ARGS{1}{29} = coder.typeof(0,[maxN+maxSpikeGen 1],[true false]); % r1
ARGS{1}{30} = coder.typeof(0,[maxN+maxSpikeGen 1],[true false]); % r2
ARGS{1}{31} = coder.typeof(0,[maxN 1],[true false]); % o1
ARGS{1}{32} = coder.typeof(0,[maxN 1],[true false]); % o2
ARGS{1}{33} = coder.typeof(0,[maxN+maxSpikeGen 1],[true false]); % A2plus
ARGS{1}{34} = coder.typeof(0,[maxN+maxSpikeGen 1],[true false]); % A3plus
ARGS{1}{35} = coder.typeof(0,[maxN 1],[true false]); % A2minus
ARGS{1}{36} = coder.typeof(0,[maxN 1],[true false]); % A3minus
ARGS{1}{37} = coder.typeof(0,[maxN+maxSpikeGen 1],[true false]); % tau_plus
ARGS{1}{38} = coder.typeof(0,[maxN+maxSpikeGen 1],[true false]); % tau_x
ARGS{1}{39} = coder.typeof(0,[maxN 1],[true false]); % tau_minus
ARGS{1}{40} = coder.typeof(0,[maxN 1],[true false]); % tau_y
ARGS{1}{41} = coder.typeof(0); % nT (# of timesteps to simulate)
ARGS{1}{42} = coder.typeof(0); % file ID for spikes

%% Invoke MATLAB Coder.
codegen loopUpdate_plasticEVLIFNetCPU -args ARGS{1}
end

