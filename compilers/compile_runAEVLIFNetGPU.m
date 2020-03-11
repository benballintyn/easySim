function [] = compile_runAEVLIFNetGPU(net,nCells2record)
% This function compiles the loopUpdateEVLIFNetGPU_fast function into CUDA
% code and produces a mex file loopUpdateEVLIFNetGPU_fast_mex to run a
% simulation on the GPU.
% compile_loopUpdateEVLIFNetGPU_fast(N,nSpikeGenCells,nCells2record)
%   INPUTS:
%       N              - # of simulated neurons in the network to be simulated
%
%       nSpikeGenCells - # of Poisson spike generator neurons in the
%                        network to be simulated
%
%       nCells2record  - # of cells whose spike times will be written to
%                        file
cfg = coder.gpuConfig('mex');
cfg.InitFltsAndDblsToZero = false;
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

N = net.nNeurons;
nSpikeGenCells = net.nSpikeGeneratorNeurons;

%% Define argument types for entry-point 'loopUpdateNetGPU'.
ARGS = cell(1,1);
ARGS{1} = cell(53,1);
ARGS{1}{1} = coder.typeof(single(0),[N   1],'Gpu',true); % V
ARGS{1}{2} = coder.typeof(single(0),[N   1],'Gpu',true); % Vreset
ARGS{1}{3} = coder.typeof(single(0),[N   1],'Gpu',true); % tau_ref
ARGS{1}{4} = coder.typeof(single(0),[N   1],'Gpu',true); % Vth
ARGS{1}{5} = coder.typeof(single(0),[N   1],'Gpu',true); % Vth0
ARGS{1}{6} = coder.typeof(single(0),[N   1],'Gpu',true); % Vth_max
ARGS{1}{7} = coder.typeof(single(0),[N   1],'Gpu',true); % Isra
ARGS{1}{8} = coder.typeof(single(0),[N   1],'Gpu',true); % tau_sra
ARGS{1}{9} = coder.typeof(single(0),[N   1],'Gpu',true); % a
ARGS{1}{10} = coder.typeof(single(0),[N   1],'Gpu',true); % b
ARGS{1}{11} = coder.typeof(single(0),[N   1],'Gpu',true); % VsynE
ARGS{1}{12} = coder.typeof(single(0),[N   1],'Gpu',true); % VsynI
ARGS{1}{13} = coder.typeof(single(0),[N   1],'Gpu',true); % GsynE
ARGS{1}{14} = coder.typeof(single(0),[N   1],'Gpu',true); % GsynI
ARGS{1}{15} = coder.typeof(single(0),[N N+nSpikeGenCells],'Gpu',true); % GsynMax
if (net.is_dynamic)
    ARGS{1}{16} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % tau_D
    ARGS{1}{17} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % tau_F
    ARGS{1}{18} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % f_fac
    ARGS{1}{19} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % D
    ARGS{1}{20} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % F
    ARGS{1}{21} = coder.typeof(false,[N+nSpikeGenCells 1],'Gpu',true); % has_facilitation
    ARGS{1}{22} = coder.typeof(false,[N+nSpikeGenCells 1],'Gpu',true); % has_depression
else
    ARGS{1}{16} = coder.typeof([]); % tau_D
    ARGS{1}{17} = coder.typeof([]); % tau_F
    ARGS{1}{18} = coder.typeof([]); % f_fac
    ARGS{1}{19} = coder.typeof([]); % D
    ARGS{1}{20} = coder.typeof([]); % F
    ARGS{1}{21} = coder.typeof([]); % has_facilitation
    ARGS{1}{22} = coder.typeof([]); % has_depression
end
ARGS{1}{23} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % p0
ARGS{1}{24} = coder.typeof(single(0),[N   1],'Gpu',true); % tau_synE
ARGS{1}{25} = coder.typeof(single(0),[N   1],'Gpu',true); % tau_synI
ARGS{1}{26} = coder.typeof(single(0),[N   1],'Gpu',true); % Cm
ARGS{1}{27} = coder.typeof(single(0),[N   1],'Gpu',true); % Gl
ARGS{1}{28} = coder.typeof(single(0),[N   1],'Gpu',true); % El
ARGS{1}{29} = coder.typeof(single(0),[N   1],'Gpu',true); % dth
ARGS{1}{30} = coder.typeof(single(0),[N   1],'Gpu',true); % Iapp
ARGS{1}{31} = coder.typeof(single(0),[N  1],'Gpu',true); % std_noise
ARGS{1}{32} = coder.typeof(single(0)); % dt
ARGS{1}{33} = coder.typeof(false,[N+nSpikeGenCells   1],'Gpu',true); % ecells
ARGS{1}{34} = coder.typeof(false,[N+nSpikeGenCells   1],'Gpu',true); % icells
if (nSpikeGenCells > 0)
    ARGS{1}{35} = coder.typeof(single(0),[nSpikeGenCells 1],'Gpu',true); % spikeGenProbs
else
    ARGS{1}{35} = coder.typeof(0,[0 0]);
end
if (nCells2record > 0)
    ARGS{1}{36} = coder.typeof(0,[nCells2record 1],'Gpu',true); % cells2record
else
    ARGS{1}{36} = coder.typeof(0,[0 0]);
end
if (net.is_plastic)
    ARGS{1}{37} = coder.typeof(false,[N N+nSpikeGenCells],'Gpu',true); % is_plastic logical matrix
    ARGS{1}{38} = coder.typeof('string',[1 20],[false true]); % plasticity_type string
    ARGS{1}{39} = coder.typeof(false,[N N+nSpikeGenCells],'Gpu',true); % Connectivity matrix
    ARGS{1}{40} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % r1
    ARGS{1}{41} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % r2
    ARGS{1}{42} = coder.typeof(single(0),[N 1],'Gpu',true); % o1
    ARGS{1}{43} = coder.typeof(single(0),[N 1],'Gpu',true); % o2
    ARGS{1}{44} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % A2plus
    ARGS{1}{45} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % A3plus
    ARGS{1}{46} = coder.typeof(single(0),[N 1],'Gpu',true); % A2minus
    ARGS{1}{47} = coder.typeof(single(0),[N 1],'Gpu',true); % A3minus
    ARGS{1}{48} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % tau_plus
    ARGS{1}{49} = coder.typeof(single(0),[N+nSpikeGenCells 1],'Gpu',true); % tau_x
    ARGS{1}{50} = coder.typeof(single(0),[N 1],'Gpu',true); % tau_minus
    ARGS{1}{51} = coder.typeof(single(0),[N 1],'Gpu',true); % tau_y
else
    ARGS{1}{37} = coder.typeof([]); % is_plastic logical matrix
    ARGS{1}{38} = coder.typeof('string',[1 20],[false true]); % plasticity_type string
    ARGS{1}{39} = coder.typeof([]); % Connectivity matrix
    ARGS{1}{40} = coder.typeof([]); % r1
    ARGS{1}{41} = coder.typeof([]); % r2
    ARGS{1}{42} = coder.typeof([]); % o1
    ARGS{1}{43} = coder.typeof([]); % o2
    ARGS{1}{44} = coder.typeof([]); % A2plus
    ARGS{1}{45} = coder.typeof([]); % A3plus
    ARGS{1}{46} = coder.typeof([]); % A2minus
    ARGS{1}{47} = coder.typeof([]); % A3minus
    ARGS{1}{48} = coder.typeof([]); % tau_plus
    ARGS{1}{49} = coder.typeof([]); % tau_x
    ARGS{1}{50} = coder.typeof([]); % tau_minus
    ARGS{1}{51} = coder.typeof([]); % tau_y
end
ARGS{1}{52} = coder.typeof(0); % nT (# of timesteps to simulate)
ARGS{1}{53} = coder.typeof(0); % file ID for spikes

%% Invoke MATLAB Coder.
codegen -config cfg runAEVLIFNetGPU -args ARGS{1}
end

