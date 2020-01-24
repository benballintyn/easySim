function [] = compile_loopUpdateAEVLIFNetGPU_fast(N,nSpikeGenCells,nCells2record)
% This function compiles the loopUpdateAEVLIFNetGPU_fast function into CUDA
% code and produces a mex file loopUpdateAEVLIFNetGPU_fast_mex to run a
% simulation on the GPU.
% compile_loopUpdateAEVLIFNetGPU_fast(N,nSpikeGenCells,nCells2record)
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

%% Define argument types for entry-point 'loopUpdateAEVLIFNetGPU_fast'.
ARGS = cell(1,1);
ARGS{1} = cell(32,1);
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
ARGS{1}{15} = coder.typeof(single(0),[N   1],'Gpu',true); % maxGsynE
ARGS{1}{16} = coder.typeof(single(0),[N   1],'Gpu',true); % maxGsynI
ARGS{1}{17} = coder.typeof(single(0),[N N],'Gpu',true); % dGsyn
ARGS{1}{18} = coder.typeof(single(0),[N   1],'Gpu',true); % tau_synE
ARGS{1}{19} = coder.typeof(single(0),[N   1],'Gpu',true); % tau_synI
ARGS{1}{20} = coder.typeof(single(0),[N   1],'Gpu',true); % Cm
ARGS{1}{21} = coder.typeof(single(0),[N   1],'Gpu',true); % Gl
ARGS{1}{22} = coder.typeof(single(0),[N   1],'Gpu',true); % El
ARGS{1}{23} = coder.typeof(single(0),[N   1],'Gpu',true); % dth
ARGS{1}{24} = coder.typeof(single(0),[N   1],'Gpu',true); % Iapp
ARGS{1}{25} = coder.typeof(single(0),[N  1],'Gpu',true); % std_noise
ARGS{1}{26} = coder.typeof(single(0)); % dt
ARGS{1}{27} = coder.typeof(false,[N   1],'Gpu',true); % ecells
ARGS{1}{28} = coder.typeof(false,[N   1],'Gpu',true); % icells
if (nSpikeGenCells > 0)
    ARGS{1}{29} = coder.typeof(single(0),[nSpikeGenCells 1],'Gpu',true); % spikeGenProbs
else
    ARGS{1}{29} = coder.typeof(0,[0 0]);
end
if (nCells2record > 0)
    ARGS{1}{30} = coder.typeof(0,[nCells2record 1],'Gpu',true); % cells2record
else
    ARGS{1}{30} = coder.typeof(0,[0 0]);
end
ARGS{1}{31} = coder.typeof(0); % nT (# of timesteps to simulate)
ARGS{1}{32} = coder.typeof(0); % file ID for spikes

%% Invoke MATLAB Coder.
codegen -config cfg loopUpdateAEVLIFNetGPU_fast -args ARGS{1}
end

