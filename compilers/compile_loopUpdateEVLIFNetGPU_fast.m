function [] = compile_loopUpdateEVLIFNetGPU_fast(N)
cfg = coder.gpuConfig('mex');
cfg.InitFltsAndDblsToZero = false;
cfg.GenerateReport = true;
cfg.ReportPotentialDifferences = false;

%% Define argument types for entry-point 'loopUpdateNetGPU'.
ARGS = cell(1,1);
ARGS{1} = cell(24,1);
ARGS{1}{1} = coder.typeof(single(0),[N   1],'Gpu',true); % V
ARGS{1}{2} = coder.typeof(single(0),[N   1],'Gpu',true); % tau_ref
ARGS{1}{3} = coder.typeof(single(0),[N   1],'Gpu',true); % Vth
ARGS{1}{4} = coder.typeof(single(0),[N   1],'Gpu',true); % Vth0
ARGS{1}{5} = coder.typeof(single(0),[N   1],'Gpu',true); % Vth_max
ARGS{1}{6} = coder.typeof(single(0),[N   1],'Gpu',true); % VsynE
ARGS{1}{7} = coder.typeof(single(0),[N   1],'Gpu',true); % VsynI
ARGS{1}{8} = coder.typeof(single(0),[N   1],'Gpu',true); % GsynE
ARGS{1}{9} = coder.typeof(single(0),[N   1],'Gpu',true); % GsynI
ARGS{1}{10} = coder.typeof(single(0),[N   1],'Gpu',true); % maxGsynE
ARGS{1}{11} = coder.typeof(single(0),[N   1],'Gpu',true); % maxGsynI
ARGS{1}{12} = coder.typeof(single(0),[N N],'Gpu',true); % dGsyn
ARGS{1}{13} = coder.typeof(single(0),[N   1],'Gpu',true); % tau_synE
ARGS{1}{14} = coder.typeof(single(0),[N   1],'Gpu',true); % tau_synI
ARGS{1}{15} = coder.typeof(single(0),[N   1],'Gpu',true); % Cm
ARGS{1}{16} = coder.typeof(single(0),[N   1],'Gpu',true); % Gl
ARGS{1}{17} = coder.typeof(single(0),[N   1],'Gpu',true); % El
ARGS{1}{18} = coder.typeof(single(0),[N   1],'Gpu',true); % dth
ARGS{1}{19} = coder.typeof(single(0),[N   1],'Gpu',true); % Iapp
ARGS{1}{20} = coder.typeof(single(0)); % dt
ARGS{1}{21} = coder.typeof(false,[N   1],'Gpu',true); % ecells
ARGS{1}{22} = coder.typeof(false,[N   1],'Gpu',true); % icells
ARGS{1}{23} = coder.typeof(0); % nT (# of timesteps to simulate)
ARGS{1}{24} = coder.typeof(0); % file ID for spikes

%% Invoke MATLAB Coder.
codegen -config cfg loopUpdateEVLIFNetGPU_fast -args ARGS{1}
end

