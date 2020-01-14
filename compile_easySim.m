function [] = compile_easySim(N)
cfg = coder.gpuConfig('mex');
cfg.GpuConfig.CompilerFlags = '--fmad=false';
cfg.GenerateReport = true;

ARGS = cell(21,1);
ARGS{1} = coder.typeof(gpuArray(0),[N 1],[0 0]); % V (membrane Voltage)
ARGS{2} = coder.typeof(gpuArray(0),[N 1],[0 0]); % Gref (refractory conductance)
ARGS{3} = coder.typeof(gpuArray(0),[N 1],[0 0]); % dGref (refractory conductance change on spike)
ARGS{4} = coder.typeof(gpuArray(0),[N 1],[0 0]); % tau_ref (refractory time_constant)
ARGS{5} = coder.typeof(gpuArray(0),[N 1],[0 0]); % Vth (spike threshold)
ARGS{6} = coder.typeof(gpuArray(0),[N 1],[0 0]); % VsynE (excitatory synaptic reversal potential)
ARGS{7} = coder.typeof(gpuArray(0),[N 1],[0 0]); % VsynI (inhibitory synaptic reversal potential)
ARGS{8} = coder.typeof(gpuArray(0),[N 1],[0 0]); % GsynE (total excitatory synaptic conductance)
ARGS{9} = coder.typeof(gpuArray(0),[N 1],[0 0]); % GsynI (total inhibitory synaptic conductance)
ARGS{10} = coder.typeof(gpuArray(0),[N N],[0 0]); % dGsyn (synaptic strength matrix)
ARGS{11} = coder.typeof(gpuArray(0),[N 1],[0 0]); % tau_synE (excitatory synaptic decay time constant)
ARGS{12} = coder.typeof(gpuArray(0),[N 1],[0 0]); % tau_synI (inhibitory synaptic decay time constant)
ARGS{13} = coder.typeof(gpuArray(0),[N 1],[0 0]); % Cm (membrane capacitance)
ARGS{14} = coder.typeof(gpuArray(0),[N 1],[0 0]); % Gl (leak conductance)
ARGS{15} = coder.typeof(gpuArray(0),[N 1],[0 0]); % El (leak reversal potential)
ARGS{16} = coder.typeof(gpuArray(0),[N 1],[0 0]); % Ek
ARGS{17} = coder.typeof(gpuArray(0),[N 1],[0 0]); % dth
ARGS{18} = coder.typeof(gpuArray(0),[N 1],[0 0]); % Iapp
ARGS{19} = coder.typeof(0,[1],[0]); % dt
ARGS{20} = coder.typeof(gpuArray(false),[N 1],[0 0]); % ecells
ARGS{21} = coder.typeof(gpuArray(false),[N 1],[0 0]); % icells

codegen updateNet -args ARGS -nargout 5

end

