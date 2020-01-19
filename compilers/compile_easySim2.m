function [] = compile_easySim2(N)
cfg = coder.gpuConfig('mex');
cfg.GpuConfig.CompilerFlags = '--fmad=false';
cfg.GenerateReport = true;

x = gpuArray(zeros(N,1,'single'));
x2 = gpuArray(zeros(N,N,'single'));
x3 = gpuArray(logical(zeros(1,N)));
x4=0;
codegen updateNet -args {x,x,x,x,x,x,x,x,x,x,x,x2,x,x,x,x,x,x,x,x,x4,x3,x3} -nargout 5

end

