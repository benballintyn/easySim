% compile updateNet.m for gpu
maxNeurons = 100000;

cfg = coder.gpuConfig('mex');
cfg.GpuConfig.CompilerFlags = '--fmad=false';
cfg.GenerateReport = true;

ARGS = cell(21,1);
ARGS{1} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{2} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{3} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{4} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{5} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{6} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{7} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{8} = coder.typeof(0,[maxNeurons maxNeurons],[1 1]);
ARGS{9} = coder.typeof(0,[maxNeurons maxNeurons],[1 1]);
ARGS{10} = coder.typeof(0,[maxNeurons maxNeurons],[1 1]);
ARGS{11} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{12} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{13} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{14} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{15} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{16} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{17} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{18} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{19} = coder.typeof(0,[1],[0]);
ARGS{20} = coder.typeof(0,[maxNeurons 1],[1 0]);
ARGS{21} = coder.typeof(0,[maxNeurons 1],[1 0]);

codegen updateNet -args ARGS -nargout 4


