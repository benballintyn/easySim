function [V,Vth,Isra,GsynE,GsynI,D,F,r1,r2,o1,o2,Iapp] = resetVars(net,useGpu)
% Total # of neurons to be simulated
N = net.nNeurons; % total # of simulated neurons

% Go through each spike generator and assign it start and end indices
totalN = N;
offset = N;
for i=1:net.nSpikeGenerators
    start_ind = offset + 1;
    end_ind = start_ind + net.spikeGeneratorInfo(i).N - 1;
    net.spikeGeneratorInfo(i).start_ind = start_ind;
    net.spikeGeneratorInfo(i).end_ind = end_ind;
    offset = offset + net.spikeGeneratorInfo(i).N;
    totalN = totalN + net.spikeGeneratorInfo(i).N;
end
nSpikeGen = totalN - N; % total number of poisson spike generator neurons

if (useGpu)
    % Variables that may change with time
    V =         gpuArray(zeros(N,1,'single')); % membrane voltage
    Vth =       gpuArray(zeros(N,1,'single')); % Spike threshold
    GsynE =     gpuArray(zeros(N,1,'single')); % total excitatory synaptic conductance
    GsynI =     gpuArray(zeros(N,1,'single')); % total inhibitory synaptic conductance
    Iapp =      gpuArray(zeros(N,1,'single')); % Applied current
    
    if (isa(net,'AEVLIFnetwork'))
        Isra =       gpuArray(zeros(N,1,'single')); % Adaptation current
    else
        Isra = single(nan); % must use nan instead of [] to avoid compilation problems
    end
    
    if (net.is_dynamic)
        D     = gpuArray(ones(totalN,1,'single')); % Depression variable
        F     = gpuArray(ones(totalN,1,'single')); % Facilitation variable
    else
        D     = single(nan); % must use 0 instead of [] to avoid compilation problems
        F     = single(nan); % must use 0 instead of [] to avoid compilation problems
    end
    
    if (net.is_plastic)
        r1 =         gpuArray(zeros(totalN,1,'single')); % presynaptic plasticity variable 1
        r2 =         gpuArray(zeros(totalN,1,'single')); % presynaptic plasticity variable 2
        o1 =         gpuArray(zeros(N,1,'single')); % postsynaptic plasticity variable 1
        o2 =         gpuArray(zeros(N,1,'single'));
    else
        r1 =         single(nan); % must use 0 instead of [] to avoid compilation problems
        r2 =         single(nan); % must use 0 instead of [] to avoid compilation problems
        o1 =         single(nan); % must use 0 instead of [] to avoid compilation problems
        o2 =         single(nan); % must use 0 instead of [] to avoid compilation problems
    end
    
else
    
    % Variables that may change with time
    V =         zeros(N,1);
    Vth =       zeros(N,1);
    GsynE =     zeros(N,1);
    GsynI =     zeros(N,1);
    Iapp =      zeros(N,1);
    
    if (isa(net,'AEVLIFnetwork'))
        Isra =          zeros(N,1);
    else
        Isra = nan; % must use 0 instead of [] to avoid compilation problems
    end
    if (net.is_dynamic)
        D     = ones(totalN,1); % Depression variable
        F     = ones(totalN,1); % Facilitation variable
    else
        D     = nan; % must use 0 instead of [] to avoid compilation problems
        F     = nan; % must use 0 instead of [] to avoid compilation problems
    end
    if (net.is_plastic)
        r1 =         zeros(totalN,1);
        r2 =         zeros(totalN,1);
        o1 =         zeros(N,1);
        o2 =         zeros(N,1);
    else
        r1 =         nan; % must use 0 instead of [] to avoid compilation problems
        r2 =         nan; % must use 0 instead of [] to avoid compilation problems
        o1 =         nan; % must use 0 instead of [] to avoid compilation problems
        o2 =         nan; % must use 0 instead of [] to avoid compilation problems
    end
end

for i=1:net.nGroups
    preStart = net.groupInfo(i).start_ind;
    preEnd = net.groupInfo(i).end_ind;
    groupN = preEnd - preStart + 1;
    % initialize all nonzero variables for current group (V(0), Vth, VsynE, VsynI, dGref,
    % tau_ref, tau_synE, tau_synI, Cm, Gl, El, Ek, dth)
    V(preStart:preEnd)         = normrnd(net.groupInfo(i).mean_V0,net.groupInfo(i).std_V0,groupN,1);
    Vth(preStart:preEnd)       = normrnd(net.groupInfo(i).mean_Vth0,net.groupInfo(i).std_Vth0,groupN,1);
end

end

