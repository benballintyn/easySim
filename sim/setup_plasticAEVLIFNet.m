function [V,Vreset,tau_ref,Vth,Vth0,Vth_max,Isra,tau_sra,a,b,VsynE,VsynI,GsynE,GsynI,GsynMax,p0,...
          tau_synE,tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record,...
          is_plastic,C,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,tau_plus,tau_x,tau_minus,tau_y] = ...
          setup_plasticAEVLIFNet(net,useGpu)
% This function initializes all of the relevant variables for simulation
% based based on whether a GPU will be used to do the simulation or not.
% The output of this function should be fed directly into one of the
% simulation functions (e.g. loopUpdateEVLIFNetGPU_fast)
% setupEVLIFNet(net,useGpu)
%   INPUTS:
%       net     - network object (e.g. EVLIFnetwork)
%
%       useGpu  - true or false. If true, variables will be initialized as
%                 gpuArrays with single precision. If false, variables will
%                 be normal arrays with double precision.

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
    offset = offset + N;
    totalN = totalN + net.spikeGeneratorInfo(i).N;
end
nSpikeGen = totalN - N; % total number of poisson spike generator neurons

% If using a GPU, initialize all variables as gpuArrays of single precision
% otherwise if using the CPU use normal arrays with double precision
if (useGpu)
    % Variables that may change with time
    V =        gpuArray(zeros(N,1,'single')); % membrane voltage
    Vreset =   gpuArray(zeros(N,1,'single')); % reset membrane voltage after spike
    GsynE =    gpuArray(zeros(N,1,'single')); % total excitatory synaptic conductance
    GsynI =    gpuArray(zeros(N,1,'single')); % total inhibitory synaptic conductance
    Iapp =     gpuArray(zeros(N,1,'single')); % Applied current
    Vth =      gpuArray(zeros(N,1,'single')); % Spike threshold
    VsynE =    gpuArray(ones(N,1,'single')); % Excitatory synaptic reversal potential
    VsynI =    gpuArray(ones(N,1,'single')); % Inhibitory synaptic reversal potential
    GsynMax =  gpuArray(zeros(N,totalN,'single')); % synaptic 'weight' matrix (change in synaptic conductance for presynaptic spike)
    r1 =       gpuArray(zeros(totalN,1,'single')); % presynaptic plasticity variable 1
    r2 =       gpuArray(zeros(totalN,1,'single')); % presynaptic plasticity variable 2
    o1 =       gpuArray(zeros(N,1,'single')); % postsynaptic plasticity variable 1
    o2 =       gpuArray(zeros(N,1,'single'));
    
    % Variables that will not change with time
    Vth0 =       gpuArray(zeros(N,1,'single')); % baseline spike threshold
    Vth_max =    gpuArray(zeros(N,1,'single')); % maximum spike threshold
    tau_ref =    gpuArray(zeros(N,1,'single')); % refactory period time constant (for variable threshold)
    Isra =       gpuArray(zeros(N,1,'single')); % Adaptation current
    tau_sra =    gpuArray(zeros(N,1,'single')); % Adaptation time constant
    a =          gpuArray(zeros(N,1,'single')); % Adaptation conductance
    b =          gpuArray(zeros(N,1,'single')); % Adaptation current increment on spike
    tau_synE =   gpuArray(zeros(N,1,'single')); % excitatory synaptic time constant
    tau_synI =   gpuArray(zeros(N,1,'single')); % inhibitory synaptic time constant
    Cm  =        gpuArray(zeros(N,1,'single')); % membrane capacitance
    Gl =         gpuArray(zeros(N,1,'single')); % leak conductance
    El =         gpuArray(zeros(N,1,'single')); % leak reversal potential
    dth =        gpuArray(zeros(N,1,'single')); % spike generation voltage range
    p0 =         gpuArray(zeros(totalN,1,'single')); % release probability
    std_noise =  gpuArray(zeros(N,1,'single')); % standard deviation of the noise current
    is_plastic = gpuArray(zeros(N,totalN,'logical')); % matrix to keep track of which synapses are plastic
    C =          gpuArray(zeros(N,totalN,'logical')); % Connectivity matrix
    A2plus =     gpuArray(zeros(totalN,1,'single')); % doublet STDP LTP factor
    A3plus =     gpuArray(zeros(totalN,1,'single')); % triplet STDP LTP factor
    A2minus =    gpuArray(zeros(N,1,'single')); % doublet STDP LTD factor
    A3minus =    gpuArray(zeros(N,1,'single')); % triplet STDP LTD factor
    tau_plus =   gpuArray(zeros(totalN,1,'single')); % r1 time constant
    tau_x =      gpuArray(zeros(totalN,1,'single')); % r2 time constant
    tau_minus =  gpuArray(zeros(N,1,'single')); % o1 time constant
    tau_y =      gpuArray(zeros(N,1,'single')); % o2 time constant
    
    ecells =    gpuArray(zeros(totalN,1)); % logical vector specifying which cells are excitatory
    icells =    gpuArray(zeros(totalN,1)); % logical vector specifying which cells are inhibitory
    if (nSpikeGen == 0)
        spikeGenProbs = [];
    else
        spikeGenProbs = gpuArray(zeros(nSpikeGen,1,'single')); % spike probability for spike generators
    end
else
    % Variables that may change with time
    V =        zeros(N,1);
    Vreset =   zeros(N,1);
    GsynE =    zeros(N,1);
    GsynI =    zeros(N,1);
    Iapp =     zeros(N,1);
    Vth =      zeros(N,1);
    VsynE =    ones(N,1);
    VsynI =    ones(N,1);
    GsynMax =  zeros(N,totalN);
    r1 = zeros(totalN,1);
    r2 = zeros(totalN,1);
    o1 = zeros(N,1);
    o2 = zeros(N,1);
    
    % Variables that will not change with time
    Vth0 =          zeros(N,1);
    Vth_max =       zeros(N,1);
    tau_ref =       zeros(N,1);
    Isra =          zeros(N,1);
    tau_sra =       zeros(N,1);
    a =             zeros(N,1);
    b =             zeros(N,1);
    tau_synE =      zeros(N,1);
    tau_synI =      zeros(N,1);
    Cm  =           zeros(N,1);
    Gl =            zeros(N,1);
    El =            zeros(N,1);
    dth =           zeros(N,1);
    p0 =            zeros(totalN,1);
    std_noise =     zeros(N,1);
    ecells =        zeros(totalN,1);
    icells =        zeros(totalN,1);
    spikeGenProbs = zeros(nSpikeGen,1);
    is_plastic = false(N,totalN);
    C =          false(N,totalN);
    A2plus =     zeros(totalN,1); % doublet STDP LTP factor
    A3plus =     zeros(totalN,1); % triplet STDP LTP factor
    A2minus =    zeros(N,1); % doublet STDP LTD factor
    A3minus =    zeros(N,1); % triplet STDP LTD factor
    tau_plus =   zeros(totalN,1); % r1 time constant
    tau_x =      zeros(totalN,1); % r2 time constant
    tau_minus =  zeros(N,1); % o1 time constant
    tau_y =      zeros(N,1); % o2 time constant
end

% Go through each group and add them to the ecells and icells vectors. Also
% determine which neurons' spikes need to be recorded
cells2record = [];
for i=1:net.nGroups
    s = net.groupInfo(i).start_ind;
    e = net.groupInfo(i).end_ind;
    if (net.groupInfo(i).isExcitatory)
        ecells(s:e) = 1;
    elseif (net.groupInfo(i).isInhibitory)
        icells(s:e) = 1;
    else
        error('Group is neither excitatory or inhibitory?')
    end
    if (net.groupInfo(i).record)
        cells2record = [cells2record s:e];
    end
end
if (useGpu)
    if (length(cells2record) > 0)
        cells2record = gpuArray(cells2record');
    end
else
    cells2record = cells2record';
end

% Go through each of the spike generator groups and add them to the ecells
% and icells vectors
for i=1:net.nSpikeGenerators
    s = net.spikeGeneratorInfo(i).start_ind;
    e = net.spikeGeneratorInfo(i).end_ind;
    if (net.spikeGeneratorInfo(i).isExcitatory)
        ecells(s:e) = 1;
    elseif (net.spikeGeneratorInfo(i).isInhibitory)
        icells(s:e) = 1;
    else
        error('spikeGenerator is neighter excitatory or inhibitory')
    end
    %spikeGenProbs((s-N):(e-N)) = net.spikeGeneratorInfo(i).firing_rate*simobj.dt;
end
ecells = logical(ecells);
icells = logical(icells);

for i=1:net.nGroups
    preStart = net.groupInfo(i).start_ind;
    preEnd = net.groupInfo(i).end_ind;
    groupN = preEnd - preStart + 1;
    % initialize all nonzero variables for current group (V(0), Vth, VsynE, VsynI, dGref,
    % tau_ref, tau_synE, tau_synI, Cm, Gl, El, Ek, dth)
    V(preStart:preEnd)         = normrnd(net.groupInfo(i).mean_V0,net.groupInfo(i).std_V0,groupN,1);
    Vreset(preStart:preEnd)    = normrnd(net.groupInfo(i).mean_Vreset,net.groupInfo(i).std_Vreset,groupN,1);
    Vth(preStart:preEnd)       = normrnd(net.groupInfo(i).mean_Vth0,net.groupInfo(i).std_Vth0,groupN,1);
    Vth0(preStart:preEnd)      = Vth(preStart:preEnd);
    Vth_max(preStart:preEnd)   = normrnd(net.groupInfo(i).mean_Vth_max,net.groupInfo(i).std_Vth_max,groupN,1);
    VsynE(preStart:preEnd)     = normrnd(net.groupInfo(i).mean_VsynE,net.groupInfo(i).std_VsynE,groupN,1);
    VsynI(preStart:preEnd)     = normrnd(net.groupInfo(i).mean_VsynI,net.groupInfo(i).std_VsynI,groupN,1);
    tau_ref(preStart:preEnd)   = normrnd(net.groupInfo(i).mean_tau_ref,net.groupInfo(i).std_tau_ref,groupN,1);
    tau_sra(preStart:preEnd)   = normrnd(net.groupInfo(i).mean_tau_sra,net.groupInfo(i).std_tau_sra,groupN,1);
    a(preStart:preEnd)         = normrnd(net.groupInfo(i).mean_a,net.groupInfo(i).std_a,groupN,1);
    b(preStart:preEnd)         = normrnd(net.groupInfo(i).mean_b,net.groupInfo(i).std_b,groupN,1);
    tau_synE(preStart:preEnd)  = normrnd(net.groupInfo(i).mean_tau_synE,net.groupInfo(i).std_tau_synE,groupN,1);
    tau_synI(preStart:preEnd)  = normrnd(net.groupInfo(i).mean_tau_synI,net.groupInfo(i).std_tau_synI,groupN,1);
    Cm(preStart:preEnd)        = normrnd(net.groupInfo(i).mean_Cm,net.groupInfo(i).std_Cm,groupN,1);
    Gl(preStart:preEnd)        = normrnd(net.groupInfo(i).mean_Gl,net.groupInfo(i).std_Gl,groupN,1);
    El(preStart:preEnd)        = normrnd(net.groupInfo(i).mean_El,net.groupInfo(i).std_El,groupN,1);
    dth(preStart:preEnd)       = normrnd(net.groupInfo(i).mean_dth,net.groupInfo(i).std_dth,groupN,1);
    p0(preStart:preEnd)        = normrnd(net.groupInfo(i).mean_p0,net.groupInfo(i).std_p0,groupN,1);
    A2plus(preStart:preEnd)    = normrnd(net.groupInfo(i).mean_A2plus,net.groupInfo(i).std_A2plus,groupN,1);
    A3plus(preStart:preEnd)    = normrnd(net.groupInfo(i).mean_A3plus,net.groupInfo(i).std_A3plus,groupN,1);
    A2minus(preStart:preEnd)   = normrnd(net.groupInfo(i).mean_A2minus,net.groupInfo(i).std_A2minus,groupN,1);
    A3minus(preStart:preEnd)   = normrnd(net.groupInfo(i).mean_A3minus,net.groupInfo(i).std_A3minus,groupN,1);
    tau_plus(preStart:preEnd)  = normrnd(net.groupInfo(i).mean_tau_plus,net.groupInfo(i).std_tau_plus,groupN,1);
    tau_x(preStart:preEnd)     = normrnd(net.groupInfo(i).mean_tau_x,net.groupInfo(i).std_tau_x,groupN,1);
    tau_minus(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_minus,net.groupInfo(i).std_tau_minus,groupN,1);
    tau_y(preStart:preEnd)     = normrnd(net.groupInfo(i).mean_tau_y,net.groupInfo(i).std_tau_y,groupN,1);
    std_noise(preStart:preEnd) = net.groupInfo(i).std_noise;
    
    % Go through all of this groups targets and call the .genConn() method
    % of the appropriate connectionType to create the dGsyn matrix (weight
    % matrix)
    targets = net.groupInfo(i).targets;
    for j=1:length(targets)
        postStart = net.groupInfo(targets(j)).start_ind;
        postEnd = net.groupInfo(targets(j)).end_ind;
        curGsynMax = net.groupInfo(i).connections(j).genConn();
        GsynMax(postStart:postEnd,preStart:preEnd) = curGsynMax;
        if (net.groupInfo(i).connectionParams{j}.is_plastic)
            is_plastic(postStart:postEnd,preStart:preEnd) = true;
        end
    end
end

% Go through all of the spike generator groups and initialize plasticity variables
% and add their connections to the dGsyn matrix
for i=1:net.nSpikeGenerators
    preStart = net.spikeGeneratorInfo(i).start_ind;
    preEnd = net.spikeGeneratorInfo(i).end_ind;
    groupN = preEnd - preStart + 1;
    p0(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_p0,net.spikeGeneratorInfo(i).std_p0,groupN,1);
    A2plus(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_A2plus,net.spikeGeneratorInfo(i).std_A2plus,groupN,1);
    A3plus(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_A3plus,net.spikeGeneratorInfo(i).std_A3plus,groupN,1);
    tau_plus(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_tau_plus,net.spikeGeneratorInfo(i).std_tau_plus,groupN,1);
    tau_x(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_tau_x,net.spikeGeneratorInfo(i).std_tau_x,groupN,1);
    
    targets = net.spikeGeneratorInfo(i).targets;
    for j=1:length(net.spikeGeneratorInfo(i).targets)
        postStart = net.groupInfo(targets(j)).start_ind;
        postEnd = net.groupInfo(targets(j)).end_ind;
        curGsynMax = net.spikeGeneratorInfo(i).connections(j).genConn();
        GsynMax(postStart:postEnd,preStart:preEnd) = curGsynMax;
        if (net.spikeGeneratorInfo(i).connectionParams{j}.is_plastic)
            is_plastic(postStart:postEnd,preStart:preEnd) = true;
        end
    end
end

% bound variables to required ranges
if (useGpu)
    minval= single(1e-40);
else
    minval=1e-100;
end
tau_ref = max(minval,tau_ref);
tau_synE = max(minval,tau_synE);
tau_synI = max(minval,tau_synI);
Cm = max(minval,Cm);
Gl = max(0,Gl);
dth = max(0,dth);
p0 = max(0,p0);
% auto-detect dt as 10x smaller than smallest time constant
dt = gather(10^(floor(log10(min(min(tau_ref),min(min(tau_synE),min(tau_synI)))/10))));
std_noise = max(0,std_noise/sqrt(dt));
GsynMax = max(0,GsynMax);
C = (GsynMax > 0);

% set spike generator spike probabilities
for i=1:net.nSpikeGenerators
    s = net.spikeGeneratorInfo(i).start_ind - N;
    e = net.spikeGeneratorInfo(i).end_ind - N;
    spikeGenProbs(s:e) = net.spikeGeneratorInfo(i).firingRate*dt;
end
end