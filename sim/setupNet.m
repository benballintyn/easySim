function [V,Vreset,Cm,Gl,El,Vth,Vth0,Vth_max,tau_ref,dth,p0,GsynE,GsynI,VsynE,VsynI,tau_synE,tau_synI,...
          Iapp,std_noise,GsynMax,Isra,tau_sra,a,b,tau_D,tau_F,f_fac,D,F,has_facilitation,has_depression,...
          ecells,icells,spikeGenProbs,cells2record,r1,r2,o1,o2,A2plus,A3plus,A2minus,A3minus,...
          tau_plus,tau_x,tau_minus,tau_y,is_plastic,C,dt] = ...
          setupNet(net,useGpu)
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
    offset = offset + net.spikeGeneratorInfo(i).N;
    totalN = totalN + net.spikeGeneratorInfo(i).N;
end
nSpikeGen = totalN - N; % total number of poisson spike generator neurons

% If using a GPU, initialize all variables as gpuArrays of single precision
% otherwise if using the CPU use normal arrays with double precision
if (useGpu)
    % Variables that may change with time
    V =         gpuArray(zeros(N,1,'single')); % membrane voltage
    Vreset =    gpuArray(zeros(N,1,'single')); % reset membrane voltage after spike
    Cm  =       gpuArray(zeros(N,1,'single')); % membrane capacitance
    Gl =        gpuArray(zeros(N,1,'single')); % leak conductance
    El =        gpuArray(zeros(N,1,'single')); % leak reversal potential
    Vth =       gpuArray(zeros(N,1,'single')); % Spike threshold
    Vth0 =      gpuArray(zeros(N,1,'single')); % baseline spike threshold
    Vth_max =   gpuArray(zeros(N,1,'single')); % maximum spike threshold
    tau_ref =   gpuArray(zeros(N,1,'single')); % refactory period time constant (for variable threshold)
    dth =       gpuArray(zeros(N,1,'single')); % spike generation voltage range
    p0 =        gpuArray(zeros(totalN,1,'single')); % release probability
    GsynE =     gpuArray(zeros(N,1,'single')); % total excitatory synaptic conductance
    GsynI =     gpuArray(zeros(N,1,'single')); % total inhibitory synaptic conductance
    VsynE =     gpuArray(ones(N,1,'single')); % Excitatory synaptic reversal potential
    VsynI =     gpuArray(ones(N,1,'single')); % Inhibitory synaptic reversal potential
    tau_synE =  gpuArray(zeros(N,1,'single')); % excitatory synaptic time constant
    tau_synI =  gpuArray(zeros(N,1,'single')); % inhibitory synaptic time constant
    Iapp =      gpuArray(zeros(N,1,'single')); % Applied current
    std_noise = gpuArray(zeros(N,1,'single')); % standard deviation of the noise current
    GsynMax =   gpuArray(zeros(N,totalN,'single')); % synaptic 'weight' matrix (change in synaptic conductance for presynaptic spike)
    if (isa(net,'AEVLIFnetwork'))
        Isra =       gpuArray(zeros(N,1,'single')); % Adaptation current
        tau_sra =    gpuArray(zeros(N,1,'single')); % Adaptation time constant
        a =          gpuArray(zeros(N,1,'single')); % Adaptation conductance
        b =          gpuArray(zeros(N,1,'single')); % Adaptation current increment on spike
    else
        Isra = single(nan); % must use 0 instead of [] to avoid compilation problems
        tau_sra = single(nan);
        a = single(nan);
        b = single(nan);
    end
    if (net.is_dynamic)
        tau_D = gpuArray(zeros(totalN,1,'single')); % synaptic depression time constant
        tau_F = gpuArray(zeros(totalN,1,'single')); % synaptic facilitation time constant
        f_fac = gpuArray(zeros(totalN,1,'single')); % facilitation factor (strength)
        D     = gpuArray(ones(totalN,1,'single')); % Depression variable
        F     = gpuArray(ones(totalN,1,'single')); % Facilitation variable
        has_depression = gpuArray(false(totalN,1));
        has_facilitation = gpuArray(false(totalN,1));
    else
        tau_D = single(nan);
        tau_F = single(nan);
        f_fac = single(nan);
        D     = single(nan); % must use 0 instead of [] to avoid compilation problems
        F     = single(nan); % must use 0 instead of [] to avoid compilation problems
        has_depression = single(nan);
        has_facilitation = single(nan);
    end
    if (net.is_plastic)
        r1 =         gpuArray(zeros(totalN,1,'single')); % presynaptic plasticity variable 1
        r2 =         gpuArray(zeros(totalN,1,'single')); % presynaptic plasticity variable 2
        o1 =         gpuArray(zeros(N,1,'single')); % postsynaptic plasticity variable 1
        o2 =         gpuArray(zeros(N,1,'single'));
        A2plus =     gpuArray(zeros(totalN,1,'single')); % doublet STDP LTP factor
        A3plus =     gpuArray(zeros(totalN,1,'single')); % triplet STDP LTP factor
        A2minus =    gpuArray(zeros(N,1,'single')); % doublet STDP LTD factor
        A3minus =    gpuArray(zeros(N,1,'single')); % triplet STDP LTD factor
        tau_plus =   gpuArray(zeros(totalN,1,'single')); % r1 time constant
        tau_x =      gpuArray(zeros(totalN,1,'single')); % r2 time constant
        tau_minus =  gpuArray(zeros(N,1,'single')); % o1 time constant
        tau_y =      gpuArray(zeros(N,1,'single')); % o2 time constant
        is_plastic = gpuArray(zeros(N,totalN,'logical')); % matrix to keep track of which synapses are plastic
        C =          gpuArray(zeros(N,totalN,'logical')); % Connectivity matrix
    else
        r1 =         single(nan); % must use 0 instead of [] to avoid compilation problems
        r2 =         single(nan); % must use 0 instead of [] to avoid compilation problems
        o1 =         single(nan); % must use 0 instead of [] to avoid compilation problems
        o2 =         single(nan); % must use 0 instead of [] to avoid compilation problems
        A2plus =     single(nan); % doublet STDP LTP factor
        A3plus =     single(nan); % triplet STDP LTP factor
        A2minus =    single(nan); % doublet STDP LTD factor
        A3minus =    single(nan); % triplet STDP LTD factor
        tau_plus =   single(nan); % r1 time constant
        tau_x =      single(nan); % r2 time constant
        tau_minus =  single(nan); % o1 time constant
        tau_y =      single(nan); % o2 time constant
        is_plastic = single(nan);
        C =          single(nan);
    end
    
    ecells =    gpuArray(zeros(totalN,1)); % logical vector specifying which cells are excitatory
    icells =    gpuArray(zeros(totalN,1)); % logical vector specifying which cells are inhibitory
    if (nSpikeGen == 0)
        spikeGenProbs = [];
    else
        spikeGenProbs = gpuArray(zeros(nSpikeGen,1,'single')); % spike probability for spike generators
    end
else
    % Variables that may change with time
    V =         zeros(N,1);
    Vreset =    zeros(N,1);
    Cm  =       zeros(N,1);
    Gl =        zeros(N,1);
    El =        zeros(N,1);
    Vth =       zeros(N,1);
    Vth0 =      zeros(N,1);
    Vth_max =   zeros(N,1);
    tau_ref =   zeros(N,1);
    dth =       zeros(N,1);
    p0 =        zeros(totalN,1);
    GsynE =     zeros(N,1);
    GsynI =     zeros(N,1);
    VsynE =     ones(N,1);
    VsynI =     ones(N,1);
    tau_synE =  zeros(N,1);
    tau_synI =  zeros(N,1);
    Iapp =      zeros(N,1);
    std_noise = zeros(N,1);
    GsynMax =   zeros(N,totalN);
    if (isa(net,'AEVLIFnetwork'))
        Isra =          zeros(N,1);
        tau_sra =       zeros(N,1);
        a =             zeros(N,1);
        b =             zeros(N,1);
    else
        Isra = nan; % must use 0 instead of [] to avoid compilation problems
        tau_sra = nan; % must use 0 instead of [] to avoid compilation problems
        a = nan;
        b = nan;
    end
    if (net.is_dynamic)
        tau_D = zeros(totalN,1); % synaptic depression time constant
        tau_F = zeros(totalN,1); % synaptic facilitation time constant
        f_fac = zeros(totalN,1); % facilitation factor (strength)
        D     = ones(totalN,1); % Depression variable
        F     = ones(totalN,1); % Facilitation variable
        has_facilitation = false(totalN,1);
        has_depression = false(totalN,1);
    else
        tau_D = nan;
        tau_F = nan;
        f_fac = nan;
        D     = nan; % must use 0 instead of [] to avoid compilation problems
        F     = nan; % must use 0 instead of [] to avoid compilation problems
        has_facilitation = nan;
        has_depression = nan;
    end
    if (net.is_plastic)
        r1 =         zeros(totalN,1);
        r2 =         zeros(totalN,1);
        o1 =         zeros(N,1);
        o2 =         zeros(N,1);
        A2plus =     zeros(totalN,1); % doublet STDP LTP factor
        A3plus =     zeros(totalN,1); % triplet STDP LTP factor
        A2minus =    zeros(N,1); % doublet STDP LTD factor
        A3minus =    zeros(N,1); % triplet STDP LTD factor
        tau_plus =   zeros(totalN,1); % r1 time constant
        tau_x =      zeros(totalN,1); % r2 time constant
        tau_minus =  zeros(N,1); % o1 time constant
        tau_y =      zeros(N,1); % o2 time constant
        is_plastic = false(N,totalN);
        C =          false(N,totalN);
    else
        r1 =         nan; % must use 0 instead of [] to avoid compilation problems
        r2 =         nan; % must use 0 instead of [] to avoid compilation problems
        o1 =         nan; % must use 0 instead of [] to avoid compilation problems
        o2 =         nan; % must use 0 instead of [] to avoid compilation problems
        A2plus =     nan; % doublet STDP LTP factor
        A3plus =     nan; % triplet STDP LTP factor
        A2minus =    nan; % doublet STDP LTD factor
        A3minus =    nan; % triplet STDP LTD factor
        tau_plus =   nan; % r1 time constant
        tau_x =      nan; % r2 time constant
        tau_minus =  nan; % o1 time constant
        tau_y =      nan; % o2 time constant
        is_plastic = nan;
        C =          nan;
    end
    
    ecells =        zeros(totalN,1);
    icells =        zeros(totalN,1);
    spikeGenProbs = zeros(nSpikeGen,1);
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
    Cm(preStart:preEnd)        = normrnd(net.groupInfo(i).mean_Cm,net.groupInfo(i).std_Cm,groupN,1);
    Gl(preStart:preEnd)        = normrnd(net.groupInfo(i).mean_Gl,net.groupInfo(i).std_Gl,groupN,1);
    El(preStart:preEnd)        = normrnd(net.groupInfo(i).mean_El,net.groupInfo(i).std_El,groupN,1);
    Vth(preStart:preEnd)       = normrnd(net.groupInfo(i).mean_Vth0,net.groupInfo(i).std_Vth0,groupN,1);
    Vth0(preStart:preEnd)      = Vth(preStart:preEnd);
    Vth_max(preStart:preEnd)   = normrnd(net.groupInfo(i).mean_Vth_max,net.groupInfo(i).std_Vth_max,groupN,1);
    tau_ref(preStart:preEnd)   = normrnd(net.groupInfo(i).mean_tau_ref,net.groupInfo(i).std_tau_ref,groupN,1);
    dth(preStart:preEnd)       = normrnd(net.groupInfo(i).mean_dth,net.groupInfo(i).std_dth,groupN,1);
    p0(preStart:preEnd)        = normrnd(net.groupInfo(i).mean_p0,net.groupInfo(i).std_p0,groupN,1);
    VsynE(preStart:preEnd)     = normrnd(net.groupInfo(i).mean_VsynE,net.groupInfo(i).std_VsynE,groupN,1);
    VsynI(preStart:preEnd)     = normrnd(net.groupInfo(i).mean_VsynI,net.groupInfo(i).std_VsynI,groupN,1);
    tau_synE(preStart:preEnd)  = normrnd(net.groupInfo(i).mean_tau_synE,net.groupInfo(i).std_tau_synE,groupN,1);
    tau_synI(preStart:preEnd)  = normrnd(net.groupInfo(i).mean_tau_synI,net.groupInfo(i).std_tau_synI,groupN,1);
    std_noise(preStart:preEnd) = net.groupInfo(i).std_noise;
    if (isa(net,'AEVLIFnetwork'))
        tau_sra(preStart:preEnd)   = normrnd(net.groupInfo(i).mean_tau_sra,net.groupInfo(i).std_tau_sra,groupN,1);
        a(preStart:preEnd)         = normrnd(net.groupInfo(i).mean_a,net.groupInfo(i).std_a,groupN,1);
        b(preStart:preEnd)         = normrnd(net.groupInfo(i).mean_b,net.groupInfo(i).std_b,groupN,1);
    end
    if (net.is_dynamic)
        tau_D(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_D,net.groupInfo(i).std_tau_D,groupN,1);
        tau_F(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_F,net.groupInfo(i).std_tau_F,groupN,1);
        f_fac(preStart:preEnd) = normrnd(net.groupInfo(i).mean_f_fac,net.groupInfo(i).std_f_fac,groupN,1);
        if (net.groupInfo(i).facilitating_synapses)
            has_facilitation(preStart:preEnd) = true;
        end
        if (net.groupInfo(i).depressed_synapses)
            has_depression(preStart:preEnd) = true;
        end
    end
    if (net.is_plastic)
        A2plus(preStart:preEnd)    = normrnd(net.groupInfo(i).mean_A2plus,net.groupInfo(i).std_A2plus,groupN,1);
        A3plus(preStart:preEnd)    = normrnd(net.groupInfo(i).mean_A3plus,net.groupInfo(i).std_A3plus,groupN,1);
        A2minus(preStart:preEnd)   = normrnd(net.groupInfo(i).mean_A2minus,net.groupInfo(i).std_A2minus,groupN,1);
        A3minus(preStart:preEnd)   = normrnd(net.groupInfo(i).mean_A3minus,net.groupInfo(i).std_A3minus,groupN,1);
        tau_plus(preStart:preEnd)  = normrnd(net.groupInfo(i).mean_tau_plus,net.groupInfo(i).std_tau_plus,groupN,1);
        tau_x(preStart:preEnd)     = normrnd(net.groupInfo(i).mean_tau_x,net.groupInfo(i).std_tau_x,groupN,1);
        tau_minus(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_minus,net.groupInfo(i).std_tau_minus,groupN,1);
        tau_y(preStart:preEnd)     = normrnd(net.groupInfo(i).mean_tau_y,net.groupInfo(i).std_tau_y,groupN,1);
    end
    
    % Go through all of this groups targets and call the .genConn() method
    % of the appropriate connectionType to create the dGsyn matrix (weight
    % matrix)
    targets = net.groupInfo(i).targets;
    for j=1:length(targets)
        postStart = net.groupInfo(targets(j)).start_ind;
        postEnd = net.groupInfo(targets(j)).end_ind;
        curGsynMax = net.groupInfo(i).connections(j).genConn();
        GsynMax(postStart:postEnd,preStart:preEnd) = curGsynMax;
        if (net.is_plastic)
            if (net.groupInfo(i).connectionParams{j}.is_plastic)
                is_plastic(postStart:postEnd,preStart:preEnd) = true;
            end
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
    if (net.is_dynamic)
        tau_D(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_tau_D,net.spikeGeneratorInfo(i).std_tau_D,groupN,1);
        tau_F(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_tau_F,net.spikeGeneratorInfo(i).std_tau_F,groupN,1);
        f_fac(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_f_fac,net.spikeGeneratorInfo(i).std_f_fac,groupN,1);
        if (net.spikeGeneratorInfo(i).facilitating_synapses)
            has_facilitation(preStart:preEnd) = true;
        end
        if (net.spikeGeneratorInfo(i).depressed_synapses)
            has_depression(preStart:preEnd) = true;
        end
    end
    if (net.is_plastic)
        A2plus(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_A2plus,net.spikeGeneratorInfo(i).std_A2plus,groupN,1);
        A3plus(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_A3plus,net.spikeGeneratorInfo(i).std_A3plus,groupN,1);
        tau_plus(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_tau_plus,net.spikeGeneratorInfo(i).std_tau_plus,groupN,1);
        tau_x(preStart:preEnd) = normrnd(net.spikeGeneratorInfo(i).mean_tau_x,net.spikeGeneratorInfo(i).std_tau_x,groupN,1);
    end
    
    targets = net.spikeGeneratorInfo(i).targets;
    for j=1:length(net.spikeGeneratorInfo(i).targets)
        postStart = net.groupInfo(targets(j)).start_ind;
        postEnd = net.groupInfo(targets(j)).end_ind;
        curGsynMax = net.spikeGeneratorInfo(i).connections(j).genConn();
        GsynMax(postStart:postEnd,preStart:preEnd) = curGsynMax;
        if (net.is_plastic)
            if (net.spikeGeneratorInfo(i).connectionParams{j}.is_plastic)
                is_plastic(postStart:postEnd,preStart:preEnd) = true;
            end
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
if (net.is_plastic)
    C = (GsynMax > 0);
end

% set spike generator spike probabilities
for i=1:net.nSpikeGenerators
    s = net.spikeGeneratorInfo(i).start_ind - N;
    e = net.spikeGeneratorInfo(i).end_ind - N;
    spikeGenProbs(s:e) = net.spikeGeneratorInfo(i).firingRate*dt;
end
end