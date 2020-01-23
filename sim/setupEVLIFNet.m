function [V,Vreset,tau_ref,Vth,Vth0,Vth_max,VsynE,VsynI,GsynE,GsynI,maxGsynE,maxGsynI,dGsyn,tau_synE,...
          tau_synI,Cm,Gl,El,dth,Iapp,std_noise,dt,ecells,icells,spikeGenProbs,cells2record] = ...
          setupEVLIFNet(net,useGpu)
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
    maxGsynE = gpuArray(zeros(N,1,'single')); % maximum total excitatory synaptic conductance
    maxGsynI = gpuArray(zeros(N,1,'single')); % maximum total inhibitory synaptic conductance
    Iapp =     gpuArray(zeros(N,1,'single')); % Applied current
    Vth =      gpuArray(zeros(N,1,'single')); % Spike threshold
    VsynE =    gpuArray(ones(N,1,'single')); % Excitatory synaptic reversal potential
    VsynI =    gpuArray(ones(N,1,'single')); % Inhibitory synaptic reversal potential
    
    % Variables that will not change with time
    Vth0 =      gpuArray(zeros(N,1,'single')); % baseline spike threshold
    Vth_max =   gpuArray(zeros(N,1,'single')); % maximum spike threshold
    tau_ref =   gpuArray(zeros(N,1,'single')); % refactory period time constant (for variable threshold)
    dGsyn =     gpuArray(zeros(N,totalN,'single')); % synaptic 'weight' matrix (change in synaptic conductance for presynaptic spike)
    tau_synE =  gpuArray(zeros(N,1,'single')); % excitatory synaptic time constant
    tau_synI =  gpuArray(zeros(N,1,'single')); % inhibitory synaptic time constant
    Cm  =       gpuArray(zeros(N,1,'single')); % membrane capacitance
    Gl =        gpuArray(zeros(N,1,'single')); % leak conductance
    El =        gpuArray(zeros(N,1,'single')); % leak reversal potential
    dth =       gpuArray(zeros(N,1,'single')); % spike generation voltage range
    std_noise = gpuArray(zeros(N,1,'single')); % standard deviation of the noise current
    
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
    maxGsynE = zeros(N,1);
    maxGsynI = zeros(N,1);
    Iapp =     zeros(N,1);
    Vth =      zeros(N,1);
    VsynE =    ones(N,1);
    VsynI =    ones(N,1);

    % Variables that will not change with time
    Vth0 =          zeros(N,1);
    Vth_max =       zeros(N,1);
    tau_ref =       zeros(N,1);
    dGsyn =         zeros(N,totalN);
    tau_synE =      zeros(N,1);
    tau_synI =      zeros(N,1);
    Cm  =           zeros(N,1);
    Gl =            zeros(N,1);
    El =            zeros(N,1);
    dth =           zeros(N,1);
    std_noise =     zeros(N,1);
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
    V(preStart:preEnd) = normrnd(net.groupInfo(i).mean_V0,net.groupInfo(i).std_V0,groupN,1);
    Vreset(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Vreset,net.groupInfo(i).std_Vreset,groupN,1);
    Vth(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Vth0,net.groupInfo(i).std_Vth0,groupN,1);
    Vth0(preStart:preEnd) = Vth(preStart:preEnd);
    Vth_max(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Vth_max,net.groupInfo(i).std_Vth_max,groupN,1);
    VsynE(preStart:preEnd) = normrnd(net.groupInfo(i).mean_VsynE,net.groupInfo(i).std_VsynE,groupN,1);
    VsynI(preStart:preEnd) = normrnd(net.groupInfo(i).mean_VsynI,net.groupInfo(i).std_VsynI,groupN,1);
    maxGsynE(preStart:preEnd) = normrnd(net.groupInfo(i).mean_max_GsynE,net.groupInfo(i).std_max_GsynE,groupN,1);
    maxGsynI(preStart:preEnd) = normrnd(net.groupInfo(i).mean_max_GsynI,net.groupInfo(i).std_max_GsynI,groupN,1);
    tau_ref(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_ref,net.groupInfo(i).std_tau_ref,groupN,1);
    tau_synE(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_synE,net.groupInfo(i).std_tau_synE,groupN,1);
    tau_synI(preStart:preEnd) = normrnd(net.groupInfo(i).mean_tau_synI,net.groupInfo(i).std_tau_synI,groupN,1);
    Cm(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Cm,net.groupInfo(i).std_Cm,groupN,1);
    Gl(preStart:preEnd) = normrnd(net.groupInfo(i).mean_Gl,net.groupInfo(i).std_Gl,groupN,1);
    El(preStart:preEnd) = normrnd(net.groupInfo(i).mean_El,net.groupInfo(i).std_El,groupN,1);
    dth(preStart:preEnd) = normrnd(net.groupInfo(i).mean_dth,net.groupInfo(i).std_dth,groupN,1);
    std_noise(preStart:preEnd) = net.groupInfo(i).std_noise;
    
    % Go through all of this groups targets and call the .genConn() method
    % of the appropriate connectionType to create the dGsyn matrix (weight
    % matrix)
    targets = net.groupInfo(i).targets;
    for j=1:length(targets)
        postStart = net.groupInfo(targets(j)).start_ind;
        postEnd = net.groupInfo(targets(j)).end_ind;
        curdGsyn = net.groupInfo(i).connections(j).genConn();
        dGsyn(postStart:postEnd,preStart:preEnd) = curdGsyn;
    end
end

% Go through all of the spike generator groups and add their connections to
% the dGsyn matrix
for i=1:net.nSpikeGenerators
    preStart = net.spikeGeneratorInfo(i).start_ind;
    preEnd = net.spikeGeneratorInfo(i).end_ind;
    targets = net.spikeGeneratorInfo(i).targets;
    for j=1:length(net.spikeGeneratorInfo(i).targets)
        postStart = net.groupInfo(targets(j)).start_ind;
        postEnd = net.groupInfo(targets(j)).end_ind;
        curdGsyn = net.spikeGeneratorInfo(i).connections(j).genConn();
        dGsyn(postStart:postEnd,preStart:preEnd) = curdGsyn;
    end
end

% bound variables to required ranges
minval= single(1e-40);
maxGsynE = max(minval,maxGsynE);
maxGsynI = max(minval,maxGsynI);
tau_ref = max(minval,tau_ref);
tau_synE = max(minval,tau_synE);
tau_synI = max(minval,tau_synI);
Cm = max(minval,Cm);
Gl = max(0,Gl);
dth = max(0,dth);
% auto-detect dt as 10x smaller than smallest time constant
dt = gather(10^(floor(log10(min(min(tau_ref),min(min(tau_synE),min(tau_synI)))/10))));
std_noise = max(0,std_noise/sqrt(dt));
dGsyn = max(0,dGsyn);

% set spike generator spike probabilities
for i=1:net.nSpikeGenerators
    s = net.spikeGeneratorInfo(i).start_ind - N;
    e = net.spikeGeneratorInfo(i).end_ind - N;
    spikeGenProbs(s:e) = net.spikeGeneratorInfo(i).firingRate*dt;
end
end