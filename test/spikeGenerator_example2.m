% spikeGenerator_example
clear all;
% Initialize the network
net = AEVLIFnetwork();

% Add a single group of 1000 excitatory neurons
net.addGroup('1',1000,'excitatory',1);

% Add a group of 1000 excitatory Poisson spikeGenerator neurons with firing
% rate 10Hz
net.addSpikeGenerator('spkE',1000,'excitatory',10);

% Add a group of 1000 inhibitory Poisson spikeGenerator neurons with firing
% rate 10Hz
net.addSpikeGenerator('spkI',1000,'inhibitory',10);

% Connect the spikeGenerator groups to the excitatory group with a random
% connection (P(connect) = .1, synaptic weights in [0 3e-8])
maxWeight = 3e-8;
weightRange = 0:1e-12:maxWeight;
uniformWeightDist = weightDistribution(weightRange,ones(1,length(weightRange))*(1/length(weightRange)));
conn_params.connProb = .1;
conn_params.weightDistribution = uniformWeightDist;
net.connect(-1,1,'random',conn_params)
net.connect(-2,1,'random',conn_params)

% Run the network for 10000 time steps on the GPU or CPU
% Save the spikes from group 1 to results/spikeGenerator_example/spikes.bin
useGpu = false;
ntimesteps = 10000;
spikefile = 'spikes.bin';
sim_dir = 'results/spikeGenerator_example';
if (useGpu)
    [outputs] = easysim(net,ntimesteps,useGpu,'sim_dir',sim_dir,'spikefile',spikefile,'recompile',true);
else
    [outputs] = easysim(net,ntimesteps,useGpu,'sim_dir',sim_dir,...
                                              'spikefile',spikefile,...
                                              'recompile',true,...
                                              'recordVars',true);
end

% Extract spikes from the data file and convert them to firing rate traces
[spikeData] = readSpikes([sim_dir '/' spikefile],outputs.cells2record);
window = .1/outputs.dt; % 100ms
downsampleFactor = 10;
frs = getFiringRates(spikeData,length(outputs.cells2record),ntimesteps,outputs.dt,downsampleFactor,window);

tvec = linspace(0,ntimesteps*outputs.dt,ntimesteps/downsampleFactor);
plot(tvec,frs)

