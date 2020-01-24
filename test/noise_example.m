% noise_example
clear all;

net = AEVLIFnetwork();

net.addGroup('1',800,'excitatory',1,'std_noise',.5e-12);

weightFunction = @(x) 1e-8*rand(size(x));
gaussConnProbFunction = @(x)(sqrt(2)/(.05*sqrt(pi)))*exp(-x.^2/(2*.05^2));
conn_1_1_params.weightFunction = weightFunction;
conn_1_1_params.connProbFunction = gaussConnProbFunction;
conn_1_1_params.useWrap = true;
net.connect(1,1,'gaussian',conn_1_1_params);

ntimesteps = 100000;
useGpu = true;
spikefile = 'spikes.bin';
sim_dir = 'results/noise_example';
[dt,cells2record,sim_dir] = easysim(net,ntimesteps,useGpu,'sim_dir',sim_dir,'spikefile',spikefile,'recompile',true);

% retrieve spike data and compute firing rates
[spikeData] = readSpikes([sim_dir '/' spikefile],cells2record);
window = 100; % 100ms
downsampleFactor = 10;
frs = getFiringRates(spikeData,length(cells2record),ntimesteps,dt,downsampleFactor,window);
