% noise_example
clear all;

net = AEVLIFnetwork();

net.addGroup('1',800,'excitatory',1,'std_noise',900e-12);
net.addGroup('2',200,'inhibitory',1,'std_noise',900e-12);

ntimesteps = 10000;
useGpu = true;
spikefile = 'spikes.bin';
sim_dir = '~/phd/easySim/results/noise_example';
[outputs] = easysim(net,ntimesteps,useGpu,...
                                    'sim_dir',sim_dir,...
                                    'spikefile',spikefile,...
                                    'recompile',true,...
                                    'recordVars',false,...
                                    'timeStepSize',.0001);

% retrieve spike data and compute firing rates
[spikeData] = readSpikes([sim_dir '/' spikefile],outputs.cells2record);
window = .1/outputs.dt; % 100ms
downsampleFactor = 10;
frs = getFiringRates(spikeData,length(outputs.cells2record),ntimesteps,outputs.dt,downsampleFactor,window);
