% noise_example
clear all;

net = EVLIFnetwork();

net.addGroup('1',1000,'excitatory',1,'std_noise',1000e-12);

ntimesteps = 100000;
useGpu = false;
spikefile = 'spikes.bin';
sim_dir = 'results/noise_example';
[outputs] = easysim(net,ntimesteps,useGpu,...
                                    'sim_dir',sim_dir,...
                                    'spikefile',spikefile,...
                                    'recompile',true,...
                                    'recordVars',true,...
                                    'timeStepSize',.0001);

% retrieve spike data and compute firing rates
[spikeData] = readSpikes([sim_dir '/' spikefile],outputs.cells2record);
window = .1/outputs.dt; % 100ms
downsampleFactor = 10;
frs = getFiringRates(spikeData,length(outputs.cells2record),ntimesteps,outputs.dt,downsampleFactor,window);
